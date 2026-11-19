"""Opt-in, per-run block-task CIFTI GLM with explicit acquisition timing metadata.

Recipes live in config/task_models; run timing metadata lives beside the raw CSV.
Uses Nilearn for HRF design, AR(1) fitting and contrasts. No temporal censoring,
extra smoothing, MGTR substitution, or silently assumed scanner alignment.
"""
import argparse
import hashlib
import html
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import warnings
import xml.etree.ElementTree as ET

import nibabel as nib
from nibabel.processing import resample_from_to
import nilearn
from nilearn.glm.first_level import make_first_level_design_matrix, run_glm
from nilearn.glm.contrasts import compute_contrast
import numpy as np
import pandas as pd
from scipy.ndimage import binary_erosion
from scipy.stats import t as student_t
from statsmodels.stats.multitest import multipletests
import yaml

ROOT = Path(__file__).resolve().parents[1]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_recipe(path):
    r = yaml.safe_load(Path(path).read_text())
    require(isinstance(r, dict) and r.get('schema_version') == 1, 'Expected recipe schema_version: 1')
    require(r.get('adapter') == 'block_events', 'Only block_events adapter is supported')
    for key in ('task', 'file_prefix', 'recipe_id'):
        require(re.fullmatch(r'[A-Za-z0-9_-]+', r.get(key, '')) is not None, 'Invalid recipe '+key)
    require(r.get('baseline') == 'implicit_rest' and r.get('hrf') == 'spm + derivative'
            and r.get('noise_model') == 'ar1' and r.get('motion_expansion') == 24,
            'Unsupported baseline/HRF/noise/motion settings')
    require(r.get('additional_smoothing_mm') == 0 and r.get('global_signal') is False,
            'This backend supports no extra smoothing and no global signal regression')
    require(r.get('input_tags') == ['OCME', 'OCME+MEICA'], 'Specify paired OCME and OCME+MEICA inputs')
    require(r.get('tissue_means') == ['white_matter', 'csf'], 'WM and CSF means are required')
    require(len(r['conditions']) == 2 and len(set(r['conditions'])) == 2, 'Two distinct conditions required')
    for name, weights in r['contrasts'].items():
        require(re.fullmatch(r'[A-Za-z0-9_-]+', name) is not None, 'Unsafe contrast name')
        require(set(weights) == set(r['conditions']) and any(float(v) != 0 for v in weights.values()),
                'Each contrast must explicitly weight both conditions')
    return r


def discover(subject, registry, tasks, start_session):
    require(isinstance(tasks, list) and tasks and len(set(tasks)) == len(tasks), 'Explicit, nonempty, unique task list required')
    known = {}
    for task in tasks:
        require(isinstance(task, str) and re.fullmatch(r'[A-Za-z0-9_-]+', task) is not None, 'Invalid task name')
        path = Path(registry)/(task+'.yaml')
        require(path.is_file(), 'No supported task-model recipe for '+task+': '+str(path))
        recipe = read_recipe(path)
        require(recipe['task'] == task, 'Recipe task must match filename: '+str(path))
        known[task] = (path, recipe)
    jobs = []
    for task in tasks:
        require(task in known, 'No supported task-model recipe for '+task)
        path, r = known[task]
        found = 0
        for raw in sorted((subject/'func/unprocessed'/task).glob('session_*/run_*')):
            m = re.fullmatch(r'session_(\d+)/run_(\d+)', str(raw.relative_to(subject/'func/unprocessed'/task)))
            if not m or int(m[1]) < start_session:
                continue
            session, run = map(int, m.groups())
            events = raw/r['events_pattern'].format(session=session, run=run)
            timing = raw/r['timing_pattern'].format(session=session, run=run)
            require(events.is_file(), 'Missing events: '+str(events))
            require(timing.is_file(), 'Missing timing metadata: '+str(timing))
            jobs.append((path, r, raw, events, timing, session, run)); found += 1
        require(found, 'No selected raw runs for '+task)
    return jobs


def load_events(path, recipe, timing, n, tr, removed=0):
    require(timing.get('alignment_status') in ('confirmed', 'assumed_authorized'),
            'Scanner-to-BOLD alignment is pending: '+str(path))
    offset = timing.get('onset_offset_seconds')
    require(isinstance(offset, (int, float)) and np.isfinite(offset), 'Explicit finite timing offset required')
    require(timing.get('slice_reference') == 'pipeline_slice_timing_provenance', 'Unsupported/unverified slice reference')
    df = pd.read_csv(path)
    require(set(['onset', 'duration', 'trial_type', 'phase', 'hand']) <= set(df), 'Missing event columns')
    require(df.hand.eq(recipe['hand']).all(), 'Event hand differs from recipe')
    df.onset = pd.to_numeric(df.onset) + offset - removed*tr
    df.duration = pd.to_numeric(df.duration)
    require(np.isfinite(df[['onset', 'duration']].values).all() and df.duration.gt(0).all(), 'Invalid event timing')
    require(df.trial_type.isin(recipe['conditions']+['rest']).all(), 'Unknown condition labels')
    require(np.all(np.diff(df.onset) >= 0), 'Events must be ordered')
    require(np.all(df.onset.values[1:] >= (df.onset+df.duration).values[:-1]-2e-6), 'Events overlap')
    require((df.onset+df.duration).max() <= n*tr+0.01, 'Events extend beyond BOLD acquisition')
    events = df[df.phase.eq(recipe['phase']) & df.trial_type.isin(recipe['conditions'])].copy()
    for condition in recipe['conditions']:
        require(events.trial_type.eq(condition).any(), 'Missing condition '+condition)
    require(events.onset.min() >= -24, 'Task events too early after discarded volumes')
    return events[['onset', 'duration', 'trial_type']], df


def build_design(n, tr, events, motion, tissues, recipe, reference_seconds):
    require(motion.shape == (n, 6) and np.isfinite(motion).all(), 'Motion must be finite n_frames x 6 FSL MCF.par')
    derivative = np.vstack([np.zeros((1, 6)), np.diff(motion, axis=0)])
    fd = np.abs(derivative[:, 3:]).sum(1) + recipe['rotation_radius_mm']*np.abs(derivative[:, :3]).sum(1)
    spikes = np.flatnonzero(fd > recipe['fd_threshold_mm'])
    nuisance = np.column_stack([motion, derivative, motion**2, derivative**2, np.eye(n)[:, spikes], tissues])
    names = ['motion_'+str(i) for i in range(24)] + ['spike_'+str(i) for i in spikes] + recipe['tissue_means']
    # Nilearn can regularize a singular design; reject warnings rather than silently accepting altered task columns.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        x = make_first_level_design_matrix(np.arange(n)*tr+reference_seconds, events, hrf_model=recipe['hrf'],
                drift_model='cosine', high_pass=1/recipe['high_pass_seconds'], oversampling=16)
    require(not any('singular' in str(w.message).lower() for w in caught), 'Singular task/drift design')
    dropped = []
    # Remove constant/redundant nuisance columns, never task columns. Raw motion and FD are retained in diagnostics.
    for name, values in zip(names, nuisance.T):
        scale = values.std()
        if not np.isfinite(scale) or scale < 1e-10:
            dropped.append(name); continue
        values = (values-values.mean())/scale
        residual = values-x.values@np.linalg.lstsq(x.values, values, rcond=None)[0]
        if np.linalg.norm(residual) < 1e-7*np.linalg.norm(values):
            dropped.append(name); continue
        x[name] = values
    rank = int(np.linalg.matrix_rank(x.values))
    require(rank == x.shape[1] and n-rank >= recipe['minimum_residual_df'], 'Insufficient rank/residual degrees of freedom')
    require(np.linalg.cond(x.values) < 1e10, 'Ill-conditioned design')
    return x, fd, {'spike_frames': spikes.tolist(), 'dropped_redundant_nuisance': dropped,
                   'rank': rank, 'residual_df': n-rank, 'condition_number': float(np.linalg.cond(x.values))}


def slice_reference(raw_json, shift_file, tr):
    source = json.loads(Path(raw_json).read_text())
    times = np.asarray(source['SliceTiming'], dtype=float)
    shifts = np.loadtxt(shift_file, ndmin=1)
    require(times.shape == shifts.shape and np.isfinite(times).all() and np.isfinite(shifts).all(),
            'SliceTiming metadata and preprocessing shifts differ')
    unique = sorted(set(times.tolist())); middle = int(np.floor(len(unique)/2+0.5))
    expected = np.asarray([(middle-unique.index(t)-1)/len(unique) for t in times])
    require(np.allclose(shifts, expected, atol=1e-7), 'Unsupported slice-timing shifts; require pipeline rank-group convention')
    effective = times+shifts*tr
    require(np.ptp(effective) < 0.05, 'Effective slice references disagree by >50 ms')
    return float(np.median(effective)), {'effective_slice_reference_min_seconds': float(effective.min()),
                                        'effective_slice_reference_max_seconds': float(effective.max())}


def tissue_signals(volume, segmentation, recipe):
    seg = resample_from_to(segmentation, (volume.shape[:3], volume.affine), order=0)
    labels = np.asarray(seg.dataobj)
    signals, counts = [], {}
    masks = []
    for name, ids in [('white_matter', [2, 41]), ('csf', [4, 43, 14, 15])]:
        mask = np.isin(labels, ids)
        iterations = int(recipe['tissue_erosion_voxels'])
        if iterations:
            mask = binary_erosion(mask, iterations=iterations)
        require(mask.sum() >= recipe['minimum_tissue_voxels'], 'Too few voxels for '+name)
        masks.append(mask); counts[name] = int(mask.sum())
    # Read one frame at a time to avoid loading a second full 4D dataset.
    for i in range(volume.shape[3]):
        frame = np.asarray(volume.dataobj[..., i], dtype=np.float32)
        signals.append([float(frame[m].mean()) for m in masks])
    signals = np.asarray(signals)
    require(np.isfinite(signals).all(), 'Nonfinite tissue signals')
    return signals, counts


def write_scalar(data, names, axis, path):
    nib.save(nib.Cifti2Image(np.asarray(data, dtype=np.float32),
             nib.Cifti2Header.from_axes((nib.cifti2.ScalarAxis(names), axis))), str(path))


def receipt_ok(folder, fingerprint):
    receipt = folder/'complete.json'
    if not receipt.is_file():
        return False
    info = json.loads(receipt.read_text())
    return info['fingerprint'] == fingerprint and all((folder/p).is_file() and digest(folder/p) == h
                 for p, h in info['output_sha256'].items())


def diagnostic_plot(x, fd, events, tr, target):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(3, 1, figsize=(14, 9), layout='constrained')
    axs[0].imshow((x-x.mean())/np.maximum(x.std(), 1e-10), aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
    axs[0].set(title='Task + nuisance design', ylabel='Frame')
    axs[1].plot(np.arange(len(fd))*tr, fd); axs[1].set(ylabel='Raw FD (mm)')
    for condition, color in zip(events.trial_type.unique(), ['#d04b36', '#3064aa']):
        rows = events[events.trial_type.eq(condition)]
        axs[2].broken_barh(list(zip(rows.onset, rows.duration)), (0 if condition == events.trial_type.iloc[0] else 1, 0.8), facecolors=color, label=condition)
    axs[2].legend(); axs[2].set(xlabel='Seconds relative to retained BOLD', yticks=[])
    fig.savefig(target, dpi=120); plt.close(fig)


def render_maps(folder, recipe, subject, brain_axis, borders):
    from pfm_ridge_scene import build_scene
    view_dir = folder/'views'; view_dir.mkdir()
    records = []
    limit = str(recipe['render']['z_limit'])
    for path in sorted(folder.glob('*.dscalar.nii')):
        data = np.asarray(nib.load(str(path)).dataobj)
        for display in ['unthresholded', 'FDR05']:
            z = np.nan_to_num(data[3].copy(), nan=0.0, posinf=0.0, neginf=0.0)
            if display == 'FDR05':
                z[~np.isfinite(data[5]) | (data[5] > 0.05)] = 0
            name = path.name.split('.')[0]+'_'+display
            scalar = view_dir/(name+'_z.dscalar.nii'); colored = view_dir/(name+'.dscalar.nii')
            scene = view_dir/(name+'.scene'); png = view_dir/(name+'.png')
            write_scalar([z], ['z'], brain_axis, scalar)
            subprocess.run(['wb_command', '-cifti-palette', str(scalar), 'MODE_USER_SCALE', str(colored),
                '-pos-user', '0', limit, '-neg-user', '0', '-'+limit, '-palette-name', 'ROY-BIG-BL',
                '-disp-zero', 'false', '-disp-pos', 'true', '-disp-neg', 'true'], check=True)
            build_scene(ROOT/'res0urces/CrawlingSeedFC/FlatMaps+Inflated.scene', colored,
                subject/'anat/T1w/fsaverage_LR32k', subject/'anat/MNINonLinear/fsaverage_LR32k',
                subject.name, scene, 0, border_paths=borders, compact=True)
            tree = ET.parse(scene)
            for node in tree.iter('Object'):
                value = node.text or ''
                if node.get('Type') == 'pathName' and value.endswith(('.surf.gii', '.border', '.dscalar.nii')):
                    require(Path(value).is_file(), 'Scene resource missing: '+value)
            display_node = tree.find('.//Object[@Class="DisplayPropertiesBorders"]')
            if display_node is not None:
                for node in display_node:
                    name_attr = node.get('Name', '')
                    value = ('0.8' if name_attr == 'm_aboveSurfaceOffset' else
                             recipe['render']['border_color'] if name_attr.startswith('m_standardColorTypeIn') else
                             str(recipe['render']['border_width']) if name_attr.startswith('m_lineWidthIn') else
                             'true' if name_attr.startswith('m_unstretchedLinesStatusIn') else None)
                    if value is not None:
                        if list(node):
                            for item in node: item.text = value
                        else: node.text = value
            tree.write(scene, encoding='UTF-8', xml_declaration=True)
            env = dict(os.environ, QT_QPA_PLATFORM='offscreen')
            capture = subprocess.run(['wb_command', '-scene-capture-image', str(scene), '1', str(png),
                            '-size-width-height', '1800', '1080'], env=env, check=True, capture_output=True, text=True)
            (view_dir/(name+'.render.log')).write_text(capture.stdout+capture.stderr)
            require('ERRORS loading scene' not in capture.stdout+capture.stderr, 'Workbench reported incomplete scene loading')
            require(png.is_file() and png.stat().st_size > 1000, 'Empty map render')
            from matplotlib.image import imread
            require(np.ptp(imread(png)[..., :3]) > 0.05, 'Blank map render')
            records.append({'contrast': path.name.split('.')[0], 'display': display,
                            'image': str(png.relative_to(folder)), 'statistics': path.name})
    return records


def fit_run(subject, job, jobs=2):
    recipe_path, r, raw, event_path, timing_path, session, run = job
    timing = json.loads(timing_path.read_text())
    require(timing.get('alignment_status') in ('confirmed', 'assumed_authorized'), 'Timing alignment pending: '+str(timing_path))
    processed = subject/'func'/r['task']/('session_'+str(session))/('run_'+str(run))
    stem = r['file_prefix']; inputs = []
    for tag in r['input_tags']:
        inputs.extend([processed/(stem+'_'+tag+'.dtseries.nii'), processed/(stem+'_'+tag+'.nii.gz')])
    raw_bold = raw/(stem+'_S'+str(session)+'_R'+str(run)+'_E1.nii.gz')
    seg_path = subject/'anat/T1w/aparc+aseg.nii.gz'
    inputs.extend([recipe_path, event_path, timing_path, processed/'MCF.par', processed/'SliceTiming.txt', raw_bold, raw_bold.with_name(raw_bold.name.replace('.nii.gz', '.json')), seg_path, Path(__file__)])
    pfm = Path(os.environ.get('PFM_OUTDIR') or subject/'pfm')
    label = Path(os.environ.get('TASK_MODEL_NETWORK_LABEL') or pfm/'RidgePFM/RidgePFM_VTX.dlabel.nii')
    borders = []
    if r['render']['enabled']:
        require(os.environ.get('AtlasSpace', 'T1w') == 'T1w', 'Task rendering/tissue backend currently requires AtlasSpace=T1w')
        inputs += [ROOT/'lib/pfm_ridge_scene.py', ROOT/'res0urces/CrawlingSeedFC/FlatMaps+Inflated.scene']
        if r['render']['borders']:
            borders = [Path(str(label).replace('.dlabel.nii', '.'+h+'.border')) for h in ['L', 'R']]
            inputs += [label] + borders
        for d in [subject/'anat/T1w/fsaverage_LR32k', subject/'anat/MNINonLinear/fsaverage_LR32k']:
            inputs += sorted(d.glob('*.surf.gii'))
    if (processed/'rmVols.txt').exists(): inputs += [processed/'rmVols.txt']
    for p in inputs: require(p.is_file(), 'Missing model input: '+str(p))
    provenance = {str(p): digest(p) for p in inputs}
    versions = {'nilearn': nilearn.__version__, 'nibabel': nib.__version__, 'numpy': np.__version__, 'python': sys.version}
    fingerprint = hashlib.sha256(json.dumps({'inputs': provenance, 'versions': versions}, sort_keys=True).encode()).hexdigest()
    parent = Path(os.environ.get('TASK_MODEL_OUTDIR') or subject/'task_models')/r['task']/('session_'+str(session))/('run_'+str(run))
    parent.mkdir(parents=True, exist_ok=True); final = parent/fingerprint[:16]
    if receipt_ok(final, fingerprint):
        print('[task-model] reuse '+str(final), flush=True); return final
    require(not final.exists(), 'Existing result is incomplete or changed; preserve it and choose a new output root: '+str(final))
    folder = Path(tempfile.mkdtemp(prefix='.incomplete-', dir=parent))
    # Scenes contain absolute paths; build directly in the unique final folder and keep incomplete status until receipt exists.
    folder.rename(final); folder = final
    (folder/'recipe.yaml').write_text(recipe_path.read_text()); (folder/'timing.json').write_text(timing_path.read_text())
    (folder/'input_provenance.json').write_text(json.dumps({'sha256': provenance, 'versions': versions}, indent=2))
    axis = None; frame_info = None; summaries = []; all_records = []
    raw_n = nib.load(str(raw_bold)).shape[3]
    raw_tr = float(nib.load(str(raw_bold)).header.get_zooms()[3])
    motion = np.loadtxt(processed/'MCF.par', ndmin=2)
    segmentation = nib.load(str(seg_path))
    for tag in r['input_tags']:
        dest = folder/tag; dest.mkdir()
        image = nib.load(str(processed/(stem+'_'+tag+'.dtseries.nii')))
        time_axis = image.header.get_axis(0); brain = image.header.get_axis(1)
        require(isinstance(time_axis, nib.cifti2.SeriesAxis) and time_axis.unit == 'SECOND' and abs(time_axis.start) < 1e-6,
                'Expected uncensored CIFTI SeriesAxis starting at zero seconds')
        n, tr = image.shape[0], float(time_axis.step)
        if borders:
            require(nib.load(str(label)).header.get_axis(1) == brain, 'Network and task grayordinate geometries differ')
        require(np.isclose(tr, raw_tr, atol=1e-5), 'Raw and processed TR differ')
        removed = int((processed/'rmVols.txt').read_text().strip()) if (processed/'rmVols.txt').exists() else 0
        require(removed >= 0 and n+removed == raw_n, 'CIFTI frame count is not raw frames minus recorded removals')
        if axis is not None: require(brain == axis and frame_info == (n, tr), 'Input geometries or frame clocks differ')
        axis = brain; frame_info = (n, tr)
        events, full_events = load_events(event_path, r, timing, n, tr, removed)
        volume = nib.load(str(processed/(stem+'_'+tag+'.nii.gz')), keep_file_open=True)
        require(volume.ndim == 4 and volume.shape[3] == n and np.isclose(volume.header.get_zooms()[3], tr, atol=1e-5), 'Volume/CIFTI frame clocks differ')
        signals, counts = tissue_signals(volume, segmentation, r)
        reference, reference_qc = slice_reference(raw_bold.with_name(raw_bold.name.replace('.nii.gz', '.json')), processed/'SliceTiming.txt', tr)
        x, fd, qc = build_design(n, tr, events, motion, signals, r, reference)
        qc.update(reference_qc, sampling_reference_seconds=reference)
        x.to_csv(dest/'design.tsv', sep='\t', index=False); events.to_csv(dest/'events.tsv', sep='\t', index=False)
        np.savetxt(dest/'raw_fd.txt', fd); diagnostic_plot(x, fd, events, tr, dest/'design_motion.png')
        y = np.asarray(image.dataobj, dtype=np.float64)
        valid = np.isfinite(y).all(0) & (y.std(0) > 1e-8)
        require(valid.any(), 'No valid grayordinates')
        labels, results = run_glm(y[:, valid], x.values, noise_model='ar1', bins=100, n_jobs=jobs, random_state=1)
        raw_residual = np.zeros_like(y[:, valid]); valid_n = int(valid.sum())
        for group, result in results.items():
            take = labels == group
            raw_residual[:, take] = y[:, valid][:, take]-x.values@result.theta
        residual_var = np.mean(raw_residual**2, axis=0)
        a = raw_residual[:-1]-raw_residual[:-1].mean(0); b = raw_residual[1:]-raw_residual[1:].mean(0)
        residual_ac = np.sum(a*b, 0)/np.maximum(np.sqrt(np.sum(a*a,0)*np.sum(b*b,0)), 1e-20)
        for name, weights in r['contrasts'].items():
            c = np.zeros(x.shape[1])
            for condition, weight in weights.items(): c[x.columns.get_loc(condition)] = weight
            con = compute_contrast(labels, results, c, stat_type='t')
            ts = con.stat(); ps = 2*student_t.sf(np.abs(ts), con.dof)
            q = multipletests(ps, method='fdr_bh')[1]
            data = np.full((6, image.shape[1]), np.nan)
            for i, arr in enumerate([con.effect, np.sqrt(con.variance), ts, con.z_score(), ps, q]): data[i, valid] = arr
            write_scalar(data, ['effect', 'standard_error', 't', 'z', 'p_two_sided', 'q_FDR_two_sided'], brain, dest/(name+'.dscalar.nii'))
        residual_maps = np.full((2, image.shape[1]), np.nan); residual_maps[0, valid] = residual_var; residual_maps[1, valid] = residual_ac
        write_scalar(residual_maps, ['raw_residual_mean_square', 'raw_residual_lag1_correlation'], brain, dest/'residual_diagnostics.dscalar.nii')
        summary = dict(input_tag=tag, n_frames=n, tr_seconds=tr, removed_volumes=removed, valid_grayordinates=valid_n,
                       tissue_voxels=counts, condition_blocks=events.trial_type.value_counts().to_dict(),
                       unlogged_tail_seconds=float(n*tr-(full_events.onset+full_events.duration).max()),
                       noise_model='Nilearn AR1, 100 coefficient bins', beta_units='input BOLD units',
                       residual_lag1_median=float(np.median(residual_ac)), **qc)
        summary['task_motion_max_absolute_correlation'] = {condition: float(np.max(np.abs(np.nan_to_num(np.corrcoef(
            np.column_stack([x[condition].values, motion]).T)[0,1:])))) for condition in r['conditions']}
        (dest/'summary.json').write_text(json.dumps(summary, indent=2)); summaries.append(summary)
        if r['render']['enabled']:
            # Diagnostics is not a contrast and must not enter the map-render loop.
            diagnostic = dest/'residual_diagnostics.dscalar.nii'; temporary = dest/'residual_diagnostics.tmp'
            diagnostic.rename(temporary)
            try: records = render_maps(dest, r, subject, brain, borders)
            finally: temporary.rename(diagnostic)
            all_records.extend([dict(input_tag=tag, **record) for record in records])
    report(folder, r, timing, all_records, summaries)
    hashes = {str(p.relative_to(folder)): digest(p) for p in folder.rglob('*') if p.is_file()}
    (folder/'complete.json').write_text(json.dumps({'fingerprint': fingerprint, 'output_sha256': hashes}, indent=2))
    print('[task-model] complete '+str(folder), flush=True)
    return folder


def report(folder, recipe, timing, records, summaries):
    (folder/'results.json').write_text(json.dumps(records, indent=2))
    payload = json.dumps(records).replace('<', '\\u003c')
    title = html.escape(recipe['task']+' · bilateral task maps')
    body = '''<!doctype html><html><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><style>
body{font:16px/1.6 system-ui;max-width:1600px;margin:auto;padding:24px;background:#f6f8fa}section{background:white;padding:20px;margin:18px 0}select{padding:8px;font:inherit}.pair{display:grid;grid-template-columns:1fr 1fr;gap:16px}img{width:100%}pre{white-space:pre-wrap}@media(max-width:850px){.pair{grid-template-columns:1fr}}</style>'''
    body += '<h1>'+title+'</h1><section><p>'+html.escape(recipe['interpretation'])+'</p><p>Canonical HRF + time derivative; motion24 + raw FD spikes; input-specific WM/CSF; cosine drift; AR1. No extra smoothing or global signal regression. Each run fitted separately.</p>'
    body += '<p>Timing status: <strong>'+html.escape(timing['alignment_status'])+'</strong>. Offset: '+str(timing['onset_offset_seconds'])+' s. See timing provenance below. Unlogged scan tail is an explicit baseline assumption.</p><p>Both hemispheres: inflated lateral/medial views and flatmaps. White borders show the supplied individual network parcellation. This parcellation includes task data and is descriptive, not an independent network-specificity test.</p></section>'
    body += '<section><select id="contrast"></select> <select id="display"><option value="unthresholded">Unthresholded z</option><option value="FDR05">Two-sided FDR q≤0.05</option></select><p>Shared z scale: −'+str(recipe['render']['z_limit'])+' to +'+str(recipe['render']['z_limit'])+'. FDR is per contrast/input over valid grayordinates; no joint correction across contrasts. Thresholded map differences do not test an OCME-versus-MEICA difference.</p><div class="pair">'
    for i, tag in enumerate(recipe['input_tags']): body += '<div><h2>'+html.escape(tag)+'</h2><img id="map'+str(i)+'" alt="Bilateral inflated and flatmap views"><a id="stats'+str(i)+'">Statistical maps</a></div>'
    body += '</div></section><section><h2>Design and timing diagnostics</h2>'
    for tag in recipe['input_tags']: body += '<h3>'+html.escape(tag)+'</h3><img src="'+tag+'/design_motion.png"><a href="'+tag+'/design.tsv">Design matrix</a> · <a href="'+tag+'/summary.json">Fit summary</a>'
    body += '<details><summary>Timing provenance</summary><pre>'+html.escape(json.dumps(timing, indent=2))+'</pre></details><a href="recipe.yaml">Recipe</a> · <a href="input_provenance.json">Input provenance</a></section>'
    body += '<script>const data='+payload+''';const c=document.getElementById('contrast'),d=document.getElementById('display');[...new Set(data.map(x=>x.contrast))].forEach(x=>c.add(new Option(x,x)));function show(){['OCME','OCME+MEICA'].forEach((tag,i)=>{let r=data.find(x=>x.input_tag===tag&&x.contrast===c.value&&x.display===d.value);if(r){document.getElementById('map'+i).src=tag+'/'+r.image;document.getElementById('stats'+i).href=tag+'/'+r.statistics;}})}c.onchange=d.onchange=show;show();</script></html>'''
    (folder/'index.html').write_text(body)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--subject-dir', type=Path, required=True)
    parser.add_argument('--registry', type=Path, default=ROOT/'config/task_models')
    parser.add_argument('--start-session', type=int, default=1)
    parser.add_argument('--task', action='append', default=[])
    parser.add_argument('--mode', choices=['check', 'fit'], default='check')
    args = parser.parse_args()
    tasks = args.task or json.loads(os.environ.get('TASK_MODEL_TASKS_JSON', '[]'))
    jobs = discover(args.subject_dir, args.registry, tasks, args.start_session)
    if args.mode == 'check':
        checks = []
        for _, r, raw, events, timing, session, run in jobs:
            j = json.loads(timing.read_text()); path = raw/(r['file_prefix']+'_S'+str(session)+'_R'+str(run)+'_E1.nii.gz')
            image = nib.load(str(path)); n, tr = image.shape[3], float(image.header.get_zooms()[3])
            ev, df = load_events(events, r, j, n, tr)
            checks.append({'task': r['task'], 'session': session, 'run': run, 'blocks': ev.trial_type.value_counts().to_dict(), 'n_frames': n, 'tr': tr})
        print(json.dumps(checks, indent=2)); return
    results = [fit_run(args.subject_dir, job, int(os.environ.get('TASK_MODEL_JOBS', '2'))) for job in jobs]
    root = Path(os.environ.get('TASK_MODEL_OUTDIR') or args.subject_dir/'task_models'); root.mkdir(parents=True, exist_ok=True)
    links = ''.join('<li><a href="'+os.path.relpath(p/'index.html', root)+'">'+html.escape(str(p.relative_to(root)))+'</a></li>' for p in results)
    (root/'index.html').write_text('<!doctype html><meta charset="utf-8"><h1>Task modeling</h1><ul>'+links+'</ul>')


if __name__ == '__main__':
    main()
