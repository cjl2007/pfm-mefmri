"""Discover explicit baseline chains without changing pipeline pointer files."""
from pathlib import Path
import json
import re
import os
import nibabel as nib
import numpy as np
from mefmri_refinement_fit import command
from mefmri_refinement_store import digest


def read_pointer(path):
    lines = Path(path).read_text().splitlines()
    if len(lines) != 1 or not lines[0].strip():
        raise ValueError('Expected one path in ' + str(path))
    target = Path(lines[0])
    if not target.is_absolute() or not target.is_file():
        raise ValueError('Missing or nonabsolute baseline pointer: ' + str(path))
    return target


def acquisition(directory):
    """Keep incompatible EPI contrast/acquisition groups out of one template."""
    directory = Path(directory)
    paths = {name: directory / name for name in ('TE.txt', 'TR.txt', 'PE.txt')}
    tes = [float(v) for v in paths['TE.txt'].read_text().split()]
    trs = [float(v) for v in paths['TR.txt'].read_text().split()]
    if not tes or len(trs) != 1 or not all(np.isfinite(v) and v > 0 for v in tes + trs):
        raise ValueError('Invalid acquisition metadata for refinement: ' + str(directory))
    pe = paths['PE.txt'].read_text().strip()
    # find_epi_params writes FSL axis names, not BIDS axis names.
    pe = {'x': 'i', '-x': 'i-', 'y': 'j', '-y': 'j-', 'z': 'k', '-z': 'k-'}.get(pe, pe)
    if pe not in ('', 'i', 'i-', 'j', 'j-', 'k', 'k-'):
        raise ValueError('Invalid phase-encoding direction: ' + str(directory))
    signature = dict(te_ms=[round(v, 6) for v in tes], tr_seconds=round(trs[0], 6), phase_encoding=pe,
                     sbref_combination=os.environ.get('SBREF_ECHO_COMBINATION', 'mean'),
                     sbref_max_te_ms=os.environ.get('SBREF_MAX_TE_MS', '60'))
    return signature, {name: str(path) for name, path in paths.items()}


def comparison_acquisition(signature):
    """Opt-in, explicit TE schedules for optimally combined references.

    Default comparisons remain exact. This permits known TE variation only;
    TR, PE, echo count, reference construction, and geometry still must match.
    Original per-run TE signatures and metadata hashes are always retained.
    """
    policy = os.environ.get('COREG_REFERENCE_TE_POLICY', 'strict')
    if policy == 'strict':
        return signature
    if policy != 'allow_list':
        raise ValueError('Invalid COREG_REFERENCE_TE_POLICY')
    if signature['sbref_combination'] not in ('t2s', 'paid'):
        raise ValueError('TE schedule allow-list requires optimally combined SBRefs (t2s or paid)')
    try:
        schedules = json.loads(os.environ.get('COREG_REFERENCE_TE_SCHEDULES_MS', ''))
        if not isinstance(schedules, list) or not schedules:
            raise ValueError('Empty schedule list')
        normalized = []
        for values in schedules:
            if not isinstance(values, list) or len(values) != len(signature['te_ms']):
                raise ValueError('Echo count differs')
            if any(isinstance(v, bool) or not isinstance(v, (int, float)) for v in values):
                raise ValueError('Invalid TE')
            values = [round(float(v), 6) for v in values]
            if not all(np.isfinite(v) and v > 0 for v in values) or any(b <= a for a, b in zip(values, values[1:])):
                raise ValueError('TEs must increase and be finite/positive')
            normalized.append(values)
        schedules = [list(values) for values in sorted(set(tuple(v) for v in normalized))]
    except (ValueError, TypeError) as exc:
        raise ValueError('Invalid COREG_REFERENCE_TE_SCHEDULES_MS: ' + str(exc)) from exc
    if signature['te_ms'] not in schedules:
        raise ValueError('Run TE schedule is not explicitly approved for refinement')
    result = dict(signature)
    result.pop('te_ms')
    result.update(te_policy='allow_list', te_schedules_ms=schedules)
    return result


def discover(subject_dir, func_dir, xfms_dir, mode, space, output_grid, work,
             start_session=1, selected_ids=None):
    subject = Path(subject_dir).resolve()
    for name in (func_dir, xfms_dir):
        if not re.fullmatch('[A-Za-z0-9_.-]+', name) or name in ('.', '..'):
            raise ValueError('Invalid functional directory name')
    if mode not in ('medic', 'topup', 'direct_b0', 'none') or space not in ('T1w', 'MNINonlinear'):
        raise ValueError('Invalid fieldmap mode or target space')
    work = Path(work)
    work.mkdir(parents=True, exist_ok=False)
    xfms = subject / 'func/xfms' / xfms_dir
    grid = xfms / 'T1w_acpc_brain_func.nii.gz'
    if not grid.is_file():
        raise ValueError('Missing ACPC functional grid')
    runs = []
    geometries = []
    acquisitions = []
    for directory in (subject / 'func' / func_dir).glob('session_*/run_*'):
        match = re.fullmatch(r'session_(\d+)/run_(\d+)', str(directory.relative_to(subject / 'func' / func_dir)))
        if match is None:
            continue
        session, number = map(int, match.groups())
        tag = f'S{session}_R{number}'
        if session < start_session or (selected_ids is not None and tag not in selected_ids):
            continue
        target = read_pointer(directory / 'IntermediateCoregTarget.txt')
        selected = read_pointer(directory / 'Intermediate2ACPCWarp.txt')
        # Resolve ACPC counterpart explicitly; never apply a rigid ACPC matrix
        # after a precomposed nonlinear standard-space field.
        if space == 'MNINonlinear':
            if '2nonlin_' not in selected.name:
                raise ValueError('Cannot identify ACPC counterpart of selected standard warp')
            baseline = selected.with_name(selected.name.replace('2nonlin_', '2acpc_', 1))
        else:
            baseline = selected
        if not baseline.is_file():
            raise ValueError('Missing baseline ACPC warp: ' + str(baseline))
        premat = None
        if mode == 'medic':
            native = read_pointer(directory / 'MEDICUnwarpedSBRef.txt')
            expected = xfms / 'AvgMEDICSBref.nii.gz'
            if target.resolve() != expected.resolve():
                raise ValueError('Unsupported MEDIC motion-target semantics')
            premat = work / (tag + '_reference_to_intermediate.mat')
            command(['flirt', '-in', native, '-ref', target, '-omat', premat,
                     '-dof', '6', '-cost', 'normcorr', '-nosearch'], work / 'commands.jsonl')
        else:
            native = directory / 'SBref.nii.gz'
            if target.resolve() != native.resolve():
                if target.resolve() != (xfms / 'AvgSBref.nii.gz').resolve():
                    raise ValueError('Unrecognized motion target; explicit baseline manifest required')
                # Follow the same QA directory contract as baseline coreg.
                qa = Path(os.environ.get('COREG_QA_DIR', str(subject / 'func' / func_dir /
                          os.environ.get('QC_OUTPUT_DIRNAME', 'QC') / 'Coreg')))
                premat = qa / ('SBref2AvgSBref_' + tag + '.mat')
                if not premat.is_file():
                    raise ValueError('Missing reference-to-average prematrix: ' + str(premat))
        image = nib.load(native)
        if len(image.shape) != 3:
            raise ValueError('Native reference must be 3D')
        geometries.append(dict(shape=list(image.shape), zooms=[round(float(v), 5) for v in image.header.get_zooms()],
                               orientation=list(nib.aff2axcodes(image.affine))))
        record = dict(id=tag, session=session, run=number, native=str(native), warp=str(baseline),
                      premat=str(premat) if premat else None, target=str(target),
                      selected_baseline_warp=str(selected), directory=str(directory))
        signature, record['metadata'] = acquisition(directory)
        record['acquisition'] = signature
        acquisitions.append(comparison_acquisition(signature))
        if mode == 'none':
            record['affine'] = str(read_pointer(directory / 'Intermediate2ACPCMat.txt'))
        runs.append(record)
    if not runs:
        raise ValueError('No eligible run references')
    if selected_ids is not None and set(selected_ids) != {r['id'] for r in runs}:
        raise ValueError('Explicit membership includes missing/excluded runs')
    if any(g != geometries[0] for g in geometries[1:]):
        raise ValueError('Mixed reference geometry requires separate refinement groups')
    if any(a != acquisitions[0] for a in acquisitions[1:]):
        raise ValueError('Mixed acquisition metadata requires separate refinement groups')
    runs.sort(key=lambda r: (r['session'], r['run']))
    anatomy = subject / 'anat/T1w/T1w_acpc_dc_restore_brain.nii.gz'
    data = dict(subject=subject.name, func_dir=func_dir, distortion_mode=mode,
                acpc_grid=str(grid), output_grid=str(Path(output_grid).resolve()), output_space=space, runs=runs,
                compatibility=dict(reference_geometry=geometries[0], acquisition=acquisitions[0], anatomy_sha256=digest(anatomy)))
    if space == 'MNINonlinear':
        data['standard_warp'] = str(subject / 'anat/MNINonLinear/xfms/acpc_dc2standard.nii.gz')
    return data
