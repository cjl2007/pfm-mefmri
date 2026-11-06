#!/usr/bin/env python3
"""Select and prepare one subject-level PFM input dataset with provenance."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    os.replace(temporary, path)


def signature(path):
    path = Path(path)
    stat = path.stat()
    return {'path': str(path.resolve()), 'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns}


def select(tasks, requested):
    if not tasks or len({task for task, _ in tasks}) != len(tasks):
        raise ValueError('Functional task names must be nonempty and unique')
    selected = []
    for value in requested:
        exact = [(task, prefix) for task, prefix in tasks if value == task]
        matches = exact or [(task, prefix) for task, prefix in tasks
                            if value.lower() in (task.lower(), prefix.lower())]
        if len(matches) != 1:
            raise ValueError('PFM input {!r} must identify exactly one configured task; matches: {}'.format(value, matches))
        task, prefix = matches[0]
        if any(row['task'] == task for row in selected):
            raise ValueError('Duplicate PFM input task: ' + task)
        selected.append({'task': task, 'prefix': prefix})
    if not selected:
        raise ValueError('PFM_INPUT_TASKS must contain at least one configured task')
    return selected


def prepare(args):
    subject = Path(args.subject_dir).resolve()
    root = Path(args.root).resolve()
    specs = json.loads(args.spec)
    if args.external_cifti:
        import nibabel as nib
        source = Path(args.external_cifti).expanduser().resolve()
        image = nib.load(str(source))
        if len(image.shape) != 2:
            raise ValueError('Explicit PFM_INPUT_CIFTI must be a time-by-grayordinate CIFTI')
        source_record = signature(source)
        fingerprint = hashlib.sha256(json.dumps({'source': source_record, 'tasks': specs}, sort_keys=True).encode()).hexdigest()
        path = root / 'dataset.json'
        previous = json.loads(path.read_text()) if path.is_file() else None
        if previous and previous.get('input_fingerprint') != fingerprint and not args.rebuild:
            raise ValueError('Explicit PFM input changed. Restart from concat to rebuild the dataset.')
        atomic_json(path, {'schema_version': 1, 'status': 'prepared', 'subject': subject.name,
                          'tasks': specs, 'input_mode': 'external_cifti_override',
                          'task_membership_verified': False, 'sources': [source_record],
                          'input_fingerprint': fingerprint, 'n_frames': int(image.shape[0]),
                          'functional_reference_task': specs[0]['task'],
                          'created_utc': datetime.now(timezone.utc).isoformat(),
                          'mapping_settings': {k: v for k, v in os.environ.items() if k.startswith('PFM_') and k != 'PFM_CONFIG_SNAPSHOT'}})
        return
    source_stats = []
    explicit_tag = os.environ.get('PFM_REQUESTED_INPUT_TAG', '')
    concat_subdir = os.environ.get('CONCAT_OUT_SUBDIR', 'ConcatenatedCiftis')
    from pfm_module import discover_runs, resolve_run_cifti
    for spec in specs:
        summary_path = subject / 'func' / spec['task'] / concat_subdir / 'PFMPrepSummary.json'
        if explicit_tag:
            spec['input_tag'] = explicit_tag
        elif summary_path.is_file():
            spec['input_tag'] = json.loads(summary_path.read_text())['input_tag']
        else:
            raise ValueError('Missing per-task concatenation metadata: {}. Run through concat first, or set PFM_INPUT_TAG explicitly.'.format(summary_path))
        runs = discover_runs(subject, spec['task'], args.start_session)
        if not runs:
            raise ValueError('No processed runs for selected PFM task: ' + spec['task'])
        for session, run, directory in runs:
            cifti = resolve_run_cifti(directory, spec['prefix'], spec['input_tag'])
            if cifti is None:
                raise ValueError('Missing selected PFM CIFTI in ' + str(directory))
            motion = directory / 'MCF.par'
            if not motion.is_file():
                motion = directory / 'mcf.par'
            source_stats.append({'task': spec['task'], 'session': session, 'run': run,
                                 'cifti': signature(cifti), 'motion': signature(motion),
                                 'tr': signature(directory / 'TR.txt')})
    knobs = {name: os.environ.get(name, default) for name, default in (
        ('PFM_FD_THRESHOLD', os.environ.get('CONCAT_FD_THRESHOLD', '0.3')),
        ('CONCAT_CENSOR_BY_FD', '1'), ('CONCAT_DEMEAN_RUNS', '1'),
        ('CONCAT_VAR_NORM_RUNS', '0'), ('CONCAT_VAR_NORM_EPS', '1e-8'))}
    payload = {'schema_version': 1, 'subject': subject.name, 'tasks': specs,
               'start_session': args.start_session, 'sources': source_stats, 'preparation': knobs}
    fingerprint = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    manifest_path = root / 'dataset.json'
    previous = json.loads(manifest_path.read_text()) if manifest_path.is_file() else None
    mapping_settings = {k: v for k, v in os.environ.items() if k.startswith('PFM_') and k != 'PFM_CONFIG_SNAPSHOT'}
    implementation = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                      for name in ('pfm_module.py', 'mefmri_pfm_dataset.py')}
    output = root / os.environ.get('PFM_CONCAT_OUT_SUBDIR', 'ConcatenatedCiftis')
    threshold = knobs['PFM_FD_THRESHOLD'].replace('.', 'p')
    censored = output / ('PFM_Combined_Concatenated+FDlt' + threshold + '.dtseries.nii')
    if previous and previous.get('input_fingerprint') == fingerprint and censored.is_file() and not args.rebuild:
        previous.update(status='prepared', mapping_settings=mapping_settings, implementation_sha256=implementation)
        previous.pop('completed_utc', None)
        previous.pop('outputs', None)
        atomic_json(manifest_path, previous)
        print('[pfm-dataset] reusing matching prepared inputs: ' + str(censored), flush=True)
        return
    if previous and previous.get('input_fingerprint') != fingerprint and not args.rebuild:
        raise ValueError('PFM input selection or source files changed. Restart from concat to rebuild the subject-level dataset; existing maps must not be resumed with different inputs.')
    if previous and previous.get('status') == 'complete' and not args.rebuild and not censored.is_file():
        raise ValueError('Prepared PFM input is missing; restart from concat instead of resuming stale maps.')
    root.mkdir(parents=True, exist_ok=True)
    running = {**payload, 'input_fingerprint': fingerprint, 'status': 'preparing',
               'created_utc': datetime.now(timezone.utc).isoformat()}
    atomic_json(manifest_path, running)
    command = [sys.executable, str(Path(__file__).with_name('pfm_module.py')),
               '--subject', subject.name, '--study-folder', str(subject.parent),
               '--start-session', str(args.start_session), '--func-dirname', specs[0]['task'],
               '--func-file-prefix', 'PFM', '--input-tag', 'Combined',
               '--input-specs-json', json.dumps(specs), '--output-dir', str(output),
               '--fd-threshold', knobs['PFM_FD_THRESHOLD'],
               '--censor-by-fd', knobs['CONCAT_CENSOR_BY_FD'],
               '--demean-runs', knobs['CONCAT_DEMEAN_RUNS'],
               '--var-norm-runs', knobs['CONCAT_VAR_NORM_RUNS'],
               '--var-norm-eps', knobs['CONCAT_VAR_NORM_EPS']]
    try:
        subprocess.run(command, check=True)
    except Exception:
        atomic_json(manifest_path, {**running, 'status': 'failed'})
        raise
    summary = json.loads((output / 'PFMPrepSummary.json').read_text())
    config = os.environ.get('PFM_CONFIG_SNAPSHOT', '')
    payload = {**running, 'status': 'prepared', 'runs': summary['run_scan_map'],
               'n_frames': summary['n_timepoints_concat'], 'n_retained': summary['n_timepoints_kept'],
               'functional_reference_task': specs[0]['task'],
               'functional_reference_xfms': os.environ.get('PFM_FUNC_XFMS_DIRNAME', specs[0]['task']),
               'atlas_space': os.environ.get('AtlasSpace', 'T1w'),
               'config_sha256': hashlib.sha256(config.encode()).hexdigest(),
               'mapping_settings': mapping_settings, 'implementation_sha256': implementation}
    atomic_json(manifest_path, payload)


def finish(args):
    root = Path(args.root)
    path = root / 'dataset.json'
    if not path.is_file():
        return  # Historical standalone callers do not yet have a dataset manifest.
    value = json.loads(path.read_text())
    outputs = []
    analyses = {}
    strategy = os.environ.get('PFM_STRATEGY', 'ridgepfm').lower()
    directories = ('RidgePFM', 'Infomap') if strategy == 'all' else ('Infomap',) if strategy == 'infomap' else ('RidgePFM',)
    for directory in directories:
        settings_file = root / directory / 'analysis.json'
        if settings_file.is_file():
            analyses[directory] = json.loads(settings_file.read_text())
        for pattern in ('*.dlabel.nii', '*.dtseries.nii', '*.dscalar.nii'):
            for file in (root / directory).glob(pattern):
                outputs.append({**signature(file), 'sha256': hashlib.sha256(file.read_bytes()).hexdigest()})
    value.update(status='complete', completed_utc=datetime.now(timezone.utc).isoformat(), outputs=outputs, analyses=analyses)
    atomic_json(path, value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest='command', required=True)
    p = subs.add_parser('select')
    p.add_argument('--task', action='append', nargs=2, required=True)
    p.add_argument('--select', action='append', required=True)
    p = subs.add_parser('prepare')
    p.add_argument('--subject-dir', required=True)
    p.add_argument('--root', required=True)
    p.add_argument('--spec', required=True)
    p.add_argument('--start-session', type=int, default=1)
    p.add_argument('--rebuild', type=int, choices=(0, 1), default=0)
    p.add_argument('--external-cifti', default='')
    p = subs.add_parser('finish')
    p.add_argument('--root', required=True)
    p = subs.add_parser('verify-selection')
    p.add_argument('--root', required=True)
    p.add_argument('--spec', required=True)
    args = parser.parse_args()
    try:
        if args.command == 'select':
            print(json.dumps(select(args.task, args.select)))
        elif args.command == 'prepare':
            prepare(args)
        elif args.command == 'finish':
            finish(args)
        else:
            path = Path(args.root) / 'dataset.json'
            if path.is_file():
                manifest = json.loads(path.read_text())
                requested = [(r['task'], r['prefix']) for r in json.loads(args.spec)]
                recorded = [(r['task'], r['prefix']) for r in manifest['tasks']]
                if requested != recorded or manifest.get('status') != 'complete':
                    raise ValueError('Requested PFM inputs do not match a complete dataset. Rebuild PFM before resuming downstream stages.')
    except (ValueError, KeyError, FileNotFoundError) as exc:
        parser.exit(2, 'ERROR: ' + str(exc) + '\n')


if __name__ == '__main__':
    main()
