"""Pipeline-facing preparation/resume of explicit subject alignment bundles."""
import argparse
import json
from pathlib import Path
import tempfile
from mefmri_refinement_manifest import discover
from mefmri_subject_refinement import execute, verify_alignment
from mefmri_refinement_consumer import selection, check_completed
from mefmri_refinement_store import digest
from mefmri_refinement_namespace import reserve


def validate_membership(subject_dir, func_dir, start_session, members):
    if not members:
        return None
    selected = Path(members).read_text().splitlines()
    if not selected or any(not x.strip() or x != x.strip() for x in selected) or len(set(selected)) != len(selected):
        raise ValueError('Membership requires unique nonempty run IDs, one per line')
    eligible = set()
    for run in (Path(subject_dir) / 'func' / func_dir).glob('session_*/run_*'):
        try:
            session = int(run.parent.name.split('session_', 1)[-1])
            number = int(run.name.split('run_', 1)[-1])
        except ValueError:
            continue
        if run.is_dir() and session >= start_session:
            eligible.add(f'S{session}_R{number}')
    if set(selected) != eligible:
        raise ValueError('Membership must match every run selected by START_SESSION; '
                         'headmotion processes whole sessions. Use an isolated subset copy for pilot subsampling.')
    return selected


def prepare(subject_dir, func_dir, xfms_dir, mode, space, grid, policy, version,
            bundle, start_session, members=None, resume=False, reference_root=None):
    subject_dir = Path(subject_dir).resolve()
    if bundle:
        result = verify_alignment(bundle, version)
        if result['compatibility']['subject'] != subject_dir.name or result['compatibility']['func_dir'] != func_dir:
            raise ValueError('Alignment bundle belongs to a different subject or functional group')
        if result['compatibility']['distortion_mode'] != mode:
            raise ValueError('Alignment bundle belongs to a different distortion path')
        if result.get('output_space') != space or result.get('output_grid_sha256') != digest(grid):
            raise ValueError('Alignment bundle target space/grid differs from requested output')
        return str(Path(bundle).resolve()), result['version'], result['xfms_namespace']
    if resume:
        raise ValueError('Resuming subject-refined headmotion requires the saved COREG_ALIGNMENT_BUNDLE')
    if policy == 'rebuild' and start_session != 1:
        raise ValueError('Rebuild requires START_SESSION=1 and explicit complete membership; use frozen for incoming sessions')
    root = Path(reference_root) if reference_root else subject_dir / 'func/xfms' / xfms_dir / 'SubjectRefinement'
    root.mkdir(parents=True, exist_ok=True)
    invocation = Path(tempfile.mkdtemp(prefix='invocation-', dir=root))
    selected_ids = validate_membership(subject_dir, func_dir, start_session, members)
    data = discover(subject_dir, func_dir, xfms_dir, mode, space, grid,
                    invocation / 'discovery', start_session, selected_ids)
    data['xfms_namespace'] = xfms_dir
    (invocation / 'inputs.json').write_text(json.dumps(data, indent=2) + '\n')
    output = invocation / 'alignment'
    result = execute(data, root / 'references', output, policy, version)
    verify_alignment(output, result['version'])
    return str(output), result['version'], xfms_dir


def check_downstream(subject_dir, func_dir, start_session, strategy, bundle, version):
    subject_dir = Path(subject_dir)
    count = 0
    selected_ids = None
    if strategy == 'subject_refined':
        current = verify_alignment(bundle, version)
        selected_ids = {r['id'] for r in current['runs']}
    for run in (subject_dir / 'func' / func_dir).glob('session_*/run_*'):
        try:
            session = int(run.parent.name.split('session_', 1)[-1])
            int(run.name.split('run_', 1)[-1])
        except ValueError:
            continue
        if session < start_session:
            continue
        if strategy == 'independent' and not (run / 'BOLDAlignment.json').exists():
            continue  # Legacy independent processing has no refinement record.
        run_bundle = bundle
        if selected_ids is not None and f'S{session}_R{int(run.name.split("run_", 1)[-1])}' not in selected_ids:
            # Frozen incoming batches have separate bundles but one immutable
            # reference version. Validate older runs against their own bundle
            # when downstream stages intentionally include all sessions.
            saved = json.loads((run / 'BOLDAlignment.json').read_text())
            if saved.get('alignment', {}).get('strategy') != strategy or saved['alignment'].get('version') != version:
                raise ValueError('Mixed alignment strategies/reference versions in downstream selection: ' + str(run))
            run_bundle = saved.get('alignment_bundle')
            if not run_bundle:
                raise ValueError('Missing prior alignment bundle for frozen downstream reuse: ' + str(run))
            prior = verify_alignment(run_bundle, version)
            if prior['compatibility'] != current['compatibility']:
                raise ValueError('Prior frozen batch has incompatible reference provenance')
        result = selection(run, strategy, run_bundle, version)
        check_completed(run, result)
        count += 1
    if not count and strategy == 'subject_refined':
        raise ValueError('No selected runs for downstream alignment validation')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'check', 'reserve'])
    parser.add_argument('--subject-dir', required=True)
    parser.add_argument('--func-dir', required=True)
    parser.add_argument('--xfms-dir', required=True)
    parser.add_argument('--mode', choices=['topup', 'direct_b0', 'medic', 'none'], required=True)
    parser.add_argument('--space', choices=['T1w', 'MNINonlinear'], required=True)
    parser.add_argument('--grid', required=True)
    parser.add_argument('--policy', choices=['rebuild', 'frozen'], default='rebuild')
    parser.add_argument('--strategy', choices=['independent', 'subject_refined'], default='independent')
    parser.add_argument('--version')
    parser.add_argument('--bundle')
    parser.add_argument('--members')
    parser.add_argument('--start-session', type=int, default=1)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--reference-root')
    args = parser.parse_args()
    if args.action == 'reserve':
        if args.bundle:
            raise ValueError('Do not rerun coregistration with an existing bundle; resume from headmotion')
        validate_membership(args.subject_dir, args.func_dir, args.start_session, args.members)
        namespace, root = reserve(args.subject_dir, args.xfms_dir, args.policy, args.version,
                                  args.start_session, args.reference_root)
        print(namespace)
        print(root)
    elif args.action == 'check':
        check_downstream(args.subject_dir, args.func_dir, args.start_session, args.strategy, args.bundle, args.version)
    elif args.strategy == 'subject_refined':
        bundle, version, namespace = prepare(args.subject_dir, args.func_dir, args.xfms_dir, args.mode, args.space,
                                  args.grid, args.policy, args.version, args.bundle, args.start_session,
                                  args.members, args.resume, args.reference_root)
        print(bundle)
        print(version)
        print(namespace)
