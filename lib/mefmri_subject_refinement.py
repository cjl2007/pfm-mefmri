"""Manifest-driven subject refinement; outputs are separate from baseline pointers.

The pipeline adapter supplies explicit ACPC baseline chains and compatibility.
This module never edits source data or previously published reference versions.
"""
import argparse
import json
from pathlib import Path
import re
import tempfile
import numpy as np
from mefmri_refinement_fit import estimate_joint, estimate_frozen, command, render
from mefmri_refinement_geometry import rigid, refined_affine, warp_commands
from mefmri_refinement_store import digest, canonical, publish, verify, publication_lock


def verify_alignment(output, expected_version=None):
    """Validate a completed transform bundle before selecting it for a run."""
    output = Path(output).resolve()
    record = json.loads((output / 'alignment.json').read_text())
    if record.get('schema') != 1 or record.get('complete') is not True or record.get('strategy') != 'subject_refined':
        raise ValueError('Incomplete or unsupported alignment bundle')
    if expected_version is not None and record['version'] != expected_version:
        raise ValueError('Requested reference version differs from saved alignment')
    verify(record['reference_store'], record['version'], record['compatibility'])
    for path, checksum in record['source_sha256'].items():
        if digest(path) != checksum:
            raise ValueError('Alignment source changed: ' + path)
    for relative, checksum in record['files'].items():
        path = output / relative
        if path.is_symlink() or (path.resolve() != output and output not in path.resolve().parents) or digest(path) != checksum:
            raise ValueError('Alignment artifact changed: ' + relative)
    ids = [r['id'] for r in record['runs']]
    if len(set(ids)) != len(ids):
        raise ValueError('Duplicate run bindings')
    for run in record['runs']:
        for key in ('warp', 'correction', 'affine', 'reference_premat'):
            if run.get(key) and run[key] not in record['files']:
                raise ValueError('Selected transform is absent from artifact hashes')
    return record


def select_run(output, run_id, expected_version=None):
    record = verify_alignment(output, expected_version)
    found = [r for r in record['runs'] if r['id'] == run_id]
    if len(found) != 1:
        raise ValueError('Run is absent from selected alignment: ' + run_id)
    result = dict(found[0], version=record['version'], strategy=record['strategy'])
    for key in ('warp', 'correction', 'affine', 'reference_premat'):
        if result.get(key):
            result[key] = str(Path(output).resolve() / result[key])
    # Stored alongside final BOLD to reject stale downstream resume.
    result['binding_sha256'] = __import__('hashlib').sha256(canonical(dict(
        version=record['version'], run=found[0], source_sha256=record['source_sha256'],
        files={k: v for k, v in record['files'].items() if k.startswith(run_id + '/')}))).hexdigest()
    return result


def inspect_manifest(data):
    required = ('subject', 'func_dir', 'distortion_mode', 'acpc_grid', 'output_grid', 'output_space', 'runs', 'compatibility')
    if any(key not in data for key in required):
        raise ValueError('Incomplete refinement manifest')
    if data['distortion_mode'] not in ('medic', 'topup', 'direct_b0', 'none'):
        raise ValueError('Unsupported distortion mode')
    if data['output_space'] not in ('T1w', 'MNINonlinear'):
        raise ValueError('Unsupported output space')
    if data['output_space'] == 'MNINonlinear' and not data.get('standard_warp'):
        raise ValueError('Standard-space output needs the anatomy-to-standard warp')
    runs = data['runs']
    if not runs or len({r['id'] for r in runs}) != len(runs):
        raise ValueError('Explicit unique run membership is required')
    inputs = {str(Path(data[k]).resolve()) for k in ('acpc_grid', 'output_grid')}
    if data.get('standard_warp'):
        inputs.add(str(Path(data['standard_warp']).resolve()))
    members = []
    for run in runs:
        if not re.fullmatch(r'S[0-9]+_R[0-9]+', run['id']):
            raise ValueError('Invalid run identifier')
        if not isinstance(run.get('session'), int) or run['session'] < 1:
            raise ValueError('Positive session number required')
        if not all(run.get(k) for k in ('native', 'warp', 'target')):
            raise ValueError('Native reference, ACPC warp and motion target are required')
        if data['distortion_mode'] == 'none' and not run.get('affine'):
            raise ValueError('No-fieldmap path requires its baseline ACPC affine')
        hashes = {}
        for key in ('native', 'warp', 'target', 'premat', 'affine'):
            if run.get(key):
                path = Path(run[key]).resolve()
                hashes[key] = digest(path)
                inputs.add(str(path))
        members.append(dict(id=run['id'], session=run['session'], inputs=hashes))
        if run.get('metadata'):
            hashes['metadata'] = {}
            for name, value in run['metadata'].items():
                path = str(Path(value).resolve())
                inputs.add(path)
                hashes['metadata'][name] = digest(path)
    hashes = {p: digest(p) for p in sorted(inputs)}
    compatibility = dict(data['compatibility'], subject=data['subject'], func_dir=data['func_dir'],
                         distortion_mode=data['distortion_mode'], acpc_grid_sha256=digest(data['acpc_grid']))
    return members, hashes, compatibility


def execute(data, store, output, policy, version=None):
    members, before, compatibility = inspect_manifest(data)
    if policy not in ('rebuild', 'frozen') or (policy == 'frozen' and not version):
        raise ValueError('Select rebuild or frozen with an explicit reference version')
    if policy == 'rebuild' and version:
        raise ValueError('Rebuild creates its reference version; do not supply one')
    output, store = Path(output), Path(store)
    if output.exists():
        raise ValueError('Output already exists; validate its binding before reuse')
    output.parent.mkdir(parents=True, exist_ok=True)
    with publication_lock(output.parent / ('.' + output.name + '-lock')):
        work = Path(tempfile.mkdtemp(prefix='.' + output.name + '-building-', dir=output.parent))
        # Failed work remains for diagnosis, but has no completed output marker.
        if policy == 'rebuild':
            result = estimate_joint(data['runs'], data['acpc_grid'], work / 'fit')
            artifacts = {name + '.nii.gz': result[name] for name in ('template', 'weight', 'normalization_mask')}
            artifacts.update({r['id'] + '.mat': result['corrections'][r['id']] for r in data['runs']})
            settings = dict(method='leave_session_out_rigid_two_pass', anchor='none',
                            implementation={p.name: digest(p) for p in Path(__file__).parent.glob('mefmri_*refinement*.py')})
            if before != {p: digest(p) for p in before}:
                raise ValueError('Input changed before reference publication')
            version = publish(store, artifacts, members, compatibility, settings)
        manifest = verify(store, version, compatibility)
        artifact = store / version
        calibration = {m['id']: m for m in manifest['members']}
        corrections = {}
        incoming = []
        for run, member in zip(data['runs'], members):
            if run['id'] in calibration:
                if canonical(member) != canonical(calibration[run['id']]):
                    raise ValueError('Calibration run changed; explicitly rebuild the reference: ' + run['id'])
                corrections[run['id']] = str(artifact / (run['id'] + '.mat'))
            else:
                incoming.append(run)
        if incoming:
            result = estimate_frozen(incoming, data['acpc_grid'], artifact, work / 'incoming')
            corrections.update(result['corrections'])
        records = []
        for run in data['runs']:
            run_dir = work / run['id']
            run_dir.mkdir()
            correction = rigid(np.loadtxt(corrections[run['id']]))
            np.savetxt(run_dir / 'correction.mat', correction)
            reference_premat = None
            if run.get('premat'):
                reference_premat = run_dir / 'reference_to_intermediate.mat'
                np.savetxt(reference_premat, np.loadtxt(run['premat']))
            acpc = run_dir / 'intermediate_to_acpc.nii.gz'
            standard_args = {}
            if data['output_space'] == 'MNINonlinear':
                standard_args = dict(standard_warp=data['standard_warp'], standard_reference=data['output_grid'],
                                     standard_output=run_dir / 'intermediate_to_standard.nii.gz')
            for args in warp_commands(run['warp'], run_dir / 'correction.mat', data['acpc_grid'], acpc, **standard_args):
                command(args, work / 'commands.jsonl')
            selected = standard_args.get('standard_output', acpc)
            affine = None
            if data['distortion_mode'] == 'none':
                affine = run_dir / 'intermediate_to_acpc.mat'
                np.savetxt(affine, refined_affine(np.loadtxt(run['affine']), correction))
            # Estimation/QC reference only; never used as final BOLD source.
            render(run, data['acpc_grid'], run_dir / 'refined_reference.nii.gz', run_dir / 'correction.mat', work / 'commands.jsonl')
            records.append(dict(id=run['id'], motion_target=run['target'],
                                reference_premat=str(reference_premat.relative_to(work)) if reference_premat else None,
                                warp=str(selected.relative_to(work)), affine=str(affine.relative_to(work)) if affine else None,
                                correction=str((run_dir / 'correction.mat').relative_to(work))))
        if before != {p: digest(p) for p in before}:
            raise ValueError('Input changed while refinement was being estimated')
        verify(store, version, compatibility)
        result = dict(schema=1, complete=True, strategy='subject_refined', policy=policy,
                      xfms_namespace=data.get('xfms_namespace'),
                      output_space=data['output_space'], output_grid_sha256=digest(data['output_grid']),
                      version=version, reference_store=str(store.resolve()), source_sha256=before,
                      compatibility=compatibility, runs=records,
                      files={str(p.relative_to(work)): digest(p) for p in work.rglob('*') if p.is_file()})
        (work / 'alignment.json').write_text(json.dumps(result, indent=2) + '\n')
        work.rename(output)
        return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--store', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--policy', choices=['rebuild', 'frozen'], required=True)
    parser.add_argument('--version')
    args = parser.parse_args()
    execute(json.loads(Path(args.manifest).read_text()), args.store, args.output, args.policy, args.version)
