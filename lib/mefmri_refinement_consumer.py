"""Select immutable run transforms and record the alignment used by final BOLD."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
from mefmri_refinement_store import canonical, digest, require_run_binding
from mefmri_subject_refinement import select_run
from mefmri_refinement_manifest import read_pointer


def selection(run_dir, strategy, bundle=None, version=None):
    run_dir = Path(run_dir)
    session = int(run_dir.parent.name.split('session_', 1)[-1])
    run = int(run_dir.name.split('run_', 1)[-1])
    target = read_pointer(run_dir / 'IntermediateCoregTarget.txt')
    if strategy == 'subject_refined':
        if not bundle or not version:
            raise ValueError('Subject refinement requires an alignment bundle and explicit version')
        result = select_run(bundle, f'S{session}_R{run}', version)
        result['alignment_bundle'] = str(Path(bundle).resolve())
        if Path(result['motion_target']).resolve() != target.resolve():
            raise ValueError('Refinement motion target no longer matches baseline')
    elif strategy == 'independent':
        warp = read_pointer(run_dir / 'Intermediate2ACPCWarp.txt')
        affine_pointer = run_dir / 'Intermediate2ACPCMat.txt'
        affine = read_pointer(affine_pointer) if affine_pointer.exists() else None
        inputs = {str(p): digest(p) for p in (target, warp)}
        if affine:
            inputs[str(affine)] = digest(affine)
        result = dict(strategy=strategy, version=None, motion_target=str(target), warp=str(warp),
                      affine=str(affine) if affine else None,
                      binding_sha256=hashlib.sha256(canonical(inputs)).hexdigest())
    else:
        raise ValueError('Unknown alignment strategy: ' + strategy)
    return result


def binding(result):
    return {k: result[k] for k in ('strategy', 'version', 'binding_sha256')}


def check_completed(run_dir, result):
    path = Path(run_dir) / 'BOLDAlignment.json'
    if not path.is_file():
        if result['strategy'] == 'subject_refined':
            raise ValueError('No completed BOLD alignment record; rerun headmotion before downstream processing')
        return  # Legacy independent outputs predate alignment records.
    saved = json.loads(path.read_text())
    require_run_binding(saved['alignment'], binding(result))
    for name, checksum in saved['outputs'].items():
        path = Path(run_dir) / name
        if Path(name).name != name or digest(path) != checksum:
            raise ValueError('Completed BOLD output is missing or changed: ' + name)


def mark_completed(run_dir, result, prefix, expected_echoes=None, expected_frames=None):
    run_dir = Path(run_dir)
    outputs = sorted(run_dir.glob(prefix + '_E*_acpc.nii.gz'))
    if not outputs:
        raise ValueError('Cannot mark alignment complete without final BOLD outputs')
    if expected_echoes is not None:
        import nibabel as nib
        import numpy as np
        expected = {f'{prefix}_E{e}_acpc.nii.gz' for e in range(1, expected_echoes + 1)}
        if expected_echoes < 1 or {p.name for p in outputs} != expected:
            raise ValueError('Final BOLD echo set differs from the expected echoes')
        images = [nib.load(p) for p in outputs]
        shape = images[0].shape
        if len(shape) != 4 or min(shape) < 1:
            raise ValueError('Final BOLD must be a nonempty 4D image')
        if expected_frames is not None and shape[3] != expected_frames:
            raise ValueError('Final BOLD frame count differs from native input after discards')
        for image in images[1:]:
            if image.shape != shape or not np.allclose(image.affine, images[0].affine, atol=1e-5):
                raise ValueError('Final BOLD echoes have different dimensions or spatial grids')
    record = dict(alignment=binding(result), outputs={p.name: digest(p) for p in outputs})
    if result.get('alignment_bundle'):
        record['alignment_bundle'] = result['alignment_bundle']
    fd, name = tempfile.mkstemp(prefix='.BOLDAlignment-', dir=run_dir)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(record, stream, indent=2)
            stream.write('\n')
        os.replace(name, run_dir / 'BOLDAlignment.json')
    finally:
        if Path(name).exists():
            Path(name).unlink()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['select', 'check', 'mark'])
    parser.add_argument('--run-dir', required=True)
    parser.add_argument('--strategy', choices=['independent', 'subject_refined'], required=True)
    parser.add_argument('--bundle')
    parser.add_argument('--version')
    parser.add_argument('--prefix', default='Rest')
    args = parser.parse_args()
    result = selection(args.run_dir, args.strategy, args.bundle, args.version)
    if args.action == 'check':
        check_completed(args.run_dir, result)
    elif args.action == 'mark':
        import nibabel as nib
        run_dir = Path(args.run_dir).resolve()
        expected_echoes = len((run_dir / 'TE.txt').read_text().split())
        raw_dir = run_dir.parents[2] / 'unprocessed' / run_dir.parents[1].name / run_dir.parent.name / run_dir.name
        raw = list(raw_dir.glob(args.prefix + '_S*_R*_E1.nii.gz'))
        if len(raw) != 1:
            raise ValueError('Cannot verify final frame count without a unique native E1')
        native_shape = nib.load(raw[0]).shape
        if len(native_shape) != 4:
            raise ValueError('Native E1 must be 4D')
        discard_file = run_dir / 'rmVols.txt'
        discard = int(discard_file.read_text().strip()) if discard_file.exists() else 0
        if discard < 0 or discard >= native_shape[3]:
            raise ValueError('Invalid number of discarded frames')
        mark_completed(args.run_dir, result, args.prefix, expected_echoes, native_shape[3] - discard)
    else:
        # Line-delimited fields are consumed as data, never evaluated as shell.
        for key in ('motion_target', 'warp', 'affine'):
            value = result.get(key) or ''
            if '\n' in value or '\r' in value:
                raise ValueError('Newlines are not supported in transform paths')
            print(value)
