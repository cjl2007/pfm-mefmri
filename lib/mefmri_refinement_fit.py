"""Subject-reference estimation, matching the ME01 two-pass pilot.

Input run records name a native reference, intermediate-to-ACPC warp and
optional reference-to-intermediate prematrix. For MEDIC, the reference is the
unwarped SBRef used only for estimation; final BOLD still uses native echoes.
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import nibabel as nib
import numpy as np
from scipy.ndimage import binary_erosion, gaussian_filter
from mefmri_refinement_geometry import rigid, center_corrections


def command(args, log):
    executable = shutil.which(args[0])
    if executable is None:
        executable = str(Path(os.environ.get('FSLDIR', '/usr/local/fsl')) / 'bin' / args[0])
    proc = subprocess.run([executable] + list(map(str, args[1:])),
                          env=dict(os.environ, FSLOUTPUTTYPE='NIFTI_GZ'),
                          capture_output=True, text=True)
    with Path(log).open('a') as stream:
        stream.write(json.dumps(dict(argv=list(map(str, args)), returncode=proc.returncode,
                                     stdout=proc.stdout, stderr=proc.stderr)) + '\n')
    if proc.returncode:
        raise RuntimeError('FSL command failed: ' + str(args[0]) + '\n' + proc.stderr)


def save(array, reference, path):
    nib.save(nib.Nifti1Image(np.asarray(array, dtype=np.float32), reference.affine, reference.header), path)


def normalized(path, reference, mask):
    image = nib.load(path)
    if image.shape != reference.shape or not np.allclose(image.affine, reference.affine, atol=1e-5):
        raise ValueError('Reference image grid mismatch: ' + str(path))
    data = image.get_fdata(dtype=np.float32)
    if not np.isfinite(data).all():
        raise ValueError('Nonfinite reference intensity: ' + str(path))
    selected = data[mask & (data > 0)]
    if selected.size < .5 * mask.sum():
        raise ValueError('Insufficient positive reference coverage')
    return data / np.median(selected)


def render(run, grid, output, correction, log):
    args = ['applywarp', '--in=' + str(run['native']), '--ref=' + str(grid),
            '--warp=' + str(run['warp']), '--out=' + str(output),
            '--interp=spline', '--datatype=float']
    if run.get('premat'):
        args.append('--premat=' + str(run['premat']))
    if correction is not None:
        args.append('--postmat=' + str(correction))
    command(args, log)


def fit(source, target, weight, identity, output, log):
    command(['flirt', '-in', source, '-ref', target, '-omat', output, '-dof', '6',
             '-cost', 'normcorr', '-nosearch', '-init', identity, '-refweight', weight], log)
    return rigid(np.loadtxt(output))


def prepare(runs, grid, work):
    work = Path(work)
    # Never inherit existence-only experimental caches.
    work.mkdir(parents=True, exist_ok=False)
    reference = nib.load(grid)
    if len(reference.shape) != 3:
        raise ValueError('ACPC reference must be 3D')
    mask = binary_erosion(reference.get_fdata() > 0, iterations=2)
    if mask.sum() < 100:
        raise ValueError('Insufficient ACPC brain mask')
    save(mask, reference, work / 'normalization_mask.nii.gz')
    save(gaussian_filter(mask.astype(float), 1), reference, work / 'weight.nii.gz')
    np.savetxt(work / 'identity.mat', np.eye(4))
    baseline = []
    for i, run in enumerate(runs):
        out = work / f'baseline_{i}.nii.gz'
        render(run, grid, out, None, work / 'commands.jsonl')
        data = normalized(out, reference, mask)
        save(data, reference, work / f'normalized_{i}.nii.gz')
        baseline.append(data)
    return reference, mask, np.stack(baseline)


def estimate_joint(runs, grid, work):
    if not runs or len({r['id'] for r in runs}) != len(runs):
        raise ValueError('Unique explicit run membership is required')
    groups = [r['session'] for r in runs]
    sessions = sorted(set(groups))
    if len(sessions) < 2:
        raise ValueError('Joint refinement requires at least two sessions')
    reference, mask, current = prepare(runs, grid, work)
    work = Path(work)
    log = work / 'commands.jsonl'
    for iteration in range(2):
        means = {s: current[[i for i, g in enumerate(groups) if g == s]].mean(axis=0) for s in sessions}
        targets = {}
        for s in sessions:
            targets[s] = work / f'target_{iteration}_{sessions.index(s)}.nii.gz'
            save(np.mean([a for group, a in means.items() if group != s], axis=0), reference, targets[s])
        raw = [fit(work / f'normalized_{i}.nii.gz', targets[r['session']],
                   work / 'weight.nii.gz', work / 'identity.mat',
                   work / f'raw_{iteration}_{i}.mat', log) for i, r in enumerate(runs)]
        corrections = center_corrections(raw, groups)
        current = []
        paths = []
        for i, (run, matrix) in enumerate(zip(runs, corrections)):
            path = work / f'correction_{iteration}_{i}.mat'
            np.savetxt(path, matrix)
            output = work / f'refined_{iteration}_{i}.nii.gz'
            render(run, grid, output, path, log)
            current.append(normalized(output, reference, mask))
            paths.append(str(path))
        current = np.stack(current)
    session_means = [current[[i for i, g in enumerate(groups) if g == s]].mean(axis=0) for s in sessions]
    save(np.mean(session_means, axis=0), reference, work / 'template.nii.gz')
    return dict(corrections=dict(zip([r['id'] for r in runs], paths)),
                template=str(work / 'template.nii.gz'), weight=str(work / 'weight.nii.gz'),
                normalization_mask=str(work / 'normalization_mask.nii.gz'))


def estimate_frozen(runs, grid, artifact, work):
    """Fit independently; never re-center new runs or refit the template pose."""
    artifact = Path(artifact)
    reference, _, _ = prepare(runs, grid, work)
    mask_image = nib.load(artifact / 'normalization_mask.nii.gz')
    if mask_image.shape != reference.shape or not np.allclose(mask_image.affine, reference.affine, atol=1e-5):
        raise ValueError('Frozen normalization mask grid mismatch')
    mask = mask_image.get_fdata() > 0
    work = Path(work)
    corrections = {}
    for i, run in enumerate(runs):
        # Reuse calibration mask, not a newly estimated incoming mask.
        data = normalized(work / f'baseline_{i}.nii.gz', reference, mask)
        save(data, reference, work / f'normalized_{i}.nii.gz')
        path = work / f'correction_{i}.mat'
        matrix = fit(work / f'normalized_{i}.nii.gz', artifact / 'template.nii.gz',
                     artifact / 'weight.nii.gz', work / 'identity.mat', path, work / 'commands.jsonl')
        np.savetxt(path, matrix)
        corrections[run['id']] = str(path)
    return dict(corrections=corrections)
