"""Rigid subject refinement in FSL scaled-mm coordinates.

These functions compose transforms only; they never resample BOLD intensities.
"""
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation


def rigid(matrix, tolerance=1e-4):
    matrix = np.asarray(matrix, dtype=float)
    if matrix.shape != (4, 4) or not np.isfinite(matrix).all():
        raise ValueError('Rigid transform must be a finite 4x4 matrix')
    if not np.allclose(matrix[3], [0, 0, 0, 1], atol=tolerance):
        raise ValueError('Invalid homogeneous transform row')
    rot = matrix[:3, :3]
    if not np.allclose(rot.T @ rot, np.eye(3), atol=tolerance) or not np.isclose(np.linalg.det(rot), 1, atol=tolerance):
        raise ValueError('Refinement must be rigid, with no scale, shear or reflection')
    return matrix.copy()


def center_corrections(matrices, sessions):
    """Equal session weighting, matching the validated two-iteration pilot."""
    if len(matrices) != len(sessions) or not matrices:
        raise ValueError('One session label is required per correction')
    if len(set(sessions)) < 2:
        raise ValueError('Joint refinement requires at least two sessions')
    mats = np.stack([rigid(m) for m in matrices])
    weights = np.array([1 / sessions.count(s) for s in sessions], dtype=float)
    weights /= weights.sum()
    mean = np.eye(4)
    mean[:3, :3] = Rotation.from_matrix(mats[:, :3, :3]).mean(weights=weights).as_matrix()
    mean[:3, 3] = np.sum(mats[:, :3, 3] * weights[:, None], axis=0)
    return [rigid(np.linalg.inv(mean) @ m) for m in mats]


def refined_affine(baseline, correction):
    """Intermediate -> ACPC. Per-volume motion is appended by the consumer."""
    baseline = np.asarray(baseline, dtype=float)
    if baseline.shape != (4, 4) or not np.isfinite(baseline).all():
        raise ValueError('Baseline affine must be finite and 4x4')
    return rigid(correction) @ baseline


def warp_commands(baseline, correction, acpc_reference, acpc_output,
                  standard_warp=None, standard_reference=None, standard_output=None):
    """Build ACPC correction first, then optional nonlinear standard mapping.

    The baseline is an intermediate-to-ACPC warp, NEVER an already-standard
    warp. No SBRef prematrix or BOLD motion matrix belongs in this shared base.
    FSL reads field conventions from the input headers; output is explicit relative.
    """
    optional = [standard_warp, standard_reference, standard_output]
    if any(x is not None for x in optional) and not all(x is not None for x in optional):
        raise ValueError('Standard warp, reference and output must be supplied together')
    paths = [baseline, correction, acpc_reference, acpc_output] + [x for x in optional if x is not None]
    if any(not str(p).strip() for p in paths):
        raise ValueError('Empty transform path')
    inputs = {Path(p).resolve() for p in [baseline, correction, acpc_reference] + [x for x in optional[:2] if x is not None]}
    outputs = [Path(acpc_output).resolve()] + ([Path(standard_output).resolve()] if standard_output else [])
    if len(set(outputs)) != len(outputs) or inputs.intersection(outputs):
        raise ValueError('Composed outputs must not overwrite inputs or each other')
    commands = [['convertwarp', '--ref=' + str(acpc_reference), '--warp1=' + str(baseline),
                 '--postmat=' + str(correction), '--out=' + str(acpc_output), '--relout']]
    if standard_warp is not None:
        commands.append(['convertwarp', '--ref=' + str(standard_reference),
                         '--warp1=' + str(acpc_output), '--warp2=' + str(standard_warp),
                         '--out=' + str(standard_output), '--relout'])
    return commands
