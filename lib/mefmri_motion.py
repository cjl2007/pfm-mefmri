"""Shared FD rules from the historical concatenate_scans.m calc_fd helper."""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from scipy.signal import butter, sosfiltfilt


FD_METHOD = "historical_concatenate_scans_sos_v1"
FD_FILTER_ORDER = 10
FD_FILTER_PADLEN = 6 * FD_FILTER_ORDER  # MATLAB: 3 * (21 coefficients - 1).


def fd_settings(tr: float) -> dict:
    """Describe the per-run calculation, including MATLAB's positive rounding."""
    if not np.isfinite(tr) or tr <= 0:
        raise ValueError(f"FD requires a finite, positive TR, got {tr}")
    filtered = tr < 2.0
    return {
        "method": FD_METHOD,
        "tr_seconds": float(tr),
        "lag_frames": int(np.floor(2.5 / tr + 0.5)),
        "rotation_radius_mm": 50.0,
        "filter_applied": filtered,
        "filter_type": "butterworth_bandstop" if filtered else "none",
        "filter_order": FD_FILTER_ORDER if filtered else None,
        "stopband_hz": [0.2, 0.5 / tr - 0.019] if filtered else None,
        "filter_implementation": "sosfiltfilt" if filtered else None,
        "filter_padtype": "odd" if filtered else None,
        "filter_padlen": FD_FILTER_PADLEN if filtered else 0,
    }


def _motion_array(motion: np.ndarray) -> np.ndarray:
    motion = np.asarray(motion, dtype=np.float64)
    if motion.ndim != 2 or motion.shape[1] != 6 or motion.shape[0] == 0:
        raise ValueError(f"Expected nonempty time-by-6 motion parameters, got {motion.shape}")
    if not np.all(np.isfinite(motion)):
        raise ValueError("FD motion parameters must all be finite")
    return motion


def filtered_motion_for_fd(
    motion: Optional[np.ndarray], tr: float
) -> Tuple[Optional[np.ndarray], str]:
    """Filter motion within a run when TR < 2 s; leave slower runs unfiltered.

    The historical 10th-order Butterworth bandstop spans 0.2 Hz through
    Nyquist - 0.019 Hz. Second-order sections preserve that transfer function
    without the numerical instability of its degree-20 direct-form polynomial.
    Odd reflection with 60 padding samples matches MATLAB filtfilt's edges.
    Like MATLAB, short/nonfinite inputs fail instead of silently using raw FD.
    """
    if motion is None:
        return None, "missing_motion"
    settings = fd_settings(tr)
    motion = _motion_array(motion)
    if not settings["filter_applied"]:
        return motion.copy(), "unfiltered_TR_ge_2s"
    if motion.shape[0] <= FD_FILTER_PADLEN:
        raise ValueError(
            f"Historical FD filtering at TR={tr:g}s requires more than "
            f"{FD_FILTER_PADLEN} frames, got {motion.shape[0]}"
        )
    nyq = 0.5 / tr
    stopband = settings["stopband_hz"]
    sos = butter(FD_FILTER_ORDER, np.asarray(stopband) / nyq, btype="bandstop", output="sos")
    filtered = sosfiltfilt(sos, motion, axis=0, padtype="odd", padlen=FD_FILTER_PADLEN)
    if not np.all(np.isfinite(filtered)):
        raise ValueError(f"Historical FD filtering produced nonfinite motion at TR={tr:g}s")
    return filtered, f"{stopband[0]:.3f}-{stopband[1]:.3f}Hz"


def compute_fd_from_motion(motion: Optional[np.ndarray], tr: float) -> Optional[np.ndarray]:
    """Calculate lagged FD from supplied (already filtered, if needed) motion."""
    if motion is None:
        return None
    lag = fd_settings(tr)["lag_frames"]
    motion = _motion_array(motion)
    delta = np.zeros_like(motion)
    # Leading frames are zero separately for each run, as in concatenate_scans.m.
    if 0 < lag < motion.shape[0]:
        delta[lag:] = np.abs(motion[lag:] - motion[:-lag])
    angular = delta[:, :3] / (2.0 * np.pi) * 100.0 * np.pi
    return np.sum(np.hstack([delta[:, 3:], angular]), axis=1)
