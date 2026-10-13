#!/usr/bin/env python3
"""MATLAB FD reference checks and concatenation/QC censoring regressions."""

from __future__ import annotations

import contextlib
import csv
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import nibabel as nib
from nibabel.cifti2.cifti2_axes import BrainModelAxis, SeriesAxis
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

from mefmri_motion import (  # noqa: E402
    FD_METHOD,
    compute_fd_from_motion,
    fd_settings,
    filtered_motion_for_fd,
)
from pfm_module import compute_fd_from_mcfpar  # noqa: E402
import mefmri_qc_motion_grayplots as qc  # noqa: E402


FIXTURE = Path(__file__).parent / "fixtures" / "historical_fd_matlab.npz"


class HistoricalFDTests(unittest.TestCase):
    def test_matches_historical_matlab_fd_and_censor_masks(self) -> None:
        with np.load(str(FIXTURE)) as reference:
            for tr, motion, expected in zip(reference["trs"], reference["motions"], reference["fd"]):
                with self.subTest(tr=tr):
                    original = motion.copy()
                    filtered, _ = filtered_motion_for_fd(motion, tr)
                    actual = compute_fd_from_motion(filtered, tr)
                    # The legacy degree-20 polynomial has small rounding errors;
                    # these references were run in MATLAB, not computed by Python.
                    np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-5)
                    np.testing.assert_array_equal(actual < 0.3, expected < 0.3)
                    np.testing.assert_array_equal(motion, original)

    def test_matlab_rounds_two_and_a_half_frames_up(self) -> None:
        motion = np.zeros((8, 6))
        motion[:, 3] = np.arange(8) * 0.1
        np.testing.assert_allclose(
            compute_fd_from_motion(motion, 1.0), [0, 0, 0, 0.3, 0.3, 0.3, 0.3, 0.3]
        )

    def test_rotation_displacements_use_fifty_mm_radius(self) -> None:
        motion = np.zeros((5, 6))
        motion[1:, :3] = [0.001, -0.002, 0.003]
        np.testing.assert_allclose(compute_fd_from_motion(motion, 2.0), [0, 0.3, 0, 0, 0])

    def test_two_second_and_slower_runs_skip_filter_even_when_short(self) -> None:
        motion = np.arange(30).reshape(5, 6)
        for tr in (2.0, 2.5):
            with self.subTest(tr=tr):
                filtered, label = filtered_motion_for_fd(motion, tr)
                np.testing.assert_array_equal(filtered, motion)
                self.assertEqual(label, "unfiltered_TR_ge_2s")
                self.assertFalse(fd_settings(tr)["filter_applied"])

    def test_filter_is_stable_just_below_two_seconds(self) -> None:
        # A constant rigid pose must not turn into apparent head movement.
        motion = np.tile([0.01, -0.02, 0.03, 1.0, -2.0, 3.0], (240, 1))
        filtered, _ = filtered_motion_for_fd(motion, 1.999)
        np.testing.assert_allclose(filtered, motion, rtol=0, atol=1e-10)
        np.testing.assert_allclose(compute_fd_from_motion(filtered, 1.999), 0, rtol=0, atol=1e-10)

    def test_invalid_or_short_motion_does_not_fall_back_to_raw_fd(self) -> None:
        for tr in (0.0, -1.0, np.nan, np.inf):
            with self.subTest(tr=tr), self.assertRaisesRegex(ValueError, "positive TR"):
                filtered_motion_for_fd(np.zeros((100, 6)), tr)
        with self.assertRaisesRegex(ValueError, "more than 60 frames"):
            filtered_motion_for_fd(np.zeros((60, 6)), 0.8)
        invalid = np.zeros((100, 6))
        invalid[40, 3] = np.nan
        with self.assertRaisesRegex(ValueError, "finite"):
            filtered_motion_for_fd(invalid, 0.8)
        with self.assertRaisesRegex(ValueError, "time-by-6"):
            filtered_motion_for_fd(np.zeros((100, 5)), 0.8)


class ConcatFDCensoringTests(unittest.TestCase):
    def test_saved_fd_cifti_qc_and_run_boundaries_agree(self) -> None:
        with tempfile.TemporaryDirectory(prefix="test-historical-fd-") as tmpdir:
            study = Path(tmpdir)
            subject_dir = study / "subject"
            brain = BrainModelAxis.from_surface(np.arange(3), 3, name="CortexLeft")
            with np.load(str(FIXTURE)) as reference:
                motion_filtered_run = reference["motions"][0].copy()
                fd_reference = reference["fd"][0].copy()
            # Exactly at, immediately below, and immediately above 0.3 mm.
            # Rounding FD.txt to six places would lose the second distinction.
            motion_boundary_run = np.zeros((4, 6))
            motion_boundary_run[:, 3] = [0, 0.3, 0.00000004, 0.30000008]
            runs = [(0.8, motion_filtered_run), (2.0, motion_boundary_run)]
            expected_blocks = []
            fds = []
            for run_number, (tr, motion) in enumerate(runs, start=1):
                run_dir = subject_dir / "func" / "rest" / "session_1" / f"run_{run_number}"
                run_dir.mkdir(parents=True)
                (run_dir / "TR.txt").write_text(str(tr))
                motion_file = run_dir / ("MCF.par" if run_number == 1 else "mcf.par")
                np.savetxt(str(motion_file), motion, fmt="%.17g")
                fds.append(compute_fd_from_mcfpar(motion_file, tr))
                frame = np.arange(len(motion), dtype=np.float32)
                data = np.column_stack([frame, frame ** 2, np.sin(frame)]).astype(np.float32)
                data += 100 * run_number
                header = nib.Cifti2Header.from_axes((SeriesAxis(0, tr, len(data)), brain))
                nib.save(nib.Cifti2Image(data, header), str(run_dir / "Rest_OCME+MEICA+MGTR.dtseries.nii"))
                expected_blocks.append(data - np.mean(data, axis=0, keepdims=True, dtype=np.float64).astype(np.float32))
            np.testing.assert_array_equal(fds[0] < 0.3, fd_reference < 0.3)
            np.testing.assert_array_equal(fds[1] < 0.3, [True, False, True, False])

            result = subprocess.run(
                [sys.executable, "-B", str(ROOT / "lib" / "pfm_module.py"),
                 "--subject", "subject", "--study-folder", str(study)],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            out = subject_dir / "func" / "rest" / "ConcatenatedCiftis"
            saved_fd = np.loadtxt(str(out / "FD.txt"))
            np.testing.assert_array_equal(saved_fd, np.concatenate(fds))
            keep = saved_fd < 0.3
            full = np.vstack(expected_blocks)
            for suffix, expected in (("Concatenated", full), ("Concatenated+FDlt0p3", full[keep])):
                img = nib.load(str(out / f"Rest_OCME+MEICA+MGTR_{suffix}.dtseries.nii"))
                np.testing.assert_array_equal(np.asanyarray(img.dataobj), expected)
            np.testing.assert_array_equal(np.loadtxt(str(out / "ScanIdx.txt")), [1] * 240 + [2] * 4)
            summary = json.loads((out / "PFMPrepSummary.json").read_text())
            self.assertEqual(summary["n_timepoints_kept"], int(np.sum(keep)))
            self.assertEqual(summary["fd_method"], FD_METHOD)
            info = json.loads((out / "ScanInfo.json").read_text())
            self.assertEqual([row["fd"]["filter_applied"] for row in info], [True, False])
            self.assertEqual([row["fd"]["lag_frames"] for row in info], [3, 1])

            # Exercise QC's actual report metrics, without spending time plotting.
            qc_out = study / "qc"
            with patch.object(qc, "plot_motion_grayplot_panel"), patch.object(qc, "plot_motion_coupling_summary"), patch.object(qc, "plot_qcfc_proxy_summary"), contextlib.redirect_stdout(io.StringIO()):
                qc.generate_report(subject_dir, "rest", "Rest", 1, qc_out, 3, 3, 0.3)
            with (qc_out / "motion_grayplot_metrics.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([int(row["n_retained"]) for row in rows], [int(np.sum(fd < 0.3)) for fd in fds])
            self.assertEqual(rows[1]["fd_filter"], "unfiltered_TR_ge_2s")
            self.assertEqual(json.loads((qc_out / "manifest.json").read_text())["fd_method"], FD_METHOD)


if __name__ == "__main__":
    unittest.main()
