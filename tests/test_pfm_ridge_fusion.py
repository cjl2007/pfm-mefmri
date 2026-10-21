#!/usr/bin/env python3
"""Focused tests for RidgePFM scoring helpers."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import nibabel as nib
import numpy as np
import scipy.io as sio
from nibabel.cifti2 import BrainModelAxis, Cifti2Header, SeriesAxis


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

from pfm_ridge_fusion import (  # noqa: E402
    clamp_prob_smooth,
    main,
    solve_target_grams,
    spatial_prior_term,
    target_specific_ridge,
)


class RidgePFMScoringTests(unittest.TestCase):
    def test_logit_prior_transform_matches_probability_log_odds(self) -> None:
        p = np.asarray([[0.25, 0.5, 0.75]], dtype=np.float32)
        got = spatial_prior_term(p, "logit")
        expected_p = clamp_prob_smooth(p)
        expected = np.log(expected_p / (1.0 - expected_p))
        self.assertTrue(np.allclose(got, expected))
        self.assertLess(got[0, 0], 0)
        self.assertAlmostEqual(float(got[0, 1]), 0.0, places=4)
        self.assertGreater(got[0, 2], 0)

    def test_log_prior_transform_is_retained_for_comparisons(self) -> None:
        p = np.asarray([[0.25, 0.5, 0.75]], dtype=np.float32)
        got = spatial_prior_term(p, "log")
        self.assertTrue(np.allclose(got, np.log(p)))

    def test_prior_smoothing_prevents_impossible_networks(self) -> None:
        p = np.asarray([[0.0, 1.0]], dtype=np.float32)
        got = clamp_prob_smooth(p)
        self.assertGreater(float(got[0, 0]), 0.0)
        self.assertLess(float(got[0, 1]), 1.0)

    def test_target_specific_ridge_matches_direct_masked_solve(self) -> None:
        design = np.asarray(
            [[1.0, 0.0], [0.5, 1.0], [-1.0, 0.5], [0.25, -0.75]],
            dtype=np.float64,
        )
        fingerprints = np.asarray(
            [[0.2, 0.0], [0.4, -0.3], [0.0, 0.6], [-0.1, 0.2]],
            dtype=np.float64,
        )
        retained = [
            np.asarray([True, True, False, True]),
            np.asarray([False, True, True, True]),
        ]
        got, got_r2 = target_specific_ridge(design, fingerprints, retained, 1.0)
        for target, keep in enumerate(retained):
            x = design[keep]
            y = fingerprints[keep, target]
            expected = np.linalg.solve(x.T @ x + np.eye(2), x.T @ y)
            residual = y - x @ expected
            expected_r2 = 1.0 - (residual @ residual) / max(
                ((y - y.mean()) ** 2).sum(), 1e-8
            )
            self.assertTrue(np.allclose(got[target], expected))
            self.assertAlmostEqual(float(got_r2[target]), float(expected_r2), places=6)

    def test_target_specific_ridge_ignores_excluded_values(self) -> None:
        design = np.asarray([[1.0], [2.0], [3.0]], dtype=np.float64)
        retained = [np.asarray([True, False, True])]
        first = np.asarray([[0.5], [999.0], [-0.25]], dtype=np.float64)
        second = np.asarray([[0.5], [-999.0], [-0.25]], dtype=np.float64)
        beta_a, r2_a = target_specific_ridge(design, first, retained, 1.0)
        beta_b, r2_b = target_specific_ridge(design, second, retained, 1.0)
        self.assertTrue(np.array_equal(beta_a, beta_b))
        self.assertTrue(np.array_equal(r2_a, r2_b))

    def test_projection_fast_path_matches_explicit_fc(self) -> None:
        rng = np.random.default_rng(11)
        features, targets, frames, networks = 13, 5, 31, 3
        x = rng.normal(size=(features, frames))
        x = (x - x.mean(1, keepdims=True)) / x.std(1, keepdims=True)
        y = rng.normal(size=(targets, frames))
        y = (y - y.mean(1, keepdims=True)) / y.std(1, keepdims=True)
        design = rng.normal(size=(features, networks))
        design -= design.mean(0, keepdims=True)
        design /= np.linalg.norm(design, axis=0, keepdims=True)
        excluded = [np.asarray([i, (i + 3) % features]) for i in range(targets)]
        fingerprints = x @ y.T / frames
        retained = []
        grams = np.broadcast_to(design.T @ design, (targets, networks, networks)).copy()
        projection = y @ (x.T @ design) / frames
        for target, local in enumerate(excluded):
            keep = np.ones(features, dtype=bool)
            keep[local] = False
            retained.append(keep)
            fingerprints[~keep, target] = 0.0
            grams[target] -= design[local].T @ design[local]
            local_fc = x[local] @ y[target] / frames
            projection[target] -= local_fc @ design[local]
        direct, _ = target_specific_ridge(design, fingerprints, retained, 1.0)
        fast = solve_target_grams(projection, grams, 1.0)
        self.assertTrue(np.allclose(fast, direct, atol=1e-12, rtol=0))

    def test_main_maps_volume_without_volume_spatial_priors(self) -> None:
        rng = np.random.default_rng(7)
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            cortex = BrainModelAxis.from_mask(
                np.ones(4, dtype=bool), name="CortexLeft"
            )
            volume_mask = np.zeros((2, 2, 2), dtype=bool)
            volume_mask[0, 0, 0] = True
            volume = BrainModelAxis.from_mask(
                volume_mask, affine=np.eye(4), name="ThalamusLeft"
            )
            brain = cortex + volume
            series = SeriesAxis(start=0.0, step=1.0, size=24)
            data = rng.normal(size=(24, brain.size)).astype(np.float32)
            cifti = tmp / "input.dtseries.nii"
            nib.save(nib.Cifti2Image(data, Cifti2Header.from_axes((series, brain))), cifti)
            distance = np.full((brain.size, brain.size), 50, dtype=np.uint8)
            np.fill_diagonal(distance, 0)
            distance_path = tmp / "distance.npy"
            np.save(distance_path, distance)
            priors = tmp / "priors.mat"
            sio.savemat(
                priors,
                {
                    "Priors": {
                        "FC": rng.normal(size=(4, 3)),
                        "Spatial": np.asarray(
                            [[0.8, 0.1, 0.1], [0.2, 0.7, 0.1],
                             [0.1, 0.2, 0.7], [0.6, 0.2, 0.2]],
                            dtype=np.float32,
                        ),
                        "NetworkLabels": np.asarray(["A", "B", "C"], dtype=object),
                        "NetworkColors": np.eye(3, dtype=np.float32),
                    }
                },
            )
            out = tmp / "out"
            argv = [
                "pfm_ridge_fusion.py", "--in-cifti", str(cifti),
                "--distance-npy", str(distance_path), "--priors-mat", str(priors),
                "--outdir", str(out), "--outfile", "Synthetic",
                "--left-surf", "unused.L.surf.gii", "--right-surf", "unused.R.surf.gii",
            ]
            with mock.patch.object(sys, "argv", argv), mock.patch(
                "pfm_ridge_fusion.subprocess.run"
            ) as run:
                self.assertEqual(main(), 0)
            self.assertEqual(run.call_count, 3)
            availability = np.asarray(
                nib.load(out / "Synthetic_SpatialPriorAvailable.dtseries.nii").dataobj
            )[0]
            self.assertTrue(np.array_equal(availability[:4], np.ones(4)))
            self.assertEqual(float(availability[4]), 0.0)
            softmax = np.asarray(
                nib.load(out / "Synthetic_FusedScoreSoftmax.dtseries.nii").dataobj
            )
            self.assertEqual(softmax.shape, (3, 5))
            self.assertTrue(np.allclose(softmax.sum(axis=0), 1.0, atol=1e-6))
            # Neutral volume priors contribute the same constant to all labels,
            # so the volume result remains an FC-only classification.
            self.assertTrue(np.isfinite(softmax[:, 4]).all())


if __name__ == "__main__":
    unittest.main()
