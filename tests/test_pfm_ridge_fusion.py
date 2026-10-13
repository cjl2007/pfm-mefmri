#!/usr/bin/env python3
"""Focused tests for RidgePFM scoring helpers."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

from pfm_ridge_fusion import clamp_prob_smooth, spatial_prior_term  # noqa: E402


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


if __name__ == "__main__":
    unittest.main()
