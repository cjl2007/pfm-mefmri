#!/usr/bin/env python3
"""Focused tests for subcortical regression neighborhood masks."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

from pfm_subcort_regress import compute_regression_masks_from_distance  # noqa: E402


class SubcorticalRegressionMaskTests(unittest.TestCase):
    def test_cortex_mode_uses_only_nearby_cortical_rows(self) -> None:
        cortex = np.asarray([True, True, True, False, False])
        distance_row = np.asarray([5, 25, 10, 0, 8])
        should_regress, mask = compute_regression_masks_from_distance(
            row=3,
            distance_row=distance_row,
            cortex=cortex,
            distance_mm=20,
        )
        self.assertTrue(should_regress)
        self.assertEqual(mask.tolist(), [True, False, True, False, False])

    def test_deep_subcortical_rows_are_skipped(self) -> None:
        cortex = np.asarray([True, True, True, False, False])
        distance_row = np.asarray([25, 30, 40, 0, 8])
        should_regress, mask = compute_regression_masks_from_distance(
            row=3,
            distance_row=distance_row,
            cortex=cortex,
            distance_mm=20,
        )
        self.assertFalse(should_regress)
        self.assertFalse(mask.any())


if __name__ == "__main__":
    unittest.main()
