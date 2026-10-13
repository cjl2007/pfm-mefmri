#!/usr/bin/env python3
"""Focused tests for the parcel ambiguity stripe helper."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

from pfm_parcel_uncertainty_stripes import (  # noqa: E402
    build_subcortical_highres_volume,
    highres_affine_and_factors,
    load_parcel_support_csv,
    multi_stripe_assignments,
    multi_voxel_stripe_assignments,
    normalized_entropy,
    parcel_network_from_label,
    stripe_assignments,
)


class ParcelUncertaintyStripeTests(unittest.TestCase):
    def test_network_name_is_recovered_from_areal_label(self) -> None:
        names = ["Default_Parietal", "Visual_V1"]
        self.assertEqual(parcel_network_from_label(17, "Default_Parietal_17", names), (0, "parcel_label_table"))

    def test_unknown_label_requests_probability_fallback(self) -> None:
        idx, source = parcel_network_from_label(2, "Unknown_2", ["Known"])
        self.assertEqual(idx, -1)
        self.assertEqual(source, "probability_argmax_fallback")

    def test_entropy_endpoints(self) -> None:
        self.assertAlmostEqual(normalized_entropy(np.array([1.0, 0.0, 0.0])), 0.0)
        self.assertAlmostEqual(normalized_entropy(np.full(4, 0.25)), 1.0)

    def test_stripes_contain_both_networks(self) -> None:
        latitude = np.linspace(-0.1, 0.1, 50)
        values = stripe_assignments(latitude, 0.025, primary_id=3, alternative_id=8)
        self.assertEqual(set(np.unique(values)), {3, 8})

    def test_narrow_parcel_fallback_contains_both_networks(self) -> None:
        latitude = np.linspace(0.001, 0.002, 20)
        values = stripe_assignments(latitude, 1.0, primary_id=1, alternative_id=2)
        self.assertEqual(set(np.unique(values)), {1, 2})
        self.assertEqual(int(np.sum(values == 1)), 10)
        self.assertEqual(int(np.sum(values == 2)), 10)

    def test_three_network_surface_stripes_contain_every_network(self) -> None:
        latitude = np.linspace(-0.1, 0.1, 60)
        values = multi_stripe_assignments(latitude, 0.025, [2, 5, 9])
        self.assertEqual(set(np.unique(values)), {2, 5, 9})

    def test_three_network_voxel_stripes_contain_every_network(self) -> None:
        voxels = np.column_stack((np.arange(12), np.zeros(12), np.zeros(12))).astype(int)
        values = multi_voxel_stripe_assignments(voxels, 1, "axial", [1, 4, 7])
        self.assertEqual(set(np.unique(values)), {1, 4, 7})

    def test_highres_subcortical_block_contains_three_networks(self) -> None:
        affine = np.diag([2.0, 2.0, 2.0, 1.0])
        stripe_maps = np.array([[1, 2], [0, 3], [0, 4]], dtype=np.int32)
        count_map = np.array([1, 3], dtype=np.int32)
        voxels = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.int32)
        mixed = np.array([False, True])
        volume, highres_affine = build_subcortical_highres_volume(
            (2, 1, 1),
            affine,
            voxels,
            stripe_maps,
            count_map,
            mixed,
            output_resolution_mm=0.5,
            stripe_width_mm=0.5,
            view="axial",
        )
        expected_affine, factors = highres_affine_and_factors(affine, 0.5)
        self.assertEqual(volume.shape, (8, 4, 4))
        self.assertTrue(np.allclose(highres_affine, expected_affine))
        self.assertEqual(set(np.unique(volume[: factors[0]])), {1})
        self.assertEqual(set(np.unique(volume[factors[0] :])), {2, 3, 4})

    def test_cv_parcel_support_csv_is_loaded(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "support.csv"
            path.write_text(
                "parcel_id,selected_network_count,selected_network_ids\n"
                "1,1,4\n"
                "2,2,5;9\n"
                "3,3,1;2;7\n"
            )
            self.assertEqual(load_parcel_support_csv(path), {2: [5, 9], 3: [1, 2, 7]})


if __name__ == "__main__":
    unittest.main()
