#!/usr/bin/env python3
"""Focused tests for subcortical regression neighborhood masks."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

from pfm_subcort_regress import compute_regression_masks_from_distance, regress_runs, load_coords_and_masks  # noqa: E402


class SubcorticalRegressionMaskTests(unittest.TestCase):
    def test_distinct_run_slopes_are_removed_separately(self):
        x=np.tile(np.arange(6,dtype=float),2);y=np.r_[2*x[:6]+3,-4*x[6:]+10]
        residual=regress_runs(y,x,[np.arange(6),np.arange(6,12)])
        np.testing.assert_allclose(residual,0,atol=1e-12)
        self.assertGreater(np.std(regress_runs(y,x,[np.arange(12)])),1)

    def test_volume_coordinates_use_zero_based_indices(self):
        import tempfile
        import nibabel as nib
        from nibabel.cifti2 import BrainModelAxis,SeriesAxis
        mask=np.zeros((3,4,5),bool);mask[1,2,3]=True
        affine=np.diag([-2.,2.,2.,1.]);affine[:3,3]=[90,-126,-72]
        axis=BrainModelAxis.from_surface([0],1,'CortexLeft')+BrainModelAxis.from_mask(mask,affine=affine,name='ThalamusLeft')
        img=nib.Cifti2Image(np.zeros((3,2)),nib.Cifti2Header.from_axes((SeriesAxis(0,1,3),axis)))
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'surface.surf.gii'
            nib.save(nib.gifti.GiftiImage(darrays=[nib.gifti.GiftiDataArray(np.zeros((1,3),np.float32),intent='NIFTI_INTENT_POINTSET')]),path)
            coords,cortex,volume=load_coords_and_masks(img,path,path)
        np.testing.assert_array_equal(coords[1],[88,-122,-66])
        self.assertTrue(volume[1]);self.assertTrue(cortex[0])

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
