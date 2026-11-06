import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
import unittest
import numpy as np
import nibabel as nib
from medic_geometry import map_to_fsl,scaled_voxel_affine,frame_indices

class GeometryTests(unittest.TestCase):
    def test_world_coordinate_equivalence(self):
        rng=np.random.default_rng(71)
        for handedness in [-1,1]:
            for theta in [0,.2,-.45]:
                c,s=np.cos(theta),np.sin(theta)
                rot=np.array([[1,0,0],[0,c,-s],[0,s,c]])
                aff=np.eye(4);aff[:3,:3]=rot@np.diag([handedness*2.,2.4,3.]);aff[:3,3]=[19,-42,11]
                ref=nib.Nifti1Image(np.zeros((9,10,11),np.float32),aff)
                values=rng.normal(size=ref.shape).astype(np.float32)
                im=nib.Nifti1Image(values,aff)
                for axis,index in [('i',0),('j-',1),('k',2)]:
                    got=map_to_fsl(im,ref,axis).get_fdata()
                    physical=np.zeros((*ref.shape,3));physical[...,index]=values
                    # Independent coordinate path: add physical displacement
                    # to a world point, map back to source voxel, then FSL mm.
                    p=np.array([3.,4.,5.,1.]);world=aff@p
                    for q in [(1,2,3),(5,6,7)]:
                        shifted=world.copy();shifted[:3]+=physical[q]
                        vox=np.linalg.solve(aff,shifted)
                        expect=(scaled_voxel_affine(ref)@(vox-p))[:3]
                        np.testing.assert_allclose(got[q],expect,atol=2e-6)
    def test_map_sign_not_reapplied(self):
        im=nib.Nifti1Image(np.ones((5,5,5),np.float32),np.diag([-2,2,2,1]))
        np.testing.assert_array_equal(map_to_fsl(im,im,'j').get_fdata(),map_to_fsl(im,im,'j-').get_fdata())
    def test_frame_counts(self):
        self.assertEqual(frame_indices(5,5,2),[(2,2),(3,3),(4,4)])
        self.assertEqual(frame_indices(5,1,2),[(2,0),(3,0),(4,0)])
        for args in [(5,4,0),(5,5,5),(5,5,-1)]:
            with self.assertRaises(ValueError):frame_indices(*args)
    def test_reject_geometry_and_data(self):
        im=nib.Nifti1Image(np.zeros((5,5,5),np.float32),np.eye(4))
        bad=nib.Nifti1Image(np.zeros((6,5,5),np.float32),np.eye(4))
        with self.assertRaises(ValueError):map_to_fsl(bad,im,'j')
        with self.assertRaises(ValueError):map_to_fsl(im,im,'bad')
        with self.assertRaises(ValueError):map_to_fsl(im,im,'j',1)
        bad=nib.Nifti1Image(np.full((5,5,5),np.nan,np.float32),np.eye(4))
        with self.assertRaises(ValueError):map_to_fsl(bad,im,'j')
        aff=np.eye(4);aff[0,1]=.1;bad=nib.Nifti1Image(np.zeros((5,5,5),np.float32),aff)
        with self.assertRaises(ValueError):map_to_fsl(bad,bad,'j')

if __name__=='__main__':unittest.main()
