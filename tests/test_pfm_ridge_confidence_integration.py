import sys,tempfile,subprocess,json
from pathlib import Path
import nibabel as nib,numpy as np
from scipy.io import savemat
from nibabel.cifti2 import BrainModelAxis,SeriesAxis

import unittest,shutil

@unittest.skipUnless(shutil.which("wb_command"), "Workbench required")
class ConfidenceIntegrationTest(unittest.TestCase):
    def test_real_workbench_session_confidence(self):
        with tempfile.TemporaryDirectory() as folder:
         p=Path(folder);rng=np.random.default_rng(3);axis=BrainModelAxis.from_surface(np.arange(4),4,name='CortexLeft')+BrainModelAxis.from_surface(np.arange(4),4,name='CortexRight');data=rng.normal(size=(80,8)).astype(np.float32)
         nib.save(nib.Cifti2Image(data,nib.Cifti2Header.from_axes((SeriesAxis(0,1,80),axis))),p/'input.dtseries.nii')
         savemat(p/'scan.mat',{'ScanIdx':np.column_stack((np.repeat(np.arange(1,5),20),np.ones(80)))})
         pri={'FC':rng.normal(size=(8,2)), 'Spatial':rng.uniform(.1,.9,size=(8,2)), 'NetworkLabels':np.array(['A','B'],dtype=object),'NetworkColors':np.array([[1,0,0],[0,0,1]])};savemat(p/'priors.mat',{'Priors':pri})
         dist=np.full((8,8),255,np.uint8);np.fill_diagonal(dist,0);np.save(p/'dist.npy',dist)
         coords=np.array([[0,0,0],[2,0,0],[0,2,0],[0,0,2]],np.float32);faces=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]],np.int32)
         for h in ['L','R']:
          g=nib.gifti.GiftiImage(darrays=[nib.gifti.GiftiDataArray(coords,intent='NIFTI_INTENT_POINTSET'),nib.gifti.GiftiDataArray(faces,intent='NIFTI_INTENT_TRIANGLE')]);g.meta['AnatomicalStructurePrimary']='CortexLeft' if h=='L' else 'CortexRight';nib.save(g,p/(h+'.surf.gii'));subprocess.run(['wb_command','-set-structure',str(p/(h+'.surf.gii')),'CORTEX_LEFT' if h=='L' else 'CORTEX_RIGHT'],check=True)
         base=['--in-cifti',str(p/'input.dtseries.nii'),'--distance-npy',str(p/'dist.npy'),'--priors-mat',str(p/'priors.mat'),'--left-surf',str(p/'L.surf.gii'),'--right-surf',str(p/'R.surf.gii')]
         subprocess.run([sys.executable,'lib/pfm_ridge_fusion.py',*base,'--outdir',str(p/'full')],check=True)
         subprocess.run([sys.executable,'lib/pfm_ridge_confidence.py',*base,'--outdir',str(p/'conf'),'--scan-idx',str(p/'scan.mat'),'--minimum-frames','10','--full-map-dlabel',str(p/'full/RidgePFM_VTX.dlabel.nii')],check=True)
         f=nib.load(p/'conf/RidgePFM_VTX_SessionAssignmentFrequency.dscalar.nii').get_fdata();l=nib.load(p/'conf/RidgePFM_VTX_SessionAssignments.dlabel.nii').get_fdata()
         np.testing.assert_allclose(f.sum(0),1)
         for k in range(2):np.testing.assert_array_equal(f[k],(l==k+1).mean(0))
         states=nib.load(p/'conf/RidgePFM_VTX_AssignmentStates.dlabel.nii').get_fdata()[0];stripe=nib.load(p/'conf/RidgePFM_VTX_NetworkStripes.dlabel.nii').get_fdata()[0];full=nib.load(p/'full/RidgePFM_VTX.dlabel.nii').get_fdata()[0]
         np.testing.assert_array_equal(stripe[states!=2],full[states!=2]);print('PASS: real Workbench end-to-end session frequencies, four states, separate stripes')
