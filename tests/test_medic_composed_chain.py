import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
"""Known-answer tests of full FSL chain, including nonidentity motion."""
import os,subprocess,tempfile,unittest
from pathlib import Path
import numpy as np
import nibabel as nib
from medic_geometry import map_to_fsl
from medic_composition_reference import check_chain
os.environ['FSLOUTPUTTYPE']='NIFTI_GZ'

class ChainTests(unittest.TestCase):
    def test_affine_and_nonlinear_target_chain(self):
        with tempfile.TemporaryDirectory(prefix='medic-chain-',dir='/tmp') as folder:
            p=Path(folder);shape=(24,26,28);coords=np.indices(shape,dtype=np.float32)
            a=np.diag([-2.,2.,2.,1.]);a[:3,3]=[20,-20,-20]
            ref=nib.Nifti1Image(np.ones(shape,np.float32),a);nib.save(ref,p/'ref.nii.gz')
            d=nib.Nifti1Image(.8*np.sin(coords[1]/8)+.1*coords[2],a)
            nativewarp=map_to_fsl(d,ref,'j-');nib.save(nativewarp,p/'d.nii.gz')
            m=np.eye(4);angle=.07;c,s=np.cos(angle),np.sin(angle);m[:2,:2]=[[c,-s],[s,c]];m[:3,3]=[1.2,-.7,.4];np.savetxt(p/'m.mat',m)
            mask=np.zeros(shape,bool);mask[5:-5,5:-5,5:-5]=True
            for nonlinear in [False,True]:
                w=np.zeros((*shape,3),np.float32);w[...,0]=.7;w[...,2]=-.6
                if nonlinear:w[...,0]+=.3*np.sin(coords[2]/4);w[...,1]=.25*np.cos(coords[0]/5)
                base=nib.Nifti1Image(w,a);nib.save(base,p/'w.nii.gz')
                subprocess.run(['/usr/local/fsl/bin/convertwarp',f'--ref={p}/ref.nii.gz',f'--premat={p}/m.mat',f'--warp1={p}/w.nii.gz',f'--out={p}/mw.nii.gz','--rel','--relout'],check=True,capture_output=True)
                subprocess.run(['/usr/local/fsl/bin/convertwarp',f'--ref={p}/ref.nii.gz',f'--warp1={p}/d.nii.gz',f'--warp2={p}/mw.nii.gz',f'--out={p}/out.nii.gz','--rel','--relout'],check=True,capture_output=True)
                combined=nib.load(p/'out.nii.gz')
                result=check_chain(ref,ref,nativewarp,base,m,combined,mask)
                self.assertLess(result['max_mm'],.001)
                with self.assertRaises(ValueError):check_chain(ref,ref,nativewarp,base,np.linalg.inv(m),combined,mask)
                with self.assertRaises(ValueError):check_chain(ref,ref,nativewarp,base,np.eye(4),combined,mask)
                print('Known chain',nonlinear,result,flush=True)
if __name__=='__main__':unittest.main()
