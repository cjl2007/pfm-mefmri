import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
import unittest,tempfile
from pathlib import Path
import numpy as np
import nibabel as nib
from scipy.ndimage import map_coordinates
from joint_sampler import fsl_temporal_weights,sample_joint_frame

class JointTests(unittest.TestCase):
    def test_each_contributor_uses_own_motion(self):
        shape=(12,26,8,40);x,y,z=np.indices(shape[:3],dtype=float)
        with tempfile.TemporaryDirectory(prefix='medic-joint-') as folder:
            weights=fsl_temporal_weights(40,np.full(8,.4),folder)
            shifts=2*np.sin(np.arange(40)*2*np.pi*.25)
            raw=np.stack([80+30*np.exp(-(y-d-13)**2/12) for d in shifts],axis=3).astype(np.float32)
            original=np.array([x.ravel(),y.ravel(),z.ravel()]);mask=(y.ravel()>7)&(y.ravel()<19)
            def pull(k):return original+np.array([0,shifts[k],0])[:,None]
            truth=80+30*np.exp(-(y.ravel()-13)**2/12)
            errors=[];bad=[]
            for t in [10,11,12]:
                got=sample_joint_frame(raw,np.zeros(shape[:3]),np.eye(4),'j',pull,weights,t)
                control=sample_joint_frame(raw,np.zeros(shape[:3]),np.eye(4),'j',pull(t),weights,t)
                errors.append(float(np.sqrt(np.mean((got[mask]-truth[mask])**2))/100))
                bad.append(float(np.sqrt(np.mean((control[mask]-truth[mask])**2))/100))
            self.assertLess(max(errors),.001);self.assertGreater(min(bad),.005)
            print('Own-frame motion RMS:',errors,'shared-target-motion control:',bad,flush=True)
    def test_fractional_slice_temporal_accuracy(self):
        shape=(8,9,12,60);x,y,z=np.indices(shape[:3],dtype=float)
        acq=np.tile([0,.5,.1,.6,.2,.7],2);offset=.5-acq
        with tempfile.TemporaryDirectory(prefix='medic-joint-') as folder:
            weights=fsl_temporal_weights(shape[3],offset,folder)
            raw=np.stack([100+10*np.sin(2*np.pi*.12*(t+acq[z.astype(int)])) for t in range(shape[3])],axis=3).astype(np.float32)
            q=np.array([x.ravel(),y.ravel(),z.ravel()+.45]);mask=(z.ravel()>2)&(z.ravel()<8)
            errors=[]
            for t in [12,18,25,31]:
                actual=sample_joint_frame(raw,np.zeros(shape[:3]),np.eye(4),'j',q,weights,t)
                truth=100+10*np.sin(2*np.pi*.12*(t+.5))
                errors.append(float(np.sqrt(np.mean((actual[mask]-truth)**2))/100))
            print('Fractional slice temporal RMS:',errors,flush=True)
            self.assertLess(max(errors),.001)
            # With no distortion, temporal correction per native slice and
            # the single spatial kernel must commute exactly by linearity.
            corrected=np.einsum('xyzk,ztk->xyzt',raw,weights)
            identity=np.broadcast_to(np.eye(shape[3],dtype=np.float32),(shape[2],shape[3],shape[3]))
            for t in [0,18,59]:
                got=sample_joint_frame(raw,np.zeros(shape[:3]),np.eye(4),'j',q,weights,t)
                reference=sample_joint_frame(corrected,np.zeros(shape[:3]),np.eye(4),'j',q,identity,t)
                np.testing.assert_allclose(got,reference,atol=3e-5)
    def test_dynamic_field_and_temporal_boundaries(self):
        shape=(12,24,8,40);x,y,z=np.indices(shape[:3],dtype=float);aff=np.eye(4)
        acq=np.array([0,.5,.125,.625,.25,.75,.375,.875]);offset=.5-acq
        with tempfile.TemporaryDirectory(prefix='medic-joint-') as folder:
            weights=fsl_temporal_weights(shape[3],offset,folder)
            np.testing.assert_allclose(weights.sum(axis=2),1,atol=2e-6)
            raw=np.empty(shape,np.float32);maps=np.empty(shape,np.float32)
            for t in range(shape[3]):
                tt=t+acq[z.astype(int)];d=.4+3*np.sin(2*np.pi*.35*tt)
                raw[...,t]=80+30*np.exp(-((y-d-12)**2)/12)
                maps[...,t]=d
            q=np.array([x.ravel(),y.ravel(),z.ravel()]);mask=((y>6)&(y<18)&(x>2)&(x<9))
            expected=80+30*np.exp(-(y-12)**2/12)
            errors=[]
            for t in [0,1,10,20,38,39]:
                actual=sample_joint_frame(raw,maps,aff,'j-',q,weights,t).reshape(shape[:3])
                error=float(np.sqrt(np.mean((actual[mask]-expected[mask])**2))/np.sqrt(np.mean(expected[mask]**2)))
                errors.append(error)
                self.assertLess(error,.001)
            print('Rapid-field joint RMS, including temporal edges:',errors,flush=True)
    def test_static_zero_offset_and_invalid_counts(self):
        rng=np.random.default_rng(7);raw=rng.normal(size=(9,10,8,5)).astype(np.float32)
        with tempfile.TemporaryDirectory(prefix='medic-joint-') as folder:
            weights=fsl_temporal_weights(5,np.zeros(8),folder)
            q=np.indices(raw.shape[:3],dtype=float).reshape(3,-1)
            got=sample_joint_frame(raw,np.zeros(raw.shape[:3]),np.eye(4),'j',q,weights,2).reshape(raw.shape[:3])
            np.testing.assert_allclose(got,raw[...,2],atol=1e-6)
            with self.assertRaises(ValueError):sample_joint_frame(raw,np.zeros((*raw.shape[:3],4)),np.eye(4),'j',q,weights,2)
if __name__=='__main__':unittest.main()
