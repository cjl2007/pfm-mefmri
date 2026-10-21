import sys, unittest, tempfile
from pathlib import Path
import numpy as np
from scipy.io import savemat
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from pfm_ridge_confidence import assignment_states, profile_correlation, load_sessions, session_splits
from pfm_ridge_fusion import target_specific_ridge, solve_target_grams

class ConfidenceTests(unittest.TestCase):
    def test_native_local_terms_match_numpy(self):
        from pfm_ridge_local import local_terms
        rng=np.random.default_rng(7);z=rng.normal(size=(8,31));a=rng.normal(size=(8,3))
        excluded=[np.array([i,(i+1)%8],int) for i in range(8)]
        with tempfile.TemporaryDirectory() as folder:
            got=local_terms(z,a,excluded,folder)
        if got is None:self.skipTest('C++ compiler unavailable')
        for i,ids in enumerate(excluded):
            np.testing.assert_allclose(got[0][i],(z[ids]@z[i]/31)@a[ids],atol=1e-12)
            np.testing.assert_allclose(got[1][i],a[ids].T@a[ids],atol=1e-12)

    def test_four_states_and_missing(self):
        labels=np.ones((6,5),int);labels[1,[1,3]]=2;labels[:,4]=0
        correlations=np.tile([.8,.8,.5,.5,np.nan],(3,1))
        np.testing.assert_array_equal(assignment_states(labels,correlations),[1,2,3,4,0])
    def test_profiles_use_retained_features_and_reject_constant(self):
        x=np.array([[1.,2],[2,2],[4,2],[100,2]])
        y=x.copy();y[-1,0]=-100
        keep=np.ones(x.shape,bool);keep[-1]=False
        r=profile_correlation(x,y,keep)
        self.assertAlmostEqual(r[0],1);self.assertTrue(np.isnan(r[1]))
    def test_session_alignment_and_disjoint_splits(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'ScanIdx.mat';savemat(p,{'ScanIdx':[[1,1],[2,1],[3,1],[4,1]]})
            s=load_sessions(p,4)
            with self.assertRaises(ValueError):load_sessions(p,3)
            for a,b in session_splits(s):
                self.assertFalse(set(a)&set(b));self.assertEqual(set(a)|set(b),{1,2,3,4})
    def test_batched_fit_matches_explicit_retained_regression(self):
        rng=np.random.default_rng(9);a=rng.normal(size=(31,5));m=rng.normal(size=(31,7))
        keep=[rng.random(31)>.3 for _ in range(7)]
        for j,k in enumerate(keep):m[~k,j]=0
        grams=np.stack([a[k].T@a[k] for k in keep]);projection=m.T@a
        beta=solve_target_grams(projection,grams,1)
        ref,r2=target_specific_ridge(a,m,keep,1)
        np.testing.assert_allclose(beta,ref,atol=1e-12)
        ss=(m*m).sum(0);total=ss-m.sum(0)**2/np.array([k.sum() for k in keep])
        residual=ss-2*(beta*projection).sum(1)+np.einsum('bi,bij,bj->b',beta,grams,beta)
        np.testing.assert_allclose(1-residual/total,r2,atol=1e-7)
if __name__=='__main__':unittest.main()
