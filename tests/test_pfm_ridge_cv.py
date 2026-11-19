import unittest,tempfile,json,sys,os
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
from pfm_ridge_cv_data import make_splits,plan_comparisons
from pfm_ridge_cv_selection import choose
from pfm_ridge_cv_metrics import prepare,mean_grid,aggregate_directions
from pfm_ridge_cv import positive_vector
class CVTests(unittest.TestCase):
 def test_partial_rotations_pair_directions_and_keep_valid_draws(self):
  def raw(effect):
   a=np.array(effect,float)
   return dict(effect=a,network_observed=np.repeat(a[:,:,None],20,axis=2),network_null=np.zeros(a.shape+(20,)))
  a=raw([[.2,.4],[.6,np.nan],[np.nan,.8]])
  b=raw([[.4,.6],[.8,.7],[.9,1.]])
  fc,net,n=aggregate_directions(a,b)
  np.testing.assert_allclose(fc,[.5,.7]);np.testing.assert_array_equal(n,[2,2]);np.testing.assert_allclose(net,np.repeat(fc[:,None],20,axis=1))
 def test_no_supported_rotation_remains_undefined(self):
  a=dict(effect=np.full((2,1),np.nan),network_observed=np.full((2,1,20),np.nan),network_null=np.zeros((2,1,20)))
  fc,net,n=aggregate_directions(a,a)
  self.assertTrue(np.isnan(fc).all());self.assertTrue(np.isnan(net).all());self.assertEqual(n[0],0)
 def test_complete_rotations_preserve_original_average(self):
  rng=np.random.default_rng(12)
  a=dict(effect=rng.normal(size=(5,3)),network_observed=rng.normal(size=(5,3,20)),network_null=rng.normal(size=(5,3,20)))
  b={k:rng.normal(size=v.shape) for k,v in a.items()}
  fc,net,n=aggregate_directions(a,b)
  np.testing.assert_allclose(fc,((a['effect']+b['effect'])/2).mean(0))
  np.testing.assert_allclose(net,((a['network_observed']-a['network_null']+b['network_observed']-b['network_null'])/2).mean(0));np.testing.assert_array_equal(n,[5,5,5])
 def split(self,runs=8,frames=100,tr=2):
  t=tempfile.TemporaryDirectory();self.addCleanup(t.cleanup);p=Path(t.name)/'scans.json';p.write_text(json.dumps([dict(n_timepoints=frames) for i in range(runs)]));return make_splits(p,np.zeros(runs*frames),runs*frames,tr,.3,4,3,.1)
 def test_temporal_crossing(self):
  split=self.split(2);pairs=plan_comparisons(split,.1,10);self.assertEqual(len(pairs),1);p=pairs[0]
  for side in ['A','B']:self.assertEqual(set(p[side]['groups']),{0,1})
  self.assertFalse(set(p['A']['rows'])&set(p['B']['rows']))
  self.assertLess(max(p['A']['raw_rows'][p['A']['groups']==0]),50)
  self.assertGreater(min(p['A']['raw_rows'][p['A']['groups']==1]),150)
 def test_even_stride_grouping_balances_time(self):
  split=self.split(8);pairs=plan_comparisons(split,.1,2);self.assertEqual(len(pairs),2)
  for p in pairs:
   early=[]
   for g in sorted(set(p['A']['groups'])):
    local=p['A']['raw_rows'][p['A']['groups']==g]-int(g)*100;early.append(np.max(local)<50)
   self.assertIn(True,early);self.assertIn(False,early)
  rows=np.concatenate([p[s]['rows'] for p in pairs for s in ['A','B']]);self.assertEqual(len(np.unique(rows)),len(rows))
 def test_short_partitions_fail(self):
  with self.assertRaises(ValueError):plan_comparisons(self.split(2),100,10)
 def test_retention_and_control(self):
  rows=[{}, {}, {},dict(control=True)];r=choose([.5,.46,.4,.9],[.6,.8,.9,1],rows,.1);self.assertEqual(r['index'],1);self.assertAlmostEqual(r['threshold'],.45)
  self.assertEqual(choose([.5,.46,.4,.9],[.6,.8,.9,1],rows,0)['index'],0)
 def test_bad_ranges_and_no_positive_fit(self):
  for v in ['0,.1','-1,1','nan','inf','']:
   with self.assertRaises(ValueError):positive_vector(v,'spatial')
  with self.assertRaises(ValueError):choose([-.1,-.2],[.9,.95],[{},{}])
 def test_exact_profile_shortcut(self):
  rng=np.random.default_rng(7);x=rng.normal(size=(60,15));y=rng.normal(size=(30,15));lab=np.repeat(np.arange(1,21),3)[None];score,net,count=mean_grid(prepare(x,y),lab);corr=np.corrcoef(x@y.T)
  expected=[]
  for k in range(1,21):
   ids=np.flatnonzero(lab[0]==k);c=corr[np.ix_(ids,ids)];expected.append((c.sum()-3)/6)
  np.testing.assert_allclose(net[0],expected,atol=1e-12);self.assertAlmostEqual(score[0],np.mean(expected))
 def test_zero_weight_missing_vs_singleton(self):
  rng=np.random.default_rng(2);metric=prepare(rng.normal(size=(60,15)),rng.normal(size=(30,15)));lab=np.ones((1,60),int);self.assertTrue(np.isfinite(mean_grid(metric,lab)[0][0]));lab[0,0]=2;self.assertTrue(np.isnan(mean_grid(metric,lab)[0][0]))
if __name__=='__main__':unittest.main()
