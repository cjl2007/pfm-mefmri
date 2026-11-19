import unittest,copy,tempfile
from pathlib import Path
from types import SimpleNamespace
from pfm_ridge_cv_dispatch import equivalent_grid
from pfm_ridge_cv_prep_cache import check,save
class CacheTests(unittest.TestCase):
 def test_only_policy_and_raw_vector_format_are_ignored(self):
  a=dict(settings=dict(max_fc_loss=.1,lambdas='.1,1',spatial_weights='.01,.1',lambda_values=[.1,1],weight_values=[.01,.1],n_rotations=100),inputs=[dict(path='x',mtime_ns=1)],code=dict(engine='same'))
  b=copy.deepcopy(a);b['settings']['max_fc_loss']=.2;b['settings']['lambdas']='0.1,1.0';self.assertTrue(equivalent_grid(a,b))
 def test_numerical_inputs_and_code_invalidate(self):
  a=dict(settings=dict(max_fc_loss=.1,lambda_values=[1],n_rotations=100),inputs=[dict(path='x',mtime_ns=1)],code=dict(engine='same'))
  for field in ['rotations','input','code']:
   b=copy.deepcopy(a)
   if field=='rotations':b['settings']['n_rotations']=32
   if field=='input':b['inputs'][0]['mtime_ns']=2
   if field=='code':b['code']['engine']='different'
   self.assertFalse(equivalent_grid(a,b))
 def fixture(self):
  t=tempfile.TemporaryDirectory();self.addCleanup(t.cleanup);r=Path(t.name)
  for name in ['source','left','right','code','output']:(r/name).write_text(name)
  return SimpleNamespace(cache=str(r/'cache.json'),source=str(r/'source'),left_surf=str(r/'left'),right_surf=str(r/'right'),scan_idx='',regress_enable=1,distance_mm=20.,sigma_mm=1.7,regression_code=str(r/'code'),pipeline_code=str(r/'code'),prepared_output=str(r/'output'))
 def test_prepared_data_requires_unchanged_context_and_output(self):
  a=self.fixture();self.assertEqual(check(a),'');save(a);self.assertEqual(check(a),str(Path(a.prepared_output).resolve()));a.sigma_mm=2.;self.assertEqual(check(a),'');a.sigma_mm=1.7;Path(a.prepared_output).write_text('changed');self.assertEqual(check(a),'')
 def test_input_change_during_preparation_rejected(self):
  a=self.fixture();check(a);Path(a.source).write_text('changed')
  with self.assertRaises(ValueError):save(a)
if __name__=='__main__':unittest.main()
