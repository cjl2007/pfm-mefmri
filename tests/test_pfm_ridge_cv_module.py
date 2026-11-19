"""Workflow contract: opt-in routing, failure propagation, fixed-mode preservation."""
import unittest,tempfile,os,subprocess,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class ModuleTests(unittest.TestCase):
 def invoke(self,enabled,fail=False):
  t=tempfile.TemporaryDirectory();self.addCleanup(t.cleanup);root=Path(t.name);subject=root/'Example';surfaces=subject/'anat/T1w/fsaverage_LR32k';surfaces.mkdir(parents=True)
  for h in ['L','R']:(surfaces/f'Example.{h}.midthickness.32k_fs_LR.surf.gii').touch()
  for name in ['input.dtseries.nii','distance.npy','priors.mat']:(root/name).touch()
  fake=root/'python';fake.write_text('''#!/usr/bin/python3
import sys,os,json
from pathlib import Path
if len(sys.argv)>1 and sys.argv[1]=='-':os.execv('/usr/bin/python3',['/usr/bin/python3']+sys.argv[1:])
with open(os.environ['CV_TEST_LOG'],'a') as f:f.write(json.dumps(sys.argv[1:])+'\\n')
if sys.argv[1].endswith('pfm_ridge_cv_dispatch.py'):
 if os.environ.get('CV_TEST_FAIL')=='1':sys.exit(2)
 p=Path(sys.argv[sys.argv.index('--outdir')+1]);p.mkdir(parents=True,exist_ok=True);(p/'selected_parameters.json').write_text(json.dumps(dict(lambda_value=3,spatial_weight=.1)))
''');fake.chmod(0o755);log=root/'commands.jsonl'
  env={k:v for k,v in os.environ.items() if not k.startswith(('PFM_','CONCAT_','NETWORK_'))}
  env.update(PFM_PYTHON=str(fake),PIPELINE_PYTHON=str(fake),CV_TEST_LOG=str(log),CV_TEST_FAIL=str(int(fail)),PFM_INPUT_CIFTI=str(root/'input.dtseries.nii'),PFM_DISTANCE_MATRIX=str(root/'distance.npy'),PFM_PRIORS_MAT=str(root/'priors.mat'),PFM_SUBCORT_PRIORS_NII='',PFM_OUTDIR=str(root/'output'),PFM_RF_CV_ENABLE=str(enabled),PFM_RF_CV_LAMBDAS='bad-range-ignored-in-fixed-mode',PFM_RF_LAMBDA='1',PFM_RF_SPATIAL_WEIGHT='.025',PFM_RF_SUBCORT_REGRESS_ENABLE='0',PFM_RF_SMOOTHING_KERNEL='0',PFM_AREAL_ENABLE='0',PFM_RF_CONFIDENCE_ENABLE='0',PFM_RF_REPEAT_SUPPORT_ENABLE='0',PFM_RF_REVIEW_SCENE_ENABLE='0',PFM_HOMOGENEITY_TEST_ENABLE='0',PFM_RF_VOLUME_ENABLE='0',PFM_RF_MULTINETWORK_EXPERIMENTAL_ENABLE='0')
  result=subprocess.run(['bash',str(ROOT/'modules/mefmri_func_pfm.sh'),'Example',str(root),str(ROOT),'1','rest','Rest'],env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
  commands=[json.loads(s) for s in log.read_text().splitlines()] if log.exists() else []
  return result,commands
 def test_fixed_mode_does_not_invoke_cv_or_change_pair(self):
  result,commands=self.invoke(0);self.assertEqual(result.returncode,0,result.stdout);self.assertFalse(any(c[0].endswith('pfm_ridge_cv_dispatch.py') for c in commands));mapper=next(c for c in commands if c[0].endswith('pfm_ridge_fusion.py'));self.assertEqual(mapper[mapper.index('--lambda')+1],'1');self.assertEqual(mapper[mapper.index('--spatial-weight')+1],'.025')
 def test_cv_pair_reaches_mapper_and_verification(self):
  result,commands=self.invoke(1);self.assertEqual(result.returncode,0,result.stdout);mapper=next(c for c in commands if c[0].endswith('pfm_ridge_fusion.py'));self.assertEqual(mapper[mapper.index('--lambda')+1],'3');self.assertEqual(mapper[mapper.index('--spatial-weight')+1],'0.1');self.assertTrue(any(c[0].endswith('pfm_ridge_cv_finalize.py') for c in commands));cv=next(c for c in commands if c[0].endswith('pfm_ridge_cv_dispatch.py'));rotations=Path(cv[cv.index('--rotations-cifti')+1]);self.assertEqual(rotations.name,'RidgePFM_CV_Rotations.dtseries.nii');self.assertTrue(rotations.is_file())
 def test_failed_cv_does_not_fall_back(self):
  result,commands=self.invoke(1,True);self.assertNotEqual(result.returncode,0);self.assertFalse(any(c[0].endswith('pfm_ridge_fusion.py') for c in commands))
if __name__=='__main__':unittest.main()
