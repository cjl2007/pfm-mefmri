"""Public templates must not depend on a developer's installation layout."""
from pathlib import Path
import importlib.util,os,re,subprocess,tempfile,unittest
from unittest.mock import patch
import yaml
ROOT=Path(__file__).resolve().parents[1]

class ReleasePathsTests(unittest.TestCase):
 def test_public_configs_have_no_literal_absolute_path_assignments(self):
  for p in (ROOT/'config').glob('*.sh'):
   for line in p.read_text().splitlines():
    if line.lstrip().startswith('#'):continue
    self.assertIsNone(re.match(r'''\s*(?:export\s+)?\w+\s*=\s*["']?/''',line),f'{p.name}: {line}')
  def walk(value):
   if isinstance(value,dict):
    for v in value.values():walk(v)
   elif isinstance(value,list):
    for v in value:walk(v)
   elif isinstance(value,str):self.assertFalse(value.startswith('/'),value)
  walk(yaml.safe_load((ROOT/'config/tans_study.example.yaml').read_text()))
 def test_hcp_setup_preserves_site_overrides_and_discovers_msm_on_path(self):
  with tempfile.TemporaryDirectory(prefix='site environment ') as folder:
   base=Path(folder);fsl=base/'fsl';wb=base/'workbench';msm=base/'msm';site=base/'site msm'
   for d in [fsl/'bin',wb,msm,site]:d.mkdir(parents=True)
   for p in [fsl/'bin/fslmaths',wb/'wb_command',msm/'msm']:
    p.write_text('#!/bin/sh\nexit 0\n');p.chmod(0o755)
   env={k:v for k,v in os.environ.items()if k not in ['MSM','MSMBINDIR','CARET7DIR','HCPPIPEDIR','MATLAB_COMPILER_RUNTIME','FSL_FIXDIR']};env.update(FSLDIR=str(fsl),CARET7DIR=str(wb),PATH=str(msm)+':'+os.environ['PATH'])
   setup=ROOT/'HCPpipelines-master/Examples/Scripts/SetUpHCPPipeline.sh'
   script='source "$1"; printf "%s\\n" "$MSMBINDIR" "$CARET7DIR" "$MATLAB_COMPILER_RUNTIME" "$FSL_FIXDIR"'
   result=subprocess.run(['bash','-eu','-c',script,'test',str(setup)],env=env,capture_output=True,text=True);self.assertEqual(result.returncode,0,result.stderr);self.assertEqual(result.stdout.splitlines(),[str(msm),str(wb),'',''])
   env.update(MSMBINDIR=str(site),MATLAB_COMPILER_RUNTIME='site-mcr',FSL_FIXDIR='site-fix')
   result=subprocess.run(['bash','-eu','-c',script,'test',str(setup)],env=env,capture_output=True,text=True);self.assertEqual(result.returncode,0,result.stderr);self.assertEqual(result.stdout.splitlines(),[str(site),str(wb),'site-mcr','site-fix'])
 def test_refinement_requires_path_or_an_explicit_fsldir(self):
  import mefmri_refinement_fit as fit
  with patch.object(fit.shutil,'which',return_value=None),patch.dict(os.environ,{},clear=True):
   with self.assertRaisesRegex(FileNotFoundError,'PATH or under FSLDIR'):fit.command(['flirt'],'unused.log')
 def test_nordic_default_is_bundled_and_generic_io_uses_matlab_path(self):
  source=(ROOT/'lib/run_nordic_mefmri.m').read_text();self.assertIn("fullfile(fileparts(mfilename('fullpath')), 'nordic')",source);self.assertNotIn('/home/',source)
  source=(ROOT/'lib/nordic/NIFTI_NORDIC.m').read_text();self.assertNotIn('/home/',source);self.assertIn("exist('load_nii', 'file')",source)

if __name__=='__main__':unittest.main()
