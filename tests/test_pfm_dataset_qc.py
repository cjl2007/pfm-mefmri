"""Check real NSI shell routing/duration without invoking the scientific backend."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import nibabel as nib
import numpy as np
from nibabel.cifti2.cifti2_axes import SeriesAxis, BrainModelAxis
ROOT=Path(__file__).resolve().parents[1]

class NSISelectedInputTests(unittest.TestCase):
    def test_actual_cifti_duration_and_external_input_without_scaninfo(self):
        with tempfile.TemporaryDirectory(prefix='selected NSI ') as temp:
            base=Path(temp);package=base/'backend/pfm_nsi';package.mkdir(parents=True)
            (package/'__init__.py').touch()
            capture=base/'arguments.json'
            (package/'cli.py').write_text("import os,sys,json\nfrom pathlib import Path\nPath(os.environ['CAPTURE']).write_text(json.dumps(sys.argv[1:]))\n")
            source=base/'external selected.dtseries.nii'
            nib.save(nib.Cifti2Image(np.arange(15,dtype=np.float32).reshape(5,3),
                nib.Cifti2Header.from_axes((SeriesAxis(0,1.5,5),BrainModelAxis.from_surface(np.arange(3),3,name='CortexLeft')))),source)
            env=os.environ.copy();env.update(NSI_ENABLE='1',NSI_USE_EXTERNAL_CLI='1',NSI_PYTHON=sys.executable,
                NSI_EXTERNAL_ROOT=str(package.parent),NSI_INPUT_CIFTI=str(source),NSI_OUTDIR=str(base/'pfm/QC/NSI'),
                NSI_SCAN_INFO_JSON=str(base/'does-not-exist.json'),CAPTURE=str(capture),NSI_USABILITY_MODEL='1',NSI_RELIABILITY_MODEL='0')
            result=subprocess.run(['bash',str(ROOT/'modules/mefmri_func_nsi.sh'),'SUB',str(base),str(ROOT),'1','effort','Effort'],env=env,text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            args=json.loads(capture.read_text());self.assertEqual(args[args.index('--cifti')+1],str(source))
            self.assertEqual(float(args[args.index('--nsi-t')+1]),.125)
            self.assertEqual(float(args[args.index('--query-t')+1]),.125)
            self.assertEqual(args[args.index('--outdir')+1],str(base/'pfm/QC/NSI'))
            self.assertIn('retained input duration',result.stdout)


    def test_disabled_backend_does_not_rebind_old_outputs_to_new_input(self):
        with tempfile.TemporaryDirectory(prefix='disabled dataset QC ') as temp:
            base=Path(temp);root=base/'SUB/pfm';out=root/'QC/NSI';out.mkdir(parents=True)
            (out/'old_result.txt').write_text('historical result')
            (root/'dataset.json').write_text(json.dumps({'input_fingerprint':'new dataset','tasks':['rest']}))
            source=base/'selected.dtseries.nii'
            nib.save(nib.Cifti2Image(np.ones((5,3),np.float32),nib.Cifti2Header.from_axes((SeriesAxis(0,2,5),BrainModelAxis.from_surface(np.arange(3),3,name='CortexLeft')))),source)
            env=os.environ.copy();env.update(PFM_DATASET_ROOT=str(root),PFM_INPUT_CIFTI=str(source),PFM_CONCAT_DIR=str(root/'ConcatenatedCiftis'),NSI_ENABLE='1',NSI_USE_EXTERNAL_CLI='0',NSI_EXTERNAL_ROOT='',NSI_PYTHON=sys.executable)
            result=subprocess.run([sys.executable,str(ROOT/'lib/mefmri_pfm_qc.py'),'nsi','SUB',str(base),str(ROOT),'1','rest','Rest'],env=env,text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            receipt=json.loads((root/'QC/nsi_dataset.json').read_text())
            self.assertEqual(receipt['status'],'skipped');self.assertEqual(receipt['outputs'],[])
            self.assertEqual((out/'old_result.txt').read_text(),'historical result')

if __name__=='__main__':unittest.main()
