import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
import os,tempfile,unittest,json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import nibabel as nib
from mefmri_medic_resample import preflight,run
os.environ['PATH']='/usr/local/fsl/bin:'+os.environ['PATH']
class InterfaceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='medic-interface-');self.p=Path(self.temp.name)
        a=np.diag([-2.,2.,2.,1.]);self.raw=np.random.default_rng(8).normal(100,2,(6,7,8,5)).astype(np.float32)
        im=nib.Nifti1Image(self.raw,a);im.header.set_zooms((2,2,2,1.6));nib.save(im,self.p/'echo.nii.gz')
        for name,data in [('map',np.zeros(self.raw.shape,np.float32)),('ref',np.ones(self.raw.shape[:3],np.float32)),('warp',np.zeros((*self.raw.shape[:3],3),np.float32))]:
            nib.save(nib.Nifti1Image(data,a),self.p/f'{name}.nii.gz')
        (self.p/'timing.txt').write_text('0\n'*8);(self.p/'mats').mkdir()
        for t in range(5):np.savetxt(self.p/'mats'/f'MAT_{t:04d}',np.eye(4))
        self.args=SimpleNamespace(inputs=[str(self.p/'echo.nii.gz')],dmap=str(self.p/'map.nii.gz'),axis='j-',timing=str(self.p/'timing.txt'),base_warp=str(self.p/'warp.nii.gz'),reference=str(self.p/'ref.nii.gz'),mask=str(self.p/'ref.nii.gz'),output_dir=str(self.p/'outputs'),motion_dir=str(self.p/'mats'),discard=0,tr=1.6,prefix='Rest',reorient=False,preflight_only=False)
    def tearDown(self):self.temp.cleanup()
    def test_identity_full_interface(self):
        run(self.args)
        result=nib.load(self.p/'outputs/Rest_E1_acpc.nii.gz').get_fdata()
        np.testing.assert_allclose(result,self.raw,atol=3e-5)
        manifest=json.loads((self.p/'outputs/Rest_MEDIC_single_resampling.json').read_text())
        self.assertTrue(manifest['complete']);self.assertEqual(len(manifest['frame_pairs']),5)
    def test_bad_counts_and_timing_preflight(self):
        self.args.discard=1
        with self.assertRaises(ValueError):preflight(self.args)
        self.args.discard=0;(self.p/'timing.txt').write_text('0\n'*2)
        with self.assertRaises(ValueError):preflight(self.args)
        self.assertFalse((self.p/'outputs').exists())
    def test_discard_and_static_field(self):
        self.args.discard=2
        for t in [3,4]:(self.p/'mats'/f'MAT_{t:04d}').unlink()
        im=nib.load(self.p/'map.nii.gz')
        nib.save(nib.Nifti1Image(np.zeros(self.raw.shape[:3],np.float32),im.affine),self.p/'map.nii.gz')
        run(self.args)
        result=nib.load(self.p/'outputs/Rest_E1_acpc.nii.gz').get_fdata()
        np.testing.assert_allclose(result,self.raw[...,2:],atol=3e-5)
        d=json.loads((self.p/'outputs/Rest_MEDIC_single_resampling.json').read_text())
        self.assertEqual(d['frame_pairs'],[[2,0],[3,0],[4,0]])
        self.assertEqual(len(list(Path(d['motion_archive']).glob('MAT_*'))),3)
    def test_nonidentical_echo_grid_rejected(self):
        im=nib.load(self.p/'echo.nii.gz');aff=im.affine.copy();aff[0,3]+=1
        nib.save(nib.Nifti1Image(self.raw,aff,im.header),self.p/'echo2.nii.gz')
        self.args.inputs.append(str(self.p/'echo2.nii.gz'))
        with self.assertRaises(ValueError):preflight(self.args)
    def test_concurrent_run_rejected_without_output_change(self):
        out=self.p/'outputs';out.mkdir();(out/'.Rest_medic.lock').write_text('fixture')
        with self.assertRaises(RuntimeError):run(self.args)
        self.assertEqual((out/'.Rest_medic.lock').read_text(),'fixture')
    def test_failure_preserves_existing_output(self):
        out=self.p/'outputs';out.mkdir();target=out/'Rest_E1_acpc.nii.gz';target.write_bytes(b'previous-output')
        im=nib.load(self.p/'echo.nii.gz');bad=self.raw.copy();bad[0,0,0,0]=np.nan;nib.save(nib.Nifti1Image(bad,im.affine,im.header),self.p/'echo.nii.gz')
        with self.assertRaises(ValueError):run(self.args)
        self.assertEqual(target.read_bytes(),b'previous-output')
        self.assertFalse((out/'Rest_MEDIC_single_resampling.json').exists())
if __name__=='__main__':unittest.main()
