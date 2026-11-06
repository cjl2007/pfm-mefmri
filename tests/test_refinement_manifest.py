import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import nibabel as nib
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lib'))
from mefmri_refinement_manifest import discover


class ManifestTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.subject = self.root / 'Subject with spaces'
        self.xfms = self.subject / 'func/xfms/rest'
        self.run = self.subject / 'func/rest/session_3/run_3'
        self.run.mkdir(parents=True)
        (self.run / 'TE.txt').write_text('15 30\n')
        (self.run / 'TR.txt').write_text('1.6\n')
        (self.run / 'PE.txt').write_text('-y\n')
        self.xfms.mkdir(parents=True)
        self.grid = self.image(self.xfms / 'T1w_acpc_brain_func.nii.gz')
        self.native = self.image(self.run / 'SBref.nii.gz')
        self.avg = self.image(self.xfms / 'AvgSBref.nii.gz')
        self.image(self.subject / 'anat/T1w/T1w_acpc_dc_restore_brain.nii.gz')
        self.warp = self.image(self.xfms / 'AvgSBref2acpc_EpiReg+BBR_warp.nii.gz')
        self.pointer('IntermediateCoregTarget.txt', self.native)
        self.pointer('Intermediate2ACPCWarp.txt', self.warp)

    def image(self, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        nib.save(nib.Nifti1Image(np.ones((12, 14, 16), np.float32), np.diag([-2., 2, 2, 1])), path)
        return path

    def pointer(self, name, target):
        (self.run / name).write_text(str(target) + '\n')

    def discover(self, mode='topup', space='T1w', **kwargs):
        return discover(self.subject, 'rest', 'rest', mode, space, self.grid,
                        self.root / ('work_' + str(len(list(self.root.glob('work_*'))))), **kwargs)

    def test_scan_specific_no_extra_prematrix(self):
        for mode in ('topup', 'direct_b0'):
            result = self.discover(mode)
            self.assertIsNone(result['runs'][0]['premat'])
            self.assertEqual(result['runs'][0]['target'], str(self.native))

    def test_same_geometry_different_acquisition_rejected(self):
        import shutil
        other = self.run.parent.parent / 'session_4/run_1'
        shutil.copytree(self.run, other)
        (other / 'IntermediateCoregTarget.txt').write_text(str(other / 'SBref.nii.gz')+'\n')
        (other / 'TE.txt').write_text('20 40\n')
        with self.assertRaisesRegex(ValueError, 'Mixed acquisition'):
            self.discover()
        (other / 'TE.txt').write_text('15 30\n')
        result = self.discover()
        self.assertEqual(len(result['runs']), 2)
        self.assertEqual(result['compatibility']['acquisition']['phase_encoding'], 'j-')

    def test_average_custom_qc_and_explicit_membership(self):
        self.pointer('IntermediateCoregTarget.txt', self.avg)
        qa = self.root / 'custom QC'
        qa.mkdir()
        matrix = qa / 'SBref2AvgSBref_S3_R3.mat'
        np.savetxt(matrix, np.eye(4))
        with patch.dict(os.environ, {'COREG_QA_DIR': str(qa)}):
            result = self.discover(selected_ids=['S3_R3'])
            self.assertEqual(result['runs'][0]['premat'], str(matrix))
            with self.assertRaises(ValueError):
                self.discover(selected_ids=['S3_R3', 'S8_R1'])

    def test_none_preserves_affine_and_standard_uses_acpc(self):
        matrix = self.xfms / 'baseline.mat'
        np.savetxt(matrix, np.eye(4))
        self.pointer('Intermediate2ACPCMat.txt', matrix)
        self.assertEqual(self.discover('none')['runs'][0]['affine'], str(matrix))
        standard = self.image(self.xfms / 'AvgSBref2nonlin_EpiReg+BBR_warp.nii.gz')
        self.pointer('Intermediate2ACPCWarp.txt', standard)
        result = self.discover(space='MNINonlinear')
        self.assertEqual(result['runs'][0]['warp'], str(self.warp))
        self.assertEqual(result['runs'][0]['selected_baseline_warp'], str(standard))

    def test_medic_estimates_reference_prematrix_only(self):
        target = self.image(self.xfms / 'AvgMEDICSBref.nii.gz')
        unwarped = self.image(self.run / 'MEDIC/reference/SBref_unwarped.nii.gz')
        self.pointer('IntermediateCoregTarget.txt', target)
        self.pointer('MEDICUnwarpedSBRef.txt', unwarped)
        def fake(args, log):
            np.savetxt(args[args.index('-omat') + 1], np.eye(4))
        with patch('mefmri_refinement_manifest.command', side_effect=fake) as call:
            result = self.discover('medic')
            args = call.call_args[0][0]
            self.assertEqual(args[args.index('-ref') + 1], target)
            self.assertEqual(args[args.index('-in') + 1], unwarped)
        self.assertEqual(result['runs'][0]['target'], str(target))
        self.assertEqual(result['runs'][0]['native'], str(unwarped))
        self.assertTrue(Path(result['runs'][0]['premat']).is_file())


if __name__ == '__main__':
    unittest.main()
