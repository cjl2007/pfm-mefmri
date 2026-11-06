import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lib'))
from mefmri_refinement_consumer import selection, mark_completed, check_completed


class ConsumerTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.run = self.root / 'session_2/run_4'
        self.run.mkdir(parents=True)
        for pointer, name in [('IntermediateCoregTarget.txt', 'target'), ('Intermediate2ACPCWarp.txt', 'warp')]:
            path = self.root / name
            path.write_text('baseline ' + name)
            (self.run / pointer).write_text(str(path) + '\n')

    def test_completed_record_rejects_changed_alignment_and_outputs(self):
        result = selection(self.run, 'independent')
        output = self.run / 'Rest_E1_acpc.nii.gz'
        output.write_bytes(b'test output')
        mark_completed(self.run, result, 'Rest')
        check_completed(self.run, result)
        changed = dict(result, binding_sha256='different')
        with self.assertRaises(ValueError):
            check_completed(self.run, changed)
        output.write_bytes(b'changed output')
        with self.assertRaises(ValueError):
            check_completed(self.run, result)

    def test_missing_refinement_does_not_fall_back(self):
        with self.assertRaises(ValueError):
            selection(self.run, 'subject_refined')
        with self.assertRaises(ValueError):
            selection(self.run, 'typo')
        with self.assertRaises(ValueError):
            check_completed(self.run, dict(strategy='subject_refined'))

    def test_missing_echo_or_wrong_frame_count_cannot_be_completed(self):
        import nibabel as nib
        import numpy as np
        result = selection(self.run, 'independent')
        first = self.run / 'Rest_E1_acpc.nii.gz'
        second = self.run / 'Rest_E2_acpc.nii.gz'
        image = nib.Nifti1Image(np.ones((4, 5, 6, 7), np.float32), np.eye(4))
        nib.save(image, first)
        with self.assertRaisesRegex(ValueError, 'echo set'):
            mark_completed(self.run, result, 'Rest', 2, 7)
        self.assertFalse((self.run / 'BOLDAlignment.json').exists())
        nib.save(image, second)
        with self.assertRaisesRegex(ValueError, 'frame count'):
            mark_completed(self.run, result, 'Rest', 2, 8)
        nib.save(nib.Nifti1Image(np.ones((4, 5, 6, 6), np.float32), np.eye(4)), second)
        with self.assertRaisesRegex(ValueError, 'different dimensions'):
            mark_completed(self.run, result, 'Rest', 2, 7)
        nib.save(image, second)
        mark_completed(self.run, result, 'Rest', 2, 7)
        check_completed(self.run, result)

    def test_shell_rejects_unknown_or_missing_configuration(self):
        lib = Path(__file__).resolve().parents[1] / 'lib'
        script = 'source "$1"; refinement_select "$2" "$3"'
        for strategy in ('typo', 'subject_refined'):
            env = dict(os.environ, COREG_ALIGNMENT_STRATEGY=strategy, PIPELINE_PYTHON=sys.executable)
            env.pop('COREG_ALIGNMENT_BUNDLE', None)
            env.pop('COREG_REFERENCE_VERSION', None)
            completed = subprocess.run(['bash', '-c', script, 'test', str(lib / 'mefmri_refinement_consumer.sh'),
                                        str(lib.parent), str(self.run)], env=env, capture_output=True)
            self.assertNotEqual(completed.returncode, 0)
            self.assertEqual((self.root / 'warp').read_text(), 'baseline warp')


if __name__ == '__main__':
    unittest.main()
