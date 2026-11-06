from pathlib import Path
import sys
import subprocess
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lib'))
from mefmri_refinement_downstream import check, stamp


class DownstreamTests(unittest.TestCase):
    def test_quiet_missing_preserves_failure_status_and_stale_diagnostics(self):
        with tempfile.TemporaryDirectory() as temp:
            run = Path(temp) / 'func/rest/session_1/run_1'
            run.mkdir(parents=True)
            bold = run / 'BOLDAlignment.json'
            bold.write_text('{"generation":1}')
            command = [sys.executable, str(Path(__file__).resolve().parents[1] /
                       'lib/mefmri_refinement_downstream.py'), 'check',
                       '--subject-dir', temp, '--func-dir', 'rest', '--stage', 'denoise']
            normal = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(normal.returncode, 1)
            self.assertIn('Missing denoise', normal.stderr)
            quiet = subprocess.run(command + ['--quiet-missing'], capture_output=True, text=True)
            self.assertEqual(quiet.returncode, 1)
            self.assertEqual(quiet.stderr, '')
            stamp(run, 'denoise')
            bold.write_text('{"generation":2}')
            stale = subprocess.run(command + ['--quiet-missing'], capture_output=True, text=True)
            self.assertEqual(stale.returncode, 1)
            self.assertIn('Stale denoise', stale.stderr)

    def test_new_bold_invalidates_all_prior_downstream_stages(self):
        with tempfile.TemporaryDirectory() as temp:
            run = Path(temp)
            bold = run / 'BOLDAlignment.json'
            bold.write_text('{"generation":1}')
            for stage in ('denoise', 'mgtr', 'surface', 'concat'):
                with self.assertRaises(ValueError):
                    check(run, stage)
                stamp(run, stage)
                check(run, stage)
            bold.write_text('{"generation":2}')
            for stage in ('denoise', 'mgtr', 'surface', 'concat'):
                with self.assertRaisesRegex(ValueError, 'Stale'):
                    check(run, stage)
            stamp(run, 'denoise')
            check(run, 'denoise')
            with self.assertRaises(ValueError):
                check(run, 'surface')


if __name__ == '__main__':
    unittest.main()
