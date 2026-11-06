from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lib'))
from mefmri_refinement_stage import validate_membership


class MembershipTests(unittest.TestCase):
    def test_partial_or_new_unlisted_runs_fail_before_coreg(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            for s, r in ((1, 1), (2, 1), (2, 2)):
                (p / f'func/rest/session_{s}/run_{r}').mkdir(parents=True)
            members = p / 'members.txt'
            members.write_text('S2_R1\n')
            with self.assertRaisesRegex(ValueError, 'whole sessions'):
                validate_membership(p, 'rest', 2, members)
            members.write_text('S2_R1\nS2_R2\n')
            self.assertEqual(validate_membership(p, 'rest', 2, members), ['S2_R1', 'S2_R2'])
            (p / 'func/rest/session_3/run_1').mkdir(parents=True)
            with self.assertRaises(ValueError):
                validate_membership(p, 'rest', 2, members)
            members.write_text('S2_R1\nS2_R1\n')
            with self.assertRaises(ValueError):
                validate_membership(p, 'rest', 2, members)


if __name__ == '__main__':
    unittest.main()
