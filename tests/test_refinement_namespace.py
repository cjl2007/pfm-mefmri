from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lib'))
from mefmri_refinement_namespace import reserve
from mefmri_refinement_store import publish, digest


class NamespaceTests(unittest.TestCase):
    def test_frozen_reserves_new_baseline_and_rejects_calibration_overlap(self):
        with tempfile.TemporaryDirectory() as temp:
            subject = Path(temp) / 'subject'
            first, root = reserve(subject, 'rest', 'rebuild')
            original = subject / 'func/xfms' / first / 'average'
            original.write_text('calibration target')
            before = digest(original)
            version = publish(Path(root)/'references', {'template': original},
                              [{'id': 'S1_R1', 'session': 1}, {'id': 'S2_R1', 'session': 2}],
                              {'subject': subject.name}, {})
            with self.assertRaises(ValueError):
                reserve(subject, 'rest', 'frozen', version, 2)
            second, same_root = reserve(subject, 'rest', 'frozen', version, 3)
            self.assertNotEqual(first, second)
            self.assertEqual(root, same_root)
            self.assertEqual(before, digest(original))
            self.assertEqual(list((subject/'func/xfms'/second).iterdir()), [])


if __name__ == '__main__':
    unittest.main()
