import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lib'))
from mefmri_refinement_store import digest, publish
from mefmri_subject_refinement import verify_alignment, select_run


class BindingTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.source = self.root / 'source'
        self.source.write_text('original source')
        self.compat = {'subject': 'test'}
        self.store = self.root / 'store'
        self.version = publish(self.store, {'template': self.source}, [{'id': 'S1_R1'}], self.compat, {})
        self.output = self.root / 'output'
        run = self.output / 'S1_R1'
        run.mkdir(parents=True)
        for name in ('warp', 'matrix'):
            (run / name).write_text('saved ' + name)
        self.record = dict(schema=1, complete=True, strategy='subject_refined', version=self.version,
                           reference_store=str(self.store), compatibility=self.compat,
                           source_sha256={str(self.source): digest(self.source)},
                           files={str(p.relative_to(self.output)): digest(p) for p in run.iterdir()},
                           runs=[dict(id='S1_R1', warp='S1_R1/warp', correction='S1_R1/matrix',
                                      affine=None, motion_target=str(self.source))])
        self.write()

    def write(self):
        (self.output / 'alignment.json').write_text(json.dumps(self.record))

    def test_selection_stable_and_explicit_version(self):
        one = select_run(self.output, 'S1_R1', self.version)
        self.assertEqual(one, select_run(self.output, 'S1_R1', self.version))
        self.assertEqual(one['motion_target'], str(self.source))
        with self.assertRaises(ValueError):
            select_run(self.output, 'S2_R1')
        with self.assertRaises(ValueError):
            select_run(self.output, 'S1_R1', '0' * 64)

    def test_source_change_rejected(self):
        self.source.write_text('changed input')
        with self.assertRaises(ValueError):
            verify_alignment(self.output)

    def test_output_change_rejected(self):
        (self.output / 'S1_R1/warp').write_text('changed warp')
        with self.assertRaises(ValueError):
            verify_alignment(self.output)

    def test_unhashed_selection_and_escape_rejected(self):
        self.record['runs'][0]['warp'] = '../source'
        self.write()
        with self.assertRaises(ValueError):
            verify_alignment(self.output)
        self.record['files']['../source'] = digest(self.source)
        self.write()
        with self.assertRaises(ValueError):
            verify_alignment(self.output)


if __name__ == '__main__':
    unittest.main()
