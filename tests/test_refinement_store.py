import sys
from pathlib import Path
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lib'))
from mefmri_refinement_store import publish, verify, publication_lock, require_run_binding


class StoreTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / 'source'
        self.source.write_bytes(b'reference bytes')
        self.store = self.root / 'versions'
        self.members = [{'id': 'S1_R1', 'source_sha256': 'example-source-hash'}]
        self.compat = {'subject': 'fixture', 'space': 'acpc'}

    def publish(self, members=None):
        return publish(self.store, {'template.nii': self.source}, members or self.members, self.compat, {'iterations': 2})

    def test_immutable_reuse_and_changed_membership(self):
        first = self.publish()
        before = (self.store / first / 'manifest.json').read_bytes()
        self.assertEqual(first, self.publish())
        second = self.publish(self.members + [{'id': 'S2_R1', 'source_sha256': 'other'}])
        self.assertNotEqual(first, second)
        self.assertEqual(before, (self.store / first / 'manifest.json').read_bytes())
        verify(self.store, first, self.compat)

    def test_corruption_and_incompatible_incoming_run(self):
        version = self.publish()
        with self.assertRaises(ValueError):
            verify(self.store, version, {'subject': 'different'})
        (self.store / version / 'template.nii').write_bytes(b'changed')
        with self.assertRaises(ValueError):
            verify(self.store, version)

    def test_lock_and_failed_publication(self):
        with publication_lock(self.store):
            with self.assertRaises(RuntimeError):
                self.publish()
        with self.assertRaises(FileNotFoundError):
            publish(self.store, {'template': self.root / 'missing'}, self.members, self.compat, {})
        self.assertEqual(list(self.store.iterdir()), [])
        self.publish()

    def test_stale_downstream_binding_and_invalid_paths(self):
        require_run_binding({'version': 'a'}, {'version': 'a'})
        with self.assertRaises(ValueError):
            require_run_binding({'version': 'a'}, {'version': 'b'})
        with self.assertRaises(ValueError):
            publish(self.store, {'../escape': self.source}, self.members, self.compat, {})
        with self.assertRaises(ValueError):
            verify(self.store, '../escape')


if __name__ == '__main__':
    unittest.main()
