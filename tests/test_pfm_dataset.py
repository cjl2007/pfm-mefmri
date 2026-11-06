"""Real CIFTI preparation and provenance tests for subject-level PFM datasets."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import nibabel as nib
from nibabel.cifti2.cifti2_axes import BrainModelAxis, SeriesAxis
import numpy as np
from scipy.io import loadmat
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'lib'))
from mefmri_pfm_dataset import select
from pfm_module import compute_fd_from_mcfpar


class DatasetTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='pooled PFM ')
        self.addCleanup(self.temp.cleanup)
        self.subject = Path(self.temp.name) / 'study with spaces' / 'SUB001'
        self.specs = [{'task': 'rest', 'prefix': 'Rest'}, {'task': 'effort', 'prefix': 'Effort'}]
        self.brain = BrainModelAxis.from_surface(np.arange(3), 3, name='CortexLeft')
        self.expected = []
        for task_index, spec in enumerate(self.specs):
            for session in (1, 2):
                run = self.subject / 'func' / spec['task'] / ('session_' + str(session)) / 'run_1'
                run.mkdir(parents=True)
                (run / 'TR.txt').write_text('2.0')
                motion = np.zeros((6, 6))
                motion[:, 3] = [0, .3, .3, .30000004, .30000004, .60000008]
                np.savetxt(run / 'MCF.par', motion, fmt='%.17g')
                data = np.arange(18, dtype=np.float32).reshape(6, 3) * (task_index + 1) + session * 10
                self.expected.append(data - data.mean(axis=0, keepdims=True))
                self.write_cifti(run / (spec['prefix'] + '_Tag.dtseries.nii'), data)
            concat = self.subject / 'func' / spec['task'] / 'ConcatenatedCiftis'
            concat.mkdir()
            (concat / 'PFMPrepSummary.json').write_text(json.dumps({'input_tag': 'Tag'}))

    def write_cifti(self, path, data, tr=2., brain=None):
        header = nib.Cifti2Header.from_axes((SeriesAxis(0, tr, len(data)), self.brain if brain is None else brain))
        nib.save(nib.Cifti2Image(data, header), path)

    def prepare(self, specs=None, rebuild=0):
        return subprocess.run([sys.executable, '-B', str(ROOT / 'lib/mefmri_pfm_dataset.py'), 'prepare',
            '--subject-dir', str(self.subject), '--root', str(self.subject / 'pfm'),
            '--spec', json.dumps(self.specs if specs is None else specs), '--rebuild', str(rebuild)],
            env={**{k: v for k, v in os.environ.items() if not k.startswith(('PFM_', 'CONCAT_'))},
                 'PFM_FD_THRESHOLD': '0.3'}, capture_output=True, text=True)

    def test_pool_exact_data_censoring_and_noncolliding_run_identity(self):
        result = self.prepare()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        output = self.subject / 'pfm/ConcatenatedCiftis'
        full = np.vstack(self.expected)
        actual = np.asanyarray(nib.load(output / 'PFM_Combined_Concatenated.dtseries.nii').dataobj)
        np.testing.assert_array_equal(actual, full)
        fd = np.loadtxt(output / 'FD.txt')
        keep = fd < .3
        censored = np.asanyarray(nib.load(output / 'PFM_Combined_Concatenated+FDlt0p3.dtseries.nii').dataobj)
        np.testing.assert_array_equal(censored, full[keep])
        ids = loadmat(output / 'ScanIdx.mat')['ScanIdx']
        for index, expected in enumerate([(1, 1), (2, 1), (1, 2), (2, 2)]):
            np.testing.assert_array_equal(ids[index * 6:(index + 1) * 6], np.tile(expected, (6, 1)))
        np.testing.assert_array_equal(loadmat(output / 'ScanIdx_Censored.mat')['ScanIdx'], ids[keep])
        manifest = json.loads((self.subject / 'pfm/dataset.json').read_text())
        self.assertEqual(manifest['status'], 'prepared')
        self.assertEqual([r['task'] for r in manifest['runs']], ['rest', 'rest', 'effort', 'effort'])
        self.assertEqual([r['run'] for r in manifest['runs']], [1, 1, 1, 1])
        self.assertEqual(manifest['n_retained'], int(keep.sum()))
        self.assertEqual(manifest['functional_reference_task'], 'rest')
        self.assertFalse((self.subject / 'func/rest/PFM').exists())

    def test_rest_only_excludes_tasks(self):
        result = self.prepare(self.specs[:1])
        self.assertEqual(result.returncode, 0, result.stderr)
        manifest = json.loads((self.subject / 'pfm/dataset.json').read_text())
        self.assertEqual(len(manifest['runs']), 2)
        self.assertEqual({r['task'] for r in manifest['runs']}, {'rest'})

    def test_changed_selection_requires_rebuild(self):
        self.assertEqual(self.prepare(self.specs[:1]).returncode, 0)
        result = self.prepare()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Restart from concat', result.stderr)
        self.assertEqual(self.prepare(rebuild=1).returncode, 0)

    def test_mismatched_grayordinate_identities_rejected(self):
        path = self.subject / 'func/effort/session_1/run_1/Effort_Tag.dtseries.nii'
        brain = BrainModelAxis.from_surface(np.array([0, 2, 3]), 4, name='CortexLeft')
        self.write_cifti(path, np.ones((6, 3), dtype=np.float32), brain=brain)
        result = self.prepare()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('BrainModelAxis mismatch', result.stderr)

    def test_mismatched_tr_rejected(self):
        path = self.subject / 'func/effort/session_1/run_1/Effort_Tag.dtseries.nii'
        self.write_cifti(path, np.ones((6, 3), dtype=np.float32), tr=2.5)
        (path.parent / 'TR.txt').write_text('2.5')
        result = self.prepare()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('TR mismatch', result.stderr)

    def test_requested_missing_task_rejected(self):
        with self.assertRaisesRegex(ValueError, 'exactly one'):
            select([('rest', 'Rest')], ['effort'])
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            select([('rest', 'Rest')], ['rest', 'Rest'])
        self.assertEqual(select([('rest', 'Rest'), ('effort', 'Effort')], ['Rest']), self.specs[:1])


class TansLocationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = ROOT.parent / 'pfm-tans/pipeline/run_pipeline.py'
        spec = importlib.util.spec_from_file_location('new_tans_driver', path)
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def test_new_root_precedence_legacy_and_explicit_override(self):
        with tempfile.TemporaryDirectory() as temp:
            subject = Path(temp)
            self.assertEqual(self.module.resolve_pfm_directory(subject, 'rest'), subject / 'pfm/RidgePFM')
            legacy = subject / 'func/rest/PFM/RidgePFM'
            legacy.mkdir(parents=True)
            self.assertEqual(self.module.resolve_pfm_directory(subject, 'rest'), legacy)
            new = subject / 'pfm/RidgePFM'
            new.mkdir(parents=True)
            self.assertEqual(self.module.resolve_pfm_directory(subject, 'effort'), new)
            self.assertEqual(self.module.resolve_pfm_directory(subject, 'rest', str(legacy)), legacy)

    def test_resume_rejects_changed_or_incomplete_pfm_provenance(self):
        with tempfile.TemporaryDirectory() as temp:
            subject = Path(temp)
            pfm = subject / 'pfm/RidgePFM'
            pfm.mkdir(parents=True)
            output = subject / 'tans'
            manifest = pfm.parent / 'dataset.json'
            manifest.write_text(json.dumps({'status': 'prepared'}))
            with self.assertRaisesRegex(ValueError, 'not complete'):
                self.module.check_pfm_dataset_provenance(output, pfm, True)
            manifest.write_text(json.dumps({'status': 'complete', 'tasks': ['rest']}))
            self.module.check_pfm_dataset_provenance(output, pfm, True)
            manifest.write_text(json.dumps({'status': 'complete', 'tasks': ['rest', 'effort']}))
            with self.assertRaisesRegex(ValueError, 'provenance changed'):
                self.module.check_pfm_dataset_provenance(output, pfm, True)

    def test_all_tans_yaml_templates_render_with_new_root(self):
        subject = Path('/example/BioSelect-23002')
        variables = {'subject': subject.name, 'subject_dir': str(subject), 'func_dir': 'rest',
                     'pfm_dir': str(subject / 'pfm/RidgePFM'), 'tans_root': str(subject / 'tans'),
                     'pfm_tans_root': str(ROOT.parent / 'pfm-tans')}
        paths = list((ROOT / 'config').glob('tans*.yaml'))
        self.assertEqual(len(paths), 3)
        for path in paths:
            with self.subTest(path=path.name):
                rendered = self.module.render_strings(yaml.safe_load(path.read_text()), variables)
                self.assertEqual(rendered['inputs']['categorical_networks'], str(subject / 'pfm/RidgePFM/RidgePFM_VTX.dlabel.nii'))


if __name__ == '__main__':
    unittest.main(verbosity=2)
