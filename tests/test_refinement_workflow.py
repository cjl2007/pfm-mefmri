"""Actual FSL reference build and frozen incoming-run workflow."""
from pathlib import Path
import sys
import tempfile
import unittest
import nibabel as nib
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lib'))
from mefmri_refinement_fit import command
from mefmri_subject_refinement import execute, select_run, verify_alignment
from mefmri_refinement_store import digest


class WorkflowTests(unittest.TestCase):
    def test_joint_then_frozen_preserves_reference_and_training_outputs(self):
        with tempfile.TemporaryDirectory(prefix='refinement-workflow-') as temp:
            p = Path(temp)
            shape = (24, 26, 28)
            x, y, z = np.indices(shape, dtype=float)
            data = (80 * np.exp(-((x-11)**2+(y-12)**2+(z-14)**2)/55)
                    + 40 * np.exp(-((x-7)**2+(y-16)**2+(z-10)**2)/12)
                    + 25 * np.exp(-((x-16)**2+(y-9)**2+(z-19)**2)/15)).astype(np.float32)
            affine = np.diag([-2., 2., 2., 1.])
            grid = p / 'grid.nii.gz'
            nib.save(nib.Nifti1Image(data, affine), grid)
            matrix = p / 'identity.mat'
            np.savetxt(matrix, np.eye(4))
            warp = p / 'warp.nii.gz'
            command(['convertwarp', '--ref='+str(grid), '--premat='+str(matrix),
                     '--out='+str(warp), '--relout'], p / 'commands.jsonl')
            def run(session):
                return dict(id=f'S{session}_R1', session=session, native=str(grid),
                            target=str(grid), warp=str(warp), premat=None)
            manifest = dict(subject='test', func_dir='rest', distortion_mode='topup',
                            acpc_grid=str(grid), output_grid=str(grid), output_space='T1w',
                            compatibility={'sequence': 'synthetic'}, runs=[run(1), run(2)])
            first = execute(manifest, p / 'store', p / 'joint', 'rebuild')
            verify_alignment(p / 'joint')
            selected = select_run(p / 'joint', 'S1_R1', first['version'])
            np.testing.assert_allclose(np.loadtxt(selected['correction']), np.eye(4), atol=.01)
            original = {str(f): digest(f) for folder in (p / 'store', p / 'joint') for f in folder.rglob('*') if f.is_file()}
            manifest['runs'] = [run(3)]
            second = execute(manifest, p / 'store', p / 'incoming', 'frozen', first['version'])
            self.assertEqual(second['version'], first['version'])
            chosen = select_run(p / 'incoming', 'S3_R1', first['version'])
            np.testing.assert_allclose(np.loadtxt(chosen['correction']), np.eye(4), atol=.05)
            self.assertEqual(original, {f: digest(f) for f in original})
            manifest['runs'] = [dict(run(1), premat=str(matrix))]
            with self.assertRaisesRegex(ValueError, 'Calibration run changed'):
                execute(manifest, p / 'store', p / 'bad', 'frozen', first['version'])


if __name__ == '__main__':
    unittest.main()
