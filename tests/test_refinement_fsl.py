"""Actual FSL field composition and native-image sampling against analytic pulls."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import nibabel as nib
import numpy as np
from scipy.spatial.transform import Rotation
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lib'))
from mefmri_refinement_geometry import warp_commands


class FSLTests(unittest.TestCase):
    def test_native_image_motion_refinement_and_nonlinear_target(self):
        env = dict(os.environ, FSLOUTPUTTYPE='NIFTI_GZ')
        def run(args):
            subprocess.run(['/usr/local/fsl/bin/' + args[0]] + args[1:], env=env,
                           check=True, capture_output=True, text=True)
        with tempfile.TemporaryDirectory(prefix='refinement-fsl-') as folder:
            p = Path(folder)
            shape = (32, 34, 36)
            affine = np.diag([-2., 2., 2., 1.])
            affine[:3, 3] = [30, -20, -20]
            ijk = np.indices(shape, dtype=float)
            xyz = 2 * ijk
            ref = p / 'ref.nii.gz'
            nib.save(nib.Nifti1Image(np.ones(shape, np.float32), affine), ref)
            ramp = 7 + .3 * xyz[0] + .7 * xyz[1] - .2 * xyz[2]
            native = p / 'native.nii.gz'
            nib.save(nib.Nifti1Image(ramp.astype(np.float32), affine), native)
            matrices = []
            for name, angles, shift in [('base', [1, -2, 1], [1, -.5, 1]),
                                        ('correction', [-2, 1, 3], [-1, 1, .5]),
                                        ('motion', [1, 2, -1], [.4, -.3, .6])]:
                m = np.eye(4)
                m[:3, :3] = Rotation.from_euler('xyz', angles, degrees=True).as_matrix()
                m[:3, 3] = shift
                np.savetxt(p / (name + '.mat'), m)
                matrices.append(m)
            base, correction, motion = matrices
            run(['convertwarp', '--ref=' + str(ref), '--premat=' + str(p / 'base.mat'),
                 '--out=' + str(p / 'base.nii.gz'), '--relout'])
            field = np.zeros((*shape, 3), np.float32)
            field[..., 0] = .4 * np.sin(ijk[1] / 6)
            field[..., 1] = .3 * np.cos(ijk[2] / 7)
            field[..., 2] = .2 * np.sin(ijk[0] / 5)
            # Give the test field FSL's relative-displacement intent/header.
            header = nib.load(p / 'base.nii.gz').header.copy()
            nib.save(nib.Nifti1Image(field, affine, header), p / 'standard.nii.gz')
            commands = warp_commands(p / 'base.nii.gz', p / 'correction.mat', ref,
                                     p / 'refined.nii.gz', p / 'standard.nii.gz', ref, p / 'final.nii.gz')
            for args in commands:
                run(args)
            run(['applywarp', '--in=' + str(native), '--ref=' + str(ref),
                 '--premat=' + str(p / 'motion.mat'), '--warp=' + str(p / 'final.nii.gz'),
                 '--out=' + str(p / 'sampled.nii.gz'), '--interp=trilinear'])
            points = np.moveaxis(xyz, 0, -1).reshape(-1, 3)
            standard_pull = points + field.reshape(-1, 3)
            homogeneous = np.c_[standard_pull, np.ones(len(points))].T
            expected_base = (np.linalg.inv(base) @ np.linalg.inv(correction) @ homogeneous)[:3].T
            actual_base = points + nib.load(p / 'final.nii.gz').get_fdata().reshape(-1, 3)
            interior = np.zeros(shape, bool)
            interior[7:-7, 7:-7, 7:-7] = True
            mask = interior.ravel()
            field_error = np.max(np.linalg.norm(actual_base[mask] - expected_base[mask], axis=1))
            self.assertLess(field_error, .001)
            native_pull = (np.linalg.inv(motion) @ np.vstack([expected_base.T, np.ones(len(points))]))[:3].T
            expected = 7 + native_pull @ np.array([.3, .7, -.2])
            actual = nib.load(p / 'sampled.nii.gz').get_fdata().ravel()
            image_error = np.max(np.abs(actual[mask] - expected[mask]))
            self.assertLess(image_error, .001)
            wrong = (np.linalg.inv(correction) @ np.linalg.inv(base) @ homogeneous)[:3].T
            self.assertGreater(np.max(np.linalg.norm(wrong[mask] - expected_base[mask], axis=1)), .01)
            print('Actual FSL refinement: max field error mm', field_error,
                  'max native-ramp intensity error', image_error)


if __name__ == '__main__':
    unittest.main()
