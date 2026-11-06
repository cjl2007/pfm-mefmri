"""Refined base warp through actual production MEDIC multi-echo interface."""
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
import nibabel as nib
import numpy as np
from scipy.spatial.transform import Rotation
LIB = Path(__file__).resolve().parents[1] / 'lib'
sys.path.insert(0, str(LIB))
# Use this checkout's MEDIC implementation, including after installation.
from mefmri_medic_resample import run
from mefmri_refinement_fit import command
from mefmri_refinement_geometry import warp_commands


class MedicTests(unittest.TestCase):
    def test_refinement_reaches_all_native_echoes_with_dynamic_fields_and_motion(self):
        with tempfile.TemporaryDirectory(prefix='refined-medic-') as temp:
            p = Path(temp)
            shape = (20, 22, 20)
            n = 12
            affine = np.diag([-2., 2., 2., 1.])
            xyz = 2 * np.indices(shape, dtype=float)
            ref = p / 'reference.nii.gz'
            nib.save(nib.Nifti1Image(np.ones(shape, np.float32), affine), ref)
            identity = p / 'identity.mat'
            np.savetxt(identity, np.eye(4))
            command(['convertwarp', '--ref='+str(ref), '--premat='+str(identity),
                     '--out='+str(p/'base.nii.gz'), '--relout'], p/'commands.jsonl')
            correction = np.eye(4)
            correction[:3, :3] = Rotation.from_euler('xyz', [1, -2, 3], degrees=True).as_matrix()
            correction[:3, 3] = [1, -.7, .4]
            np.savetxt(p/'correction.mat', correction)
            for args in warp_commands(p/'base.nii.gz', p/'correction.mat', ref, p/'refined.nii.gz'):
                command(args, p/'commands.jsonl')
            fields = np.zeros((*shape, n), np.float32)
            raw = np.zeros_like(fields)
            (p/'motion').mkdir()
            for k in range(n):
                d = 1.1 * np.sin(k * 1.7)
                m = .8 * np.cos(k * 1.3)
                fields[..., k] = d
                motion = np.eye(4)
                motion[1, 3] = m
                np.savetxt(p/'motion'/f'MAT_{k:04d}', motion)
                raw[..., k] = 100 + .3 * xyz[0] + .7 * (xyz[1] - d + m) - .2 * xyz[2]
            inputs = []
            for echo, scale in enumerate((1., .6), 1):
                image = nib.Nifti1Image(raw * scale, affine)
                image.header.set_zooms((2, 2, 2, 1.6))
                path = p/f'echo{echo}.nii.gz'
                nib.save(image, path)
                inputs.append(str(path))
            nib.save(nib.Nifti1Image(fields, affine), p/'fields.nii.gz')
            np.savetxt(p/'timing.txt', np.linspace(-.4, .4, shape[2]))
            args = SimpleNamespace(inputs=inputs, dmap=str(p/'fields.nii.gz'), axis='j-',
                timing=str(p/'timing.txt'), base_warp=str(p/'refined.nii.gz'), reference=str(ref),
                mask=str(ref), output_dir=str(p/'outputs'), motion_dir=str(p/'motion'),
                discard=0, tr=1.6, prefix='Rest', reorient=False, preflight_only=False)
            run(args)
            points = np.vstack([xyz.reshape(3, -1), np.ones(np.prod(shape))])
            q = (np.linalg.inv(correction) @ points)[:3].reshape((3, *shape))
            truth = 100 + .3*q[0] + .7*q[1] - .2*q[2]
            mask = np.zeros(shape, bool)
            mask[6:-6, 6:-6, 6:-6] = True
            for echo, scale in enumerate((1., .6), 1):
                actual = nib.load(p/f'outputs/Rest_E{echo}_acpc.nii.gz').get_fdata()
                error = float(np.max(np.abs(actual[mask] - scale*truth[mask, None])))
                print('Refined MEDIC echo', echo, 'max native-ramp error', error)
                self.assertLess(error, .01)
                unrefined = 100 + .3*xyz[0] + .7*xyz[1] - .2*xyz[2]
                self.assertGreater(float(np.mean(np.abs(actual[mask] - scale*unrefined[mask, None]))), .1)


if __name__ == '__main__':
    unittest.main()
