import sys
from pathlib import Path
import unittest
import numpy as np
from scipy.spatial.transform import Rotation
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lib'))
from mefmri_refinement_geometry import rigid, center_corrections, refined_affine, warp_commands


def transform(degrees, shift):
    m = np.eye(4)
    m[:3, :3] = Rotation.from_euler('xyz', degrees, degrees=True).as_matrix()
    m[:3, 3] = shift
    return m


class GeometryTests(unittest.TestCase):
    def test_motion_then_baseline_then_refinement(self):
        motion = transform([2, -1, 3], [1, 3, -2])
        base = transform([-3, 2, 1], [-5, 7, 2])
        correction = transform([1, 4, -2], [3, -2, 1])
        point = np.array([20, -30, 12, 1.])
        expected = correction @ (base @ (motion @ point))
        np.testing.assert_allclose(refined_affine(base, correction) @ motion @ point, expected)
        self.assertGreater(np.linalg.norm(base @ correction @ motion @ point - expected), .1)
        np.testing.assert_array_equal(refined_affine(base, np.eye(4)), base)

    def test_session_balance_not_run_count(self):
        a = transform([0, 0, 0], [2, 0, 0])
        b = transform([0, 0, 0], [-2, 0, 0])
        centered = center_corrections([a, a, a, b], [1, 1, 1, 2])
        np.testing.assert_allclose(centered[0], a)
        np.testing.assert_allclose(centered[-1], b)
        with self.assertRaises(ValueError):
            center_corrections([a], [1])

    def test_reject_nonrigid(self):
        for bad in [np.diag([2., 1, 1, 1]), np.diag([-1., 1, 1, 1]), np.full((4, 4), np.nan)]:
            with self.assertRaises(ValueError):
                rigid(bad)

    def test_standard_composed_after_acpc_and_inputs_preserved(self):
        commands = warp_commands('base', 'correction', 'acpc', 'refined', 'standard', 'mni', 'final')
        self.assertIn('--postmat=correction', commands[0])
        self.assertIn('--warp1=refined', commands[1])
        self.assertIn('--warp2=standard', commands[1])
        self.assertFalse(any('postmat' in arg for arg in commands[1]))
        with self.assertRaises(ValueError):
            warp_commands('base', 'correction', 'acpc', 'base')
        with self.assertRaises(ValueError):
            warp_commands('base', 'correction', 'acpc', 'out', standard_warp='standard')


if __name__ == '__main__':
    unittest.main()
