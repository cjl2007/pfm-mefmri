"""Run module routing with stub scientific programs; never process subject data."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class RidgePFMModuleTests(unittest.TestCase):
    def run_module(self, **overrides):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        subject = root/'study space/SUB001'
        surfaces = subject/'anat/T1w/fsaverage_LR32k'
        surfaces.mkdir(parents=True)
        for hemi in ('L', 'R'):
            (surfaces/f'SUB001.{hemi}.midthickness.32k_fs_LR.surf.gii').touch()
        for name in ('input.dtseries.nii', 'distance.npy', 'priors.mat'):
            (root/name).touch()
        capture = root/'calls.jsonl'
        stub = root/'stub-python'
        stub.write_text(f'#!{sys.executable}\nimport json, os, sys\n'
                        'with open(os.environ["CAPTURE"], "a") as f: f.write(json.dumps(sys.argv[1:])+"\\n")\n')
        stub.chmod(0o755)
        wb = root/'wb_command'
        wb.write_text('#!/bin/sh\nexit 0\n')
        wb.chmod(0o755)
        env = {k:v for k,v in os.environ.items() if not k.startswith(('PFM_', 'CONCAT_', 'NETWORK_'))}
        env.update(PATH=str(root)+os.pathsep+env.get('PATH',''), CAPTURE=str(capture),
                   PFM_PYTHON=str(stub), PFM_INPUT_CIFTI=str(root/'input.dtseries.nii'),
                   PFM_DISTANCE_MATRIX=str(root/'distance.npy'), PFM_PRIORS_MAT=str(root/'priors.mat'),
                   PFM_RF_SUBCORT_REGRESS_ENABLE='0', PFM_RF_SMOOTHING_KERNEL='0',
                   PFM_RF_REVIEW_SCENE_ENABLE='0', PFM_AREAL_ENABLE='0', PFM_HOMOGENEITY_TEST_ENABLE='0')
        env.update(overrides)
        result = subprocess.run(['bash', str(ROOT/'modules/mefmri_func_pfm.sh'),
                                 'SUB001', str(subject.parent), str(ROOT), '1', 'rest', 'Rest'],
                                env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        calls = [json.loads(line) for line in capture.read_text().splitlines()]
        return calls, subject

    def test_default_and_legacy_strategies_use_new_output_and_weights(self):
        for strategy in (None, 'ridgepfm', 'RidgePFM', 'ridge_fusion', 'RidgeFusion', 'ridge-fusion'):
            with self.subTest(strategy=strategy):
                calls, subject = self.run_module(**({'PFM_STRATEGY':strategy} if strategy else {}))
                self.assertEqual(len(calls), 1)
                args = calls[0]
                self.assertEqual(Path(args[0]).name, 'pfm_ridge_fusion.py')
                for flag, expected in [('--outdir', str(subject/'func/rest/PFM/RidgePFM')),
                                       ('--outfile', 'RidgePFM_VTX'), ('--fc-weight', '1.0'),
                                       ('--spatial-weight', '0.05')]:
                    self.assertEqual(args[args.index(flag)+1], expected)

    def test_explicit_study_weights_and_prefix_are_preserved(self):
        calls, _ = self.run_module(PFM_RF_SPATIAL_WEIGHT='0.2', PFM_RF_FC_WEIGHT='0.8',
                                  PFM_RF_OUTFILE='StudyMap')
        for flag, expected in [('--spatial-weight','0.2'), ('--fc-weight','0.8'), ('--outfile','StudyMap')]:
            self.assertEqual(calls[0][calls[0].index(flag)+1], expected)


if __name__ == '__main__':
    unittest.main()
