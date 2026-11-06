"""QC must use selected reference chains, including the saved MEDIC prematrix."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import nibabel as nib
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'lib'))
from mefmri_subject_refinement import inspect_manifest, execute
from mefmri_refinement_store import publish
import mefmri_qc_coreg_sdc_report as qc


class QCTests(unittest.TestCase):
    def test_refined_report_without_baseline_qa_and_saved_medic_premat(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            run = p / 'func/rest/session_1/run_1'
            xfms = p / 'func/xfms/rest'
            run.mkdir(parents=True)
            xfms.mkdir(parents=True)
            xyz = np.indices((18, 20, 22), dtype=float)
            array = np.exp(-sum((xyz[k]-v)**2 for k, v in enumerate((8, 9, 10)))/30).astype(np.float32)
            affine = np.diag([-2., 2., 2., 1.])
            native = run / 'SBref.nii.gz'
            grid = xfms / 'T1w_acpc_func.nii.gz'
            for path in (native, grid):
                nib.save(nib.Nifti1Image(array, affine), path)
            identity = p / 'identity.mat'
            np.savetxt(identity, np.eye(4))
            correction = p / 'correction.mat'
            c = np.eye(4)
            c[0, 3] = 1.2
            np.savetxt(correction, c)
            warp = p / 'base.nii.gz'
            env = dict(os.environ, PATH='/usr/local/fsl/bin:'+os.environ['PATH'], FSLOUTPUTTYPE='NIFTI_GZ')
            subprocess.run(['convertwarp', '--ref='+str(grid), '--premat='+str(identity), '--out='+str(warp), '--relout'], env=env, check=True)
            data = dict(subject=p.name, func_dir='rest', distortion_mode='medic', acpc_grid=str(grid),
                        output_grid=str(grid), output_space='T1w', compatibility={},
                        runs=[dict(id='S1_R1', session=1, native=str(native), target=str(native),
                                   warp=str(warp), premat=str(identity))])
            members, _, compatibility = inspect_manifest(data)
            version = publish(p/'store', {'S1_R1.mat': correction, 'template.nii.gz': native}, members, compatibility, {})
            bundle = p / 'bundle'
            execute(data, p/'store', bundle, 'frozen', version)
            env.update(COREG_ALIGNMENT_STRATEGY='subject_refined', COREG_ALIGNMENT_BUNDLE=str(bundle), COREG_REFERENCE_VERSION=version)
            with patch.dict(os.environ, env):
                runs = qc.discover_runs(p, 'rest', 'rest', 'topup')
                self.assertEqual(len(runs), 1, 'Refined QC must not require obsolete baseline QA images')
                self.assertEqual(runs[0].corrected_acpc, bundle/'S1_R1/refined_reference.nii.gz')
                reportlets = p / 'reportlets'
                reportlets.mkdir()
                panel = qc.write_coreg(p, 'rest', reportlets)
                self.assertTrue((reportlets / panel['coreg_panel']).is_file())
                for name, value in [('MEDICUnwarpedSBRef.txt', native), ('IntermediateCoregTarget.txt', native), ('Intermediate2ACPCWarp.txt', warp)]:
                    (run/name).write_text(str(value)+'\n')
                info = runs[0]
                info.unwarped = native
                info.method = 'medic'
                actual_run = subprocess.run
                with patch.object(qc.subprocess, 'run', wraps=actual_run) as calls:
                    pair = qc.resample_medic_sdc_pair_to_acpc(info, grid, p/'sdc')
                self.assertTrue(all(call.args[0][0] != 'flirt' for call in calls.call_args_list),
                                'QC must not re-estimate a different MEDIC reference alignment')
                self.assertEqual(len(calls.call_args_list), 2)
                expected = nib.load(bundle/'S1_R1/refined_reference.nii.gz').get_fdata()
                np.testing.assert_allclose(nib.load(pair[1]).get_fdata(), expected, atol=1e-5)


if __name__ == '__main__':
    unittest.main()
