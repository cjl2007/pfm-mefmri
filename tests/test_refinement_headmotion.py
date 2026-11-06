"""Execute the actual non-MEDIC headmotion module, with real FSL commands."""
import json
import os
from pathlib import Path
import subprocess
import shutil
import sys
import tempfile
import unittest
import nibabel as nib
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'lib'))
from mefmri_subject_refinement import inspect_manifest, execute
from mefmri_refinement_store import publish, digest
from mefmri_refinement_consumer import selection, check_completed
WARPKIT = os.environ.get('MEFMRI_TEST_WARPKIT_BIN') or shutil.which('wk-apply-warp') or str(Path.home()/'miniconda3/envs/warpkit_env/bin/wk-apply-warp')


class HeadmotionTests(unittest.TestCase):
    def exercise(self, mode, space, fail_final=False):
        with tempfile.TemporaryDirectory(prefix='headmotion-chain-') as temp:
            p = Path(temp)
            subject = p / 'TEST'
            proc = subject / 'func/rest/session_1/run_1'
            raw = subject / 'func/unprocessed/rest/session_1/run_1'
            xfms = subject / 'func/xfms/rest'
            for folder in (proc, raw, xfms):
                folder.mkdir(parents=True)
            env = dict(os.environ, PATH=str(Path(sys.executable).parent)+':/usr/local/fsl/bin:'+os.environ['PATH'],
                       FSLDIR='/usr/local/fsl', FSLOUTPUTTYPE='NIFTI_GZ', OMP_NUM_THREADS='1')
            def cmd(args):
                return subprocess.run([str(a) for a in args], env=env, check=True, capture_output=True, text=True)
            shape = (24, 26, 28)
            ijk = np.indices(shape, dtype=float)
            x, y, z = ijk
            brain = (100*np.exp(-((x-11)**2+(y-12)**2+(z-13)**2)/45)
                     + 40*np.exp(-((x-7)**2+(y-16)**2+(z-9)**2)/12)).astype(np.float32)
            affine = np.diag([-2., 2., 2., 1.])
            ref = p / 'grid.nii.gz'
            native = proc / 'SBref.nii.gz'
            for dest in (ref, native):
                nib.save(nib.Nifti1Image(brain, affine), dest)
            mask_path = xfms / ('T1w_nonlin_brain_func_mask.nii.gz' if space == 'MNINonlinear' else 'T1w_acpc_brain_func_mask.nii.gz')
            nib.save(nib.Nifti1Image(np.ones(shape, np.float32), affine), mask_path)
            n = 4
            for echo, scale in ((1, 1.), (2, .6)):
                image = nib.Nifti1Image(np.repeat((brain*scale)[..., None], n, axis=3), affine)
                image.header.set_zooms((2, 2, 2, 1.6))
                nib.save(image, raw / f'Rest_S1_R1_E{echo}.nii.gz')
            (proc / 'TE.txt').write_text('15 30\n')
            (proc / 'TR.txt').write_text('1.6\n')
            np.savetxt(proc / 'SliceTiming.txt', np.zeros(shape[2]))
            base = np.eye(4)
            base[:3, 3] = [.7, -.3, .4]
            correction = np.eye(4)
            angle = np.deg2rad(2)
            correction[:2, :2] = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
            correction[:3, 3] = [-.5, 1.1, .2]
            np.savetxt(p/'base.mat', base)
            np.savetxt(p/'C.mat', correction)
            warp = p / 'base_warp.nii.gz'
            cmd(['convertwarp', '--ref='+str(ref), '--premat='+str(p/'base.mat'), '--out='+str(warp), '--relout'])
            field = np.zeros((*shape, 3), np.float32)
            if space == 'MNINonlinear':
                field[..., 0] = 1.5 * np.sin(y / 6)
                field[..., 1] = .8 * np.cos(z / 7)
            header = nib.load(warp).header.copy()
            standard = p / 'standard.nii.gz'
            nib.save(nib.Nifti1Image(field, affine, header), standard)
            run = dict(id='S1_R1', session=1, native=str(native), target=str(native), warp=str(warp), affine=str(p/'base.mat'))
            data = dict(subject='TEST', func_dir='rest', distortion_mode=mode,
                        output_space=space, acpc_grid=str(ref), output_grid=str(ref),
                        standard_warp=str(standard) if space == 'MNINonlinear' else None,
                        compatibility={}, runs=[run])
            members, _, compatibility = inspect_manifest(data)
            store = p / 'store'
            version = publish(store, {'S1_R1.mat': p/'C.mat'}, members, compatibility,
                              {'scope': 'known correction for consumer integration test'})
            bundle = p / 'bundle'
            execute(data, store, bundle, 'frozen', version)
            for name, target in [('IntermediateCoregTarget.txt', native),
                                 ('Intermediate2ACPCWarp.txt', warp),
                                 ('Intermediate2ACPCMat.txt', p/'base.mat')]:
                (proc / name).write_text(str(target)+'\n')
            if mode == 'medic':
                dmap = p / 'fields.nii.gz'
                nib.save(nib.Nifti1Image(np.zeros(shape+(n,), np.float32), affine), dmap)
                (proc/'MEDICDisplacementMap.txt').write_text(str(dmap)+'\n')
                (proc/'MEDICUnwarpedSBRef.txt').write_text(str(native)+'\n')
                (proc/'MEDICPhaseEncodingAxis.txt').write_text('j-\n')
            # Forwarding wrappers log actual final sampling calls; scientific
            # executables are still FSL, not mocked stage implementations.
            shim = p / 'bin'
            shim.mkdir()
            trace = p / 'sampling.jsonl'
            for name in ('applywarp', 'flirt'):
                path = shim / name
                path.write_text('#!'+sys.executable+'\nimport json,os,sys\n'
                                'with open(os.environ["SAMPLING_TRACE"],"a") as f: f.write(json.dumps(sys.argv)+"\\n")\n'
                                'if os.environ.get("FAIL_FINAL_SAMPLING") == "1": sys.exit(17)\n'
                                'os.execv("/usr/local/fsl/bin/'+name+'",sys.argv)\n')
                path.chmod(0o755)
            if mode == 'medic':
                wk = shim/'wk-apply-warp'
                wk.write_text('#!'+sys.executable+'\nimport json,os,sys\n'
                              'with open(os.environ["SAMPLING_TRACE"],"a") as f: f.write(json.dumps(sys.argv)+"\\n")\n'
                              'if os.environ.get("FAIL_FINAL_SAMPLING") == "1": sys.exit(17)\n'
                              'os.execv('+repr(WARPKIT)+',sys.argv)\n')
                wk.chmod(0o755)
                env.update(WARPKIT_APPLY_WARP_BIN=str(wk), WARPKIT_ACTIVATE_MODE='direct', MEDIC_RESAMPLE_PYTHON=sys.executable)
            env.update(PATH=str(shim)+':'+env['PATH'], SAMPLING_TRACE=str(trace), PIPELINE_PYTHON=sys.executable,
                       FAIL_FINAL_SAMPLING='1' if fail_final else '0',
                       COREG_ALIGNMENT_STRATEGY='subject_refined', COREG_ALIGNMENT_BUNDLE=str(bundle),
                       COREG_REFERENCE_VERSION=version, FUNC_XFMS_DIRNAME='rest',
                       FUNC_NOFIELDMAP_MODE='1' if mode == 'none' else '0',
                       DISTORTION_CORRECTION_MODE=mode, APPLY_N4_BIAS='0',
                       HEADMOTION_KEEP_MCF='1', FUNC_REORIENT_TO_STD='0',
                       SBREF_ECHO_COMBINATION='mean', SBREF_MAX_TE_MS='60')
            module = 'mefmri_func_headmotion_medic.sh' if mode == 'medic' else 'mefmri_func_headmotion.sh'
            completed = subprocess.run(['bash', str(ROOT/'modules'/module),
                        str(ROOT), 'TEST', str(p), str(ref), '6', '1', '1', space, 'rest', 'Rest'],
                        env=env, capture_output=True, text=True, timeout=180)
            if fail_final:
                self.assertNotEqual(completed.returncode, 0, 'A failed final warp must stop the module')
                self.assertFalse((proc/'BOLDAlignment.json').exists())
                self.assertEqual(len(trace.read_text().splitlines()), 1, 'Stop at the first failed final warp')
                return
            self.assertEqual(completed.returncode, 0, completed.stdout+completed.stderr)
            selected = selection(proc, 'subject_refined', bundle, version)
            check_completed(proc, selected)
            calls = [json.loads(line) for line in trace.read_text().splitlines()]
            if mode == 'medic':
                self.assertEqual(len(calls), 2, 'Warpkit resamples only two estimation scouts')
                self.assertTrue(all(Path(call[0]).name == 'wk-apply-warp' for call in calls))
                sampling = json.loads((proc/'Rest_MEDIC_single_resampling.json').read_text())
                self.assertEqual(sampling['method'], 'joint_native_slice_temporal_spatial_v2')
                self.assertEqual(sampling['input_sha256'][selected['warp']], digest(selected['warp']))
                for echo in (1, 2):
                    source = raw/f'Rest_S1_R1_E{echo}.nii.gz'
                    self.assertEqual(sampling['input_sha256'][str(source)], digest(source))
                    actual = nib.load(proc/f'Rest_E{echo}_acpc.nii.gz').get_fdata()
                    self.assertEqual(actual.shape, shape+(n,))
                    self.assertTrue(np.isfinite(actual).all())
                    self.assertGreater(np.linalg.norm(actual), 0)
                return
            self.assertEqual(len(calls), 2*n, 'Exactly one final spatial interpolation per echo/frame')
            expected_tool = 'flirt' if mode == 'none' and space == 'T1w' else 'applywarp'
            self.assertTrue(all(Path(call[0]).name == expected_tool for call in calls))
            # Independent pull-field calculation: inverse(B) inverse(C) (q+D).
            points = np.moveaxis(2*ijk, 0, -1).reshape(-1, 3)
            mapped = (np.linalg.inv(base) @ np.linalg.inv(correction) @
                      np.c_[points+field.reshape(-1, 3), np.ones(len(points))].T)[:3].T
            expected_warp = p / 'analytic_warp.nii.gz'
            nib.save(nib.Nifti1Image((mapped-points).reshape((*shape, 3)).astype(np.float32), affine, header), expected_warp)
            env['PATH'] = str(Path(sys.executable).parent)+':/usr/local/fsl/bin:'+os.environ['PATH']
            for echo, scale in ((1, 1.), (2, .6)):
                frame = p / f'frame{echo}.nii.gz'
                nib.save(nib.Nifti1Image(brain*scale, affine), frame)
                actual = nib.load(proc/f'Rest_E{echo}_acpc.nii.gz').get_fdata()
                self.assertEqual(actual.shape, shape+(n,))
                for t in range(n):
                    dest = p / 'expected.nii.gz'
                    cmd(['applywarp', '--in='+str(frame), '--ref='+str(ref), '--warp='+str(expected_warp),
                         '--premat='+str(proc/'MCF'/f'MAT_{t:04d}'), '--out='+str(dest), '--interp=spline'])
                    expected = nib.load(dest).get_fdata()
                    interior = np.s_[5:-5, 5:-5, 5:-5]
                    error = np.linalg.norm((actual[..., t]-expected)[interior])/np.linalg.norm(expected[interior])
                    self.assertLess(error, 1e-4, (mode, space, echo, t, error))

    def test_topup(self):
        self.exercise('topup', 'T1w')

    def test_direct_b0(self):
        self.exercise('direct_b0', 'T1w')

    def test_no_fieldmap_anatomical(self):
        self.exercise('none', 'T1w')

    def test_no_fieldmap_nonlinear(self):
        self.exercise('none', 'MNINonlinear')

    def test_final_warp_failure_cannot_be_marked_complete(self):
        self.exercise('topup', 'T1w', fail_final=True)

    @unittest.skipUnless(Path(WARPKIT).is_file(), 'warpkit required')
    def test_medic_real_module_native_sources_and_selected_warp(self):
        self.exercise('medic', 'T1w')

    @unittest.skipUnless(Path(WARPKIT).is_file(), 'warpkit required')
    def test_medic_failed_scout_stops_before_final_sampling(self):
        self.exercise('medic', 'T1w', fail_final=True)


if __name__ == '__main__':
    unittest.main()
