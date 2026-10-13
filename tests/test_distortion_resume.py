"""Exercise real pipeline routing with synthetic inputs and scientific-stage stubs."""
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(os.environ.get('MEFMRI_TEST_ROOT', Path(__file__).resolve().parents[1]))
TOPUP_FILES = ('Avg_FM_rads_acpc.nii.gz', 'Avg_FM_mag_acpc.nii.gz',
               'Avg_FM_mag_acpc_brain.nii.gz', 'acqparams.txt')


class DistortionResumeTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='distortion resume ')
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        self.subject = self.base / 'study with spaces' / 'SUB001'
        self.events = self.base / 'events.txt'
        self.setup_script = self.base / 'environment.sh'
        self.setup_script.write_text('tkregister() { return 0; }\n')
        self.license = self.base / 'fixture-license.txt'
        self.license.write_text('No real FreeSurfer command is executed.\n')
        for modality in ('T1w', 'T2w'):
            path = self.subject / f'anat/{modality}/{modality}.nii.gz'
            path.parent.mkdir(parents=True)
            path.write_bytes(b'fixture')
        self.stubs = {}
        for name in ('validate', 'fieldmaps', 'medic', 'coreg', 'coreg_medic',
                     'qc', 'headmotion', 'headmotion_medic', 'denoise'):
            path = self.base / f'{name}.sh'
            path.write_text('#!/usr/bin/env bash\n'
                            f'printf "{name}:%s\\n" "$DISTORTION_CORRECTION_MODE" >> "$EVENTS"\n')
            path.chmod(0o755)
            self.stubs[name] = path
        self.runs = [(s, r) for s in range(1, 14) for r in (1, 2)]
        self.task = 'rest'
        self.fm_rel = 'func/rest/field_maps'
        self.make_runs()

    def make_runs(self):
        for session, run in self.runs:
            proc = self.subject / f'func/{self.task}/session_{session}/run_{run}'
            raw = self.subject / f'func/unprocessed/{self.task}/session_{session}/run_{run}'
            proc.mkdir(parents=True, exist_ok=True)
            raw.mkdir(parents=True, exist_ok=True)
            (proc / 'TE.txt').write_text('0.015 0.030 0.045\n')
            for echo in (1, 2, 3):
                (raw / f'Rest_S{session}_R{run}_E{echo}.nii.gz').write_bytes(b'fixture')

    def make_topup(self, mode='measured', missing=None):
        folder = self.subject / self.fm_rel
        folder.mkdir(parents=True, exist_ok=True)
        (folder / 'FIELDMAP_MODE.txt').write_text(mode + '\n')
        for name in TOPUP_FILES:
            if name != missing:
                (folder / name).write_bytes(b'fixture')

    def make_medic(self, runs=None, empty=False):
        for session, run in self.runs if runs is None else runs:
            folder = self.subject / f'func/{self.task}/session_{session}/run_{run}/MEDIC'
            folder.mkdir(parents=True, exist_ok=True)
            (folder / f'Rest_S{session}_R{run}_displacementmaps.nii').write_bytes(
                b'' if empty else b'fixture')

    def make_raw_phases(self):
        for path in (self.subject / f'func/unprocessed/{self.task}').glob('session_*/run_*/*_E*.nii.gz'):
            path.with_name(path.name[:-7] + '_phase.nii.gz').write_bytes(b'fixture')

    def run_pipeline(self, **overrides):
        settings = dict(
            MEDIR=ROOT, EnvironmentScript=self.setup_script, FS_LICENSE_FILE=self.license,
            START_FROM_MODULE='coreg', STOP_AFTER_MODULE='coreg', START_SESSION=1,
            DISTORTION_CORRECTION_MODE='medic', FUNC_NOFIELDMAP_MODE=0,
            FUNC_DIRNAME=self.task, FUNC_FILE_PREFIX='Rest', FUNC_XFMS_DIRNAME=self.task,
            FM_OUT_DIR_REL=self.fm_rel, PROCESSING_MODE='multi_echo',
            MULTI_ECHO_DENOISE_METHOD='meica', PIPELINE_PYTHON=sys.executable,
            TANS_ENABLE=0, PFM_ENABLE=0, CONCAT_ENABLE=0, NSI_ENABLE=0,
            CRAWLING_SEED_FC_ENABLE=0, RUN_CONFIG_SNAPSHOT=1,
            COREG_SDC_QC_ENABLE=1, MOTION_GRAYPLOTS_QC_ENABLE=0,
            PIPELINE_QUIET_MODULE_OUTPUT=0,
            VALIDATE_MODULE=self.stubs['validate'], FUNC_FIELDMAPS_MODULE=self.stubs['fieldmaps'],
            FUNC_MEDIC_MODULE=self.stubs['medic'], FUNC_COREG_MODULE=self.stubs['coreg'],
            FUNC_MEDIC_COREG_MODULE=self.stubs['coreg_medic'], FUNC_QC_COREG_SDC_MODULE=self.stubs['qc'],
            FUNC_HEADMOTION_MODULE=self.stubs['headmotion'],
            FUNC_MEDIC_HEADMOTION_MODULE=self.stubs['headmotion_medic'],
            FUNC_MEICA_MODULE=self.stubs['denoise'],
        )
        settings.update(overrides)
        config = self.base / 'study config.sh'
        config.write_text('\n'.join(f'{k}={shlex.quote(str(v))}' for k, v in settings.items()) + '\n')
        env = {k: v for k, v in os.environ.items()
               if not k.startswith(('PIPELINE_', 'FUNC_', 'FM_', 'MEDIC_', 'TANS_', 'PFM_'))
               and k not in ('BASH_ENV', 'ENV')}
        env['EVENTS'] = str(self.events)
        return subprocess.run(['bash', str(ROOT / 'bin/mefmri_pipeline.sh'), str(self.subject), str(config)],
                              env=env, capture_output=True, text=True, timeout=30)

    def assert_events(self, result, expected):
        self.assertEqual(result.returncode, 0, result.stdout + '\n' + result.stderr)
        actual = self.events.read_text().splitlines() if self.events.exists() else []
        self.assertEqual(actual, expected)

    def assert_preflight_failure(self, result):
        self.assertEqual(result.returncode, 2, result.stdout + '\n' + result.stderr)
        self.assertIn('START_FROM_MODULE=fieldmaps', result.stderr)
        self.assertFalse(self.events.exists(), result.stdout)

    def test_coreg_resume_recovers_topup_without_raw_phase_or_raw_fieldmaps(self):
        self.make_topup()
        # A failed attempt's startup metadata must not override the saved products.
        metadata = self.subject / 'func/rest/QC/RunMetadata'
        metadata.mkdir(parents=True)
        (metadata / 'pipeline_run_20000101_000000.txt').write_text('distortion_correction_mode=medic\n')
        result = self.run_pipeline()
        self.assert_events(result, ['coreg:topup', 'qc:topup'])
        snapshots = sorted(metadata.glob('pipeline_run_*.txt'))
        text = snapshots[-1].read_text()
        self.assertIn('distortion_correction_mode=topup\n', text)
        self.assertIn('distortion_correction_mode_requested=medic\n', text)

    def test_coreg_then_headmotion_both_use_topup(self):
        self.make_topup()
        result = self.run_pipeline(STOP_AFTER_MODULE='headmotion')
        self.assert_events(result, ['coreg:topup', 'qc:topup', 'headmotion:topup'])

    def test_headmotion_only_resume_recovers_topup(self):
        self.make_topup()
        result = self.run_pipeline(START_FROM_MODULE='headmotion', STOP_AFTER_MODULE='headmotion')
        self.assert_events(result, ['headmotion:topup'])

    def test_complete_medic_outputs_preserve_medic_even_with_saved_topup(self):
        self.make_medic()
        self.make_topup()
        self.assert_events(self.run_pipeline(), ['coreg_medic:medic', 'qc:medic'])

    def test_complete_medic_does_not_require_raw_phase_companions(self):
        self.make_medic()
        result = self.run_pipeline(START_FROM_MODULE='headmotion', STOP_AFTER_MODULE='headmotion')
        self.assert_events(result, ['headmotion_medic:medic'])

    def test_partial_medic_is_not_silently_replaced_by_topup(self):
        self.make_medic(runs=[(1, 1)])
        self.make_topup()
        self.assert_preflight_failure(self.run_pipeline())

    def test_empty_medic_output_is_not_counted_as_complete(self):
        self.make_medic(empty=True)
        self.make_topup()
        self.assert_preflight_failure(self.run_pipeline())

    def test_no_processed_outputs_fails_before_scientific_stages(self):
        self.assert_preflight_failure(self.run_pipeline())

    def test_missing_topup_average_is_not_reusable(self):
        self.make_topup(missing='Avg_FM_mag_acpc_brain.nii.gz')
        self.assert_preflight_failure(self.run_pipeline())

    def test_none_placeholders_are_not_treated_as_topup(self):
        self.make_topup(mode='none')
        self.assert_preflight_failure(self.run_pipeline())

    def test_direct_b0_outputs_are_not_relabelled_topup(self):
        self.make_topup(mode='direct_b0')
        self.assert_preflight_failure(self.run_pipeline())

    def test_start_session_only_requires_selected_medic_runs(self):
        self.make_medic(runs=[pair for pair in self.runs if pair[0] >= 10])
        self.assert_events(self.run_pipeline(START_SESSION=10), ['coreg_medic:medic', 'qc:medic'])

    def test_custom_task_and_fieldmap_directory(self):
        self.task = 'taskA'
        self.fm_rel = 'func/custom fieldmaps'
        self.make_runs()
        self.make_topup()
        self.assert_events(self.run_pipeline(), ['coreg:topup', 'qc:topup'])

    def test_fresh_fieldmaps_stage_still_prefers_raw_medic_phases(self):
        self.make_topup()
        self.make_raw_phases()
        result = self.run_pipeline(START_FROM_MODULE='fieldmaps', STOP_AFTER_MODULE='fieldmaps')
        self.assert_events(result, ['medic:medic'])

    def test_fresh_fieldmaps_stage_still_falls_back_to_raw_topup(self):
        raw = self.subject / 'func/unprocessed/rest/field_maps'
        raw.mkdir()
        for name in ('AP_S1_R1.nii.gz', 'PA_S1_R1.nii.gz'):
            (raw / name).write_bytes(b'fixture')
        result = self.run_pipeline(START_FROM_MODULE='fieldmaps', STOP_AFTER_MODULE='fieldmaps')
        self.assert_events(result, ['fieldmaps:topup'])

    def test_later_denoising_resume_does_not_require_distortion_inputs(self):
        result = self.run_pipeline(START_FROM_MODULE='denoise', STOP_AFTER_MODULE='denoise')
        self.assert_events(result, ['denoise:medic'])

    def test_stop_before_fieldmaps_does_not_require_distortion_inputs(self):
        result = self.run_pipeline(START_FROM_MODULE='validate', STOP_AFTER_MODULE='validate')
        self.assert_events(result, ['validate:medic'])


if __name__ == '__main__':
    unittest.main()
