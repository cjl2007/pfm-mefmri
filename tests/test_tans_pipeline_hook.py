"""Exercise the real shell runner/wrapper, with scientific computation stubbed."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class TansPipelineHookTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='tans hook ')
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        self.subject = self.base/'study with spaces'/'SUB001'
        for modality in ('T1w', 'T2w'):
            path = self.subject/'anat'/modality/f'{modality}.nii.gz'
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b'fixture')
        (self.subject/'func/rest/PFM/RidgePFM').mkdir(parents=True)
        self.events = self.base/'events.txt'
        self.capture = self.base/'args.json'
        self.env_file = self.base/'environment.sh'
        self.env_file.write_text('tkregister() { return 0; }\n')
        self.code = self.base/'tans code'
        self.code.mkdir()
        (self.code/'run_pipeline.py').write_text(
            'import os,json,sys\nfrom pathlib import Path\n'
            'Path(os.environ["ARG_CAPTURE"]).write_text(json.dumps(sys.argv[1:]))\n'
            'with open(os.environ["EVENTS"], "a") as f: f.write("tans\\n")\n')
        self.yaml = self.base/'study config.yaml'
        self.yaml.write_text('subject: "{subject}"\n')
        self.stubs = {}
        for stage in ('pfm', 'pfm_update'):
            path = self.base/f'{stage}.sh'
            path.write_text('#!/usr/bin/env bash\n'+f'printf "{stage}\\n" >> "$EVENTS"\n')
            self.stubs[stage] = path

    def run_pipeline(self, **overrides):
        settings = dict(
            MEDIR=ROOT, EnvironmentScript=self.env_file,
            START_FROM_MODULE='pfm', STOP_AFTER_MODULE='tans',
            FUNC_DIRNAME='rest', FUNC_FILE_PREFIX='Rest',
            TANS_ENABLE='1', TANS_CODE_ROOT=self.code, TANS_CONFIG_FILE=self.yaml,
            TANS_PYTHON=sys.executable, TANS_TARGET_NETWORK='CinguloOpercular',
            TANS_OUTPUT_DIRNAME='tans', TANS_INCLUDE_SM_HAND_CALIBRATION='0',
            TANS_SM_HAND_COIL_MODEL='', TANS_TREATMENT_COIL_MODEL='',
            TANS_SM_HAND_DIDT_CALIBRATION='', TANS_TREATMENT_DIDT_CALIBRATION='',
            TANS_RESUME='1', PFM_ENABLE='1', PFM_INFOMAP_UPDATE_ENABLE='1',
            FUNC_PFM_MODULE=self.stubs['pfm'], FUNC_PFM_INFOMAP_UPDATE_MODULE=self.stubs['pfm_update'],
            FUNC_TANS_MODULE=ROOT/'modules/mefmri_func_tans.sh',
            CONCAT_ENABLE='0', NSI_ENABLE='0', CRAWLING_SEED_FC_ENABLE='0',
            COREG_SDC_QC_ENABLE='0', MOTION_GRAYPLOTS_QC_ENABLE='0',
            RUN_CONFIG_SNAPSHOT='0', PIPELINE_QUIET_MODULE_OUTPUT='0',
        )
        settings.update(overrides)
        config = self.base/'pipeline config.sh'
        config.write_text('\n'.join(f'{key}={shlex.quote(str(value))}' for key,value in settings.items())+'\n')
        environment = {k:v for k,v in os.environ.items() if not k.startswith(('TANS_', 'PIPELINE_', 'PFM_'))}
        environment.update(EVENTS=str(self.events), ARG_CAPTURE=str(self.capture))
        return subprocess.run(['bash',str(ROOT/'bin/mefmri_pipeline.sh'),str(self.subject),str(config)],
                              capture_output=True, text=True, env=environment)

    def test_tans_runs_after_pfm_and_manual_update(self):
        result = self.run_pipeline()
        self.assertEqual(result.returncode, 0, result.stdout+'\n'+result.stderr)
        self.assertEqual(self.events.read_text().splitlines(), ['pfm','pfm_update','tans'])
        arguments = json.loads(self.capture.read_text())
        self.assertEqual(arguments[arguments.index('--subject-dir')+1], str(self.subject))
        self.assertEqual(arguments[arguments.index('--output')+1], str(self.subject/'tans'))
        self.assertIn('--resume', arguments)
        self.assertIn('--skip-sm-hand-calibration', arguments)

    def test_tans_only_resume_skips_upstream_modules(self):
        result = self.run_pipeline(START_FROM_MODULE='tans')
        self.assertEqual(result.returncode, 0, result.stdout+'\n'+result.stderr)
        self.assertEqual(self.events.read_text().splitlines(), ['tans'])

    def test_disabled_tans_does_not_invoke_driver(self):
        result = self.run_pipeline(TANS_ENABLE='0', STOP_AFTER_MODULE='')
        self.assertEqual(result.returncode, 0, result.stdout+'\n'+result.stderr)
        self.assertEqual(self.events.read_text().splitlines(), ['pfm','pfm_update'])
        self.assertFalse(self.capture.exists())

    def test_explicit_tans_route_requires_enable(self):
        for field in ('START_FROM_MODULE', 'STOP_AFTER_MODULE'):
            options = {'TANS_ENABLE':'0','START_FROM_MODULE':'pfm','STOP_AFTER_MODULE':''}
            options[field] = 'tans'
            result = self.run_pipeline(**options)
            self.assertEqual(result.returncode, 2, result.stdout+'\n'+result.stderr)
            self.assertIn('requires TANS_ENABLE=1', result.stderr)
            self.assertFalse(self.capture.exists())

    def test_invalid_enable_is_not_a_silent_noop(self):
        result = self.run_pipeline(TANS_ENABLE='yes')
        self.assertEqual(result.returncode, 2)
        self.assertIn('TANS_ENABLE must be 0 or 1', result.stderr)

    def test_missing_tans_config_fails_before_scientific_driver(self):
        result = self.run_pipeline(START_FROM_MODULE='tans', TANS_CONFIG_FILE=self.base/'missing.yaml')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('missing TANS config', result.stdout+result.stderr)
        self.assertFalse(self.capture.exists())


if __name__ == '__main__':
    unittest.main()
