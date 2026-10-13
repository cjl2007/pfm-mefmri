"""Exercise copied configs through real launchers without scientific computation."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class PortableConfigTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='portable config ')
        self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name)
        self.subject = self.base/'study space'/'SUB001'
        self.subject.mkdir(parents=True)
        self.config = self.base/'study config.sh'
        self.capture = self.base/'capture.json'
        self.env = {k: v for k, v in os.environ.items()
                    if not k.startswith(('PFM_', 'MEFMRI_', 'PIPELINE_', 'TANS_', 'IMPORT_'))}
        self.env['CAPTURE'] = str(self.capture)

    def test_copied_processing_config_keeps_repository_resources(self):
        setup = self.base/'environment.sh'
        setup.write_text('tkregister() { return 0; }\n')
        stub = self.base/'validate.sh'
        stub.write_text(f'#!{sys.executable}\nimport json, os\nfrom pathlib import Path\n'
                        'Path(os.environ["CAPTURE"]).write_text(json.dumps({k:os.environ[k] '
                        'for k in ("MEDIR", "AtlasTemplate", "NSI_EXTERNAL_ROOT")}))\n')
        # The runner executes validation with bash, so invoke the Python stub explicitly.
        wrapper = self.base/'validation wrapper.sh'
        wrapper.write_text(f'{shlex.quote(sys.executable)} {shlex.quote(str(stub))}\n')
        self.config.write_text((ROOT/'config/mefmri_wrapper_config.sh').read_text() + '\n' +
            '\n'.join(f'{k}={shlex.quote(str(v))}' for k,v in {
                'EnvironmentScript':setup, 'START_FROM_MODULE':'validate',
                'STOP_AFTER_MODULE':'validate', 'VALIDATE_MODULE':wrapper,
                'RUN_CONFIG_SNAPSHOT':'0', 'PIPELINE_QUIET_MODULE_OUTPUT':'0',
                'DISTORTION_CORRECTION_MODE':'none', 'FUNC_NOFIELDMAP_MODE':'1',
            }.items())+'\n')
        result = subprocess.run(['bash', str(ROOT/'bin/mefmri_pipeline.sh'),
                                 str(self.subject), str(self.config)], env=self.env,
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        values = json.loads(self.capture.read_text())
        self.assertEqual(values['MEDIR'], str(ROOT))
        self.assertEqual(values['AtlasTemplate'], str(ROOT/'res0urces/MNI152_T1_2mm.nii.gz'))
        self.assertEqual(values['NSI_EXTERNAL_ROOT'], str(ROOT/'lib/pfm-nsi'))

    def test_unconfigured_import_stops_before_converter(self):
        self.config.write_text((ROOT/'config/mefmri_import_config.sh').read_text())
        raw = self.base/'raw'; raw.mkdir()
        result = subprocess.run(['bash', str(ROOT/'bin/mefmri_import_raw.sh'), str(raw),
                                 str(self.base/'new subject'), str(self.config), '--session','1', '--dry-run'],
                                env=self.env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 2, result.stdout+result.stderr)
        self.assertIn('IMPORT_CONFIG_READY=1', result.stderr)
        self.assertFalse((self.base/'new subject').exists())

    def test_copied_import_config_keeps_backend_in_repository(self):
        text = (ROOT/'config/mefmri_import_config.sh').read_text().replace('IMPORT_CONFIG_READY=0', 'IMPORT_CONFIG_READY=1')
        for key,value in {'TASK_RUNS_PER_SESSION':1,'ECHOES_PER_RUN':3,
                          'SBREF_PER_RUN':0,'FMAP_AP_PER_SESSION':0,'FMAP_PA_PER_SESSION':0}.items():
            text = text.replace(f'IMPORT_EXPECT_{key}=""',f'IMPORT_EXPECT_{key}={value}')
        text = text.replace("IMPORT_TASK_REGEX='(?!)'", "IMPORT_TASK_REGEX='^EXAMPLE_REST$'")
        self.config.write_text(text)
        raw = self.base/'raw'; raw.mkdir()
        stub = self.base/'python stub'
        stub.write_text(f'#!{sys.executable}\nimport json, os, sys\nfrom pathlib import Path\n'
                        'Path(os.environ["CAPTURE"]).write_text(json.dumps(sys.argv[1:]))\n')
        stub.chmod(0o755)
        self.env['PIPELINE_PYTHON'] = str(stub)
        result = subprocess.run(['bash', str(ROOT/'bin/mefmri_import_raw.sh'),str(raw),
                                 str(self.base/'new subject'),str(self.config),'--session','1','--dry-run'],
                                env=self.env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        args = json.loads(self.capture.read_text())
        self.assertEqual(args[0], str(ROOT/'lib/mefmri_raw_import.py'))
        self.assertEqual(args[args.index('--config-file')+1], str(self.config))

    def test_copied_container_config_uses_installation_template(self):
        self.config.write_text((ROOT/'config/mefmri_wrapper_config_container.sh').read_text())
        result = subprocess.run(['bash','-eu','-c',
            'PFM_MEFMRI_ROOT="$1"; source "$2"; printf "%s\\n" "$MEDIR" "$TEDANA_ACTIVATE_MODE"',
            'test',str(ROOT),str(self.config)],env=self.env,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(result.stdout.splitlines(),[str(ROOT),'direct'])


if __name__ == '__main__':
    unittest.main()
