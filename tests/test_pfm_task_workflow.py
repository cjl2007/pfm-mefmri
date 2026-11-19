"""Exercise the real vector wrapper, concatenation and TANS shell adapter."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import unittest

import test_pfm_dataset as fixtures
ROOT = fixtures.ROOT


class TaskWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.DatasetTests('test_rest_only_excludes_tasks')
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.subject = self.fixture.subject
        self.base = self.subject.parent.parent
        self.tans_config = self.base / 'tans.yaml'
        self.tans_config.write_text('{}\n')
        self.capture = self.base / 'tans-arguments.json'
        self.events = self.base / 'events.jsonl'
        for spec in self.fixture.specs:
            for session in (1, 2):
                run = self.subject / 'func/unprocessed' / spec['task'] / ('session_' + str(session)) / 'run_1'
                run.mkdir(parents=True)
                for echo in range(1, 5):
                    (run / '{}_S{}_R1_E{}.nii.gz'.format(spec['prefix'], session, echo)).touch()
                (self.subject / 'func' / spec['task'] / ('session_' + str(session)) / 'run_1/TE.txt').write_text('0.01 0.02 0.03 0.04')
        for modality in ('T1w', 'T2w'):
            file = self.subject / 'anat' / modality / (modality + '.nii.gz')
            file.parent.mkdir(parents=True)
            file.touch()
        self.setup = self.base / 'environment.sh'
        self.setup.write_text('tkregister() { return 0; }\n')
        self.pfm = self.base / 'pfm.sh'
        self.pfm.write_text('''#!/usr/bin/env bash
"$PIPELINE_PYTHON" - "$@" <<'PY'
import json,os,sys
from pathlib import Path
root=Path(os.environ.get('PFM_OUTDIR') or (str(Path(sys.argv[2])/sys.argv[1]/'pfm')))
(root/'RidgePFM').mkdir(parents=True,exist_ok=True)
(root/'RidgePFM/RidgePFM_VTX.dlabel.nii').write_bytes(b'scientific computation stub')
with open(os.environ['EVENTS'],'a') as f:
 f.write(json.dumps({'stage':'pfm','task':sys.argv[5], 'input':os.environ['PFM_INPUT_CIFTI'],
                    'reference':os.environ['PFM_FUNC_XFMS_DIRNAME']})+'\\n')
PY
''')
        self.tans_code = self.base / 'tans code'
        self.tans_code.mkdir()
        (self.tans_code / 'run_pipeline.py').write_text('''import json,os,sys
from pathlib import Path
Path(os.environ['CAPTURE']).write_text(json.dumps(sys.argv[1:]))
with open(os.environ['EVENTS'],'a') as f: f.write(json.dumps({'stage':'tans'})+'\\n')
''')

    def invoke(self, **overrides):
        settings = dict(MEDIR=ROOT, EnvironmentScript=self.setup,
            FUNC_DIRNAME=['rest', 'effort'], FUNC_FILE_PREFIX=['Rest', 'Effort'],
            PFM_INPUT_TASKS=['rest'], START_SESSION=1, START_FROM_MODULE='concat', STOP_AFTER_MODULE='tans',
            PIPELINE_SKIP_ANAT=1, PIPELINE_PYTHON=sys.executable, PFM_PYTHON=sys.executable,
            CONCAT_PYTHON=sys.executable, CONCAT_INPUT_TAG='Tag',
            DISTORTION_CORRECTION_MODE='none', FUNC_NOFIELDMAP_MODE=1,
            PFM_ENABLE=1, PFM_RF_CONFIDENCE_ENABLE=0, FUNC_PFM_MODULE=self.pfm,
            TANS_ENABLE=1, TANS_CODE_ROOT=self.tans_code,
            TANS_CONFIG_FILE=self.tans_config,
            TANS_PYTHON=sys.executable, TANS_INCLUDE_SM_HAND_CALIBRATION=0,
            NSI_ENABLE=0, CONCAT_ENABLE=1, CRAWLING_SEED_FC_ENABLE=0,
            RUN_CONFIG_SNAPSHOT=0, PIPELINE_QUIET_MODULE_OUTPUT=0)
        settings.update(overrides)
        config = self.base / 'workflow config.sh'
        def shell(value):
            if isinstance(value, list):
                return '(' + ' '.join(shlex.quote(str(v)) for v in value) + ')'
            return shlex.quote(str(value))
        config.write_text('\n'.join(k+'='+shell(v) for k,v in settings.items())+'\n')
        environment = {k:v for k,v in os.environ.items() if not k.startswith(
            ('PIPELINE_', 'PFM_', 'TANS_', 'FUNC_', 'CONCAT_', 'NSI_')) and k not in ('BASH_ENV','ENV')}
        environment.update(EVENTS=str(self.events), CAPTURE=str(self.capture))
        return subprocess.run(['bash', str(ROOT/'bin/mefmri_pipeline.sh'), str(self.subject), str(config)],
                              env=environment, text=True, capture_output=True, timeout=40)

    def assert_success(self, result):
        self.assertEqual(result.returncode, 0, result.stdout+'\n'+result.stderr)
        events = [json.loads(line) for line in self.events.read_text().splitlines()]
        self.assertEqual([row['stage'] for row in events], ['pfm','tans'])
        manifest = json.loads((self.subject/'pfm/dataset.json').read_text())
        self.assertEqual(manifest['status'], 'complete')
        arguments = json.loads(self.capture.read_text())
        self.assertEqual(arguments[arguments.index('--pfm-dir')+1], str(self.subject/'pfm/RidgePFM'))
        return events, manifest, arguments

    def test_process_all_map_rest_only_once_and_tans_uses_new_root(self):
        events, manifest, args = self.assert_success(self.invoke())
        self.assertEqual([s['task'] for s in manifest['tasks']], ['rest'])
        self.assertEqual(manifest['n_frames'], 12)
        self.assertIn('/pfm/ConcatenatedCiftis/', events[0]['input'])
        self.assertEqual(events[0]['reference'], 'rest')
        for task in ('rest', 'effort'):
            self.assertTrue((self.subject/'func'/task/'ConcatenatedCiftis/PFMPrepSummary.json').exists())

    def test_pool_requested_tasks_and_keep_reference_separate_from_output(self):
        events, manifest, args = self.assert_success(self.invoke(PFM_INPUT_TASKS=['Rest','Effort']))
        self.assertEqual(manifest['n_frames'], 24)
        self.assertEqual([s['task'] for s in manifest['tasks']], ['rest','effort'])
        self.assertEqual(args[args.index('--func-dir')+1], 'rest')

    def test_effort_only_pfm_reference_is_effort(self):
        events, manifest, args = self.assert_success(self.invoke(PFM_INPUT_TASKS=['Effort']))
        self.assertEqual(events[0]['reference'], 'effort')
        self.assertEqual(args[args.index('--func-dir')+1], 'effort')

    def test_arrays_must_have_matching_lengths_before_processing(self):
        result = self.invoke(FUNC_FILE_PREFIX=['Rest'])
        self.assertEqual(result.returncode, 2)
        self.assertIn('same length', result.stderr)
        self.assertFalse(self.events.exists())

    def test_two_sessions_from_multiple_tasks_do_not_become_four_sessions(self):
        result = self.invoke(PFM_INPUT_TASKS=['rest','effort'], PFM_RF_CONFIDENCE_ENABLE=1)
        self.assertEqual(result.returncode, 2, result.stdout+result.stderr)
        self.assertIn('found 2 raw sessions', result.stderr)
        self.assertFalse(self.events.exists())

    def test_tans_only_resume_does_not_rerun_mapping(self):
        self.assert_success(self.invoke())
        self.events.unlink()
        result = self.invoke(START_FROM_MODULE='tans')
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertEqual([json.loads(l)['stage'] for l in self.events.read_text().splitlines()], ['tans'])

    def test_tans_resume_cannot_silently_change_input_selection(self):
        self.assert_success(self.invoke())
        self.events.unlink()
        result = self.invoke(START_FROM_MODULE='tans', PFM_INPUT_TASKS=['rest','effort'])
        self.assertEqual(result.returncode, 2, result.stdout+result.stderr)
        self.assertIn('do not match', result.stderr)
        self.assertFalse(self.events.exists())

    def test_explicit_cifti_override_is_honored_and_recorded(self):
        source = self.subject / 'func/rest/session_1/run_1/Rest_Tag.dtseries.nii'
        events, manifest, args = self.assert_success(self.invoke(START_FROM_MODULE='pfm', PFM_INPUT_CIFTI=source))
        self.assertEqual(events[0]['input'], str(source))
        self.assertEqual(manifest['input_mode'], 'external_cifti_override')
        self.assertFalse(manifest['task_membership_verified'])

    def test_custom_pfm_root_is_also_the_tans_input(self):
        root = self.subject / 'custom mapping root'
        result = self.invoke(PFM_OUTDIR=root)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        arguments = json.loads(self.capture.read_text())
        self.assertEqual(arguments[arguments.index('--pfm-dir')+1], str(root / 'RidgePFM'))
        self.assertEqual(json.loads((root / 'dataset.json').read_text())['status'], 'complete')

    def test_all_inputs_are_validated_before_shared_anatomy(self):
        stubs = {}
        for stage in ('validate','anat_hcp','anat_charm'):
            path = self.base / (stage + '.sh')
            path.write_text('printf \'{"stage":"%s","task":"%s"}\\n\' ' + stage + ' "$FUNC_DIRNAME" >> "$EVENTS"\n')
            stubs[stage] = path
        license = self.base / 'fixture-license.txt'
        license.touch()
        result = self.invoke(START_FROM_MODULE='validate', STOP_AFTER_MODULE='anat_charm',
            PIPELINE_SKIP_ANAT=0, FS_LICENSE_FILE=license,
            VALIDATE_MODULE=stubs['validate'], ANAT_HCP_MODULE=stubs['anat_hcp'], ANAT_CHARM_MODULE=stubs['anat_charm'])
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        rows = [json.loads(l) for l in self.events.read_text().splitlines()]
        self.assertEqual([r['stage'] for r in rows], ['validate','validate','anat_hcp','anat_charm'])
        self.assertEqual([r['task'] for r in rows[:2]], ['rest','effort'])


    def qc_stubs(self):
        modules = {}
        for stage in ('nsi', 'crawling_seed_fc'):
            path = self.base / (stage + '.sh')
            path.write_text("\"$PIPELINE_PYTHON\" - <<'PYQC'\n" +
                "import os,json\nfrom pathlib import Path\nimport nibabel as nib\n" +
                "stage=" + repr(stage) + "\n" +
                "source=os.environ['NSI_INPUT_CIFTI' if stage=='nsi' else 'CRAWLING_SEED_FC_INPUT_CIFTI']\n" +
                "out=Path(os.environ['NSI_OUTDIR' if stage=='nsi' else 'CRAWLING_SEED_FC_OUTDIR']);out.mkdir(parents=True,exist_ok=True)\n" +
                "row={'stage':stage,'input':source,'frames':nib.load(source).shape[0],'out':str(out)}\n" +
                "(out/'result.json').write_text(json.dumps(row))\n" +
                "with open(os.environ['EVENTS'],'a') as f:f.write(json.dumps(row)+'\\n')\nPYQC\n")
            modules['FUNC_NSI_MODULE' if stage == 'nsi' else 'FUNC_CRAWLING_SEED_FC_MODULE'] = path
        return dict(NSI_ENABLE=1, CRAWLING_SEED_FC_ENABLE=1, **modules)

    def test_pooled_qc_uses_selected_input_once_before_mapping(self):
        result = self.invoke(PFM_INPUT_TASKS=['rest','effort'], **self.qc_stubs())
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        rows = [json.loads(x) for x in self.events.read_text().splitlines()]
        self.assertEqual([x['stage'] for x in rows], ['crawling_seed_fc','nsi','pfm','tans'])
        self.assertEqual([x['frames'] for x in rows[:2]], [16,16])
        self.assertTrue(all(x['input'] == rows[2]['input'] for x in rows[:2]))
        self.assertEqual(rows[1]['out'], str(self.subject/'pfm/QC/NSI'))
        self.assertFalse((self.subject/'func/rest/QC/NSI').exists())
        self.assertIn('pfm_prepare', result.stdout)

    def test_qc_only_and_scalar_selection_use_rest_dataset(self):
        result = self.invoke(FUNC_DIRNAME='rest', FUNC_FILE_PREFIX='Rest',
            PFM_ENABLE=0, TANS_ENABLE=0, STOP_AFTER_MODULE='nsi', **self.qc_stubs())
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        rows = [json.loads(x) for x in self.events.read_text().splitlines()]
        self.assertEqual([x['stage'] for x in rows], ['crawling_seed_fc','nsi'])
        self.assertEqual([x['frames'] for x in rows], [8,8])
        self.assertFalse((self.subject/'pfm/RidgePFM').exists())

    def test_qc_receipts_reuse_and_missing_outputs_regenerate(self):
        settings = dict(STOP_AFTER_MODULE='nsi', **self.qc_stubs())
        first = self.invoke(**settings)
        self.assertEqual(first.returncode, 0, first.stdout+first.stderr)
        self.events.unlink()
        again = self.invoke(START_FROM_MODULE='pfm_prepare', **settings)
        self.assertEqual(again.returncode, 0, again.stdout+again.stderr)
        self.assertFalse(self.events.exists())
        self.assertIn('reusing matching nsi', again.stdout)
        (self.subject/'pfm/QC/CrawlingSeedFC/result.json').unlink()
        missing = self.invoke(START_FROM_MODULE='pfm_prepare', **settings)
        self.assertEqual(missing.returncode, 0, missing.stdout+missing.stderr)
        self.assertEqual([json.loads(x)['stage'] for x in self.events.read_text().splitlines()], ['crawling_seed_fc'])

    def test_resume_nsi_rejects_different_dataset_selection(self):
        settings = dict(STOP_AFTER_MODULE='nsi', **self.qc_stubs())
        first = self.invoke(**settings)
        self.assertEqual(first.returncode, 0, first.stdout+first.stderr)
        self.events.unlink()
        changed = self.invoke(START_FROM_MODULE='nsi', PFM_INPUT_TASKS=['rest','effort'], **settings)
        self.assertNotEqual(changed.returncode, 0)
        self.assertIn('input selection or source files changed', changed.stderr)
        self.assertFalse(self.events.exists())

    def test_prepare_stop_and_custom_qc_root(self):
        settings = self.qc_stubs()
        root = self.subject/'alternative PFM root'
        first = self.invoke(STOP_AFTER_MODULE='pfm_prepare', PFM_OUTDIR=root, **settings)
        self.assertEqual(first.returncode, 0, first.stdout+first.stderr)
        self.assertFalse(self.events.exists())
        second = self.invoke(START_FROM_MODULE='crawling_seed_fc', STOP_AFTER_MODULE='crawling_seed_fc', PFM_OUTDIR=root, **settings)
        self.assertEqual(second.returncode, 0, second.stdout+second.stderr)
        rows = [json.loads(x) for x in self.events.read_text().splitlines()]
        self.assertEqual([x['stage'] for x in rows], ['crawling_seed_fc'])
        self.assertEqual(rows[0]['out'], str(root/'QC/CrawlingSeedFC'))

    def test_nsi_setting_change_invalidates_nsi_receipt(self):
        settings = dict(STOP_AFTER_MODULE='nsi', **self.qc_stubs())
        first = self.invoke(NSI_EXTERNAL_MORANS=0, **settings)
        self.assertEqual(first.returncode, 0, first.stdout+first.stderr)
        self.events.unlink()
        second = self.invoke(START_FROM_MODULE='nsi', NSI_EXTERNAL_MORANS=1, **settings)
        self.assertEqual(second.returncode, 0, second.stdout+second.stderr)
        self.assertEqual([json.loads(x)['stage'] for x in self.events.read_text().splitlines()], ['nsi'])


    def test_concat_stop_does_not_run_selected_dataset_qc(self):
        result = self.invoke(STOP_AFTER_MODULE='concat', **self.qc_stubs())
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertFalse(self.events.exists())
        self.assertFalse((self.subject/'pfm/dataset.json').exists())

    def test_prepare_can_run_when_mapping_and_qc_are_disabled(self):
        result = self.invoke(PFM_ENABLE=0, TANS_ENABLE=0, NSI_ENABLE=0,
            CRAWLING_SEED_FC_ENABLE=0, STOP_AFTER_MODULE='pfm_prepare')
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertTrue((self.subject/'pfm/dataset.json').is_file())
        self.assertFalse(self.events.exists())

    def test_external_cifti_qc_uses_override_without_scaninfo(self):
        source = self.subject/'func/effort/session_1/run_1/Effort_Tag.dtseries.nii'
        result = self.invoke(START_FROM_MODULE='pfm_prepare', STOP_AFTER_MODULE='nsi',
            PFM_INPUT_CIFTI=source, **self.qc_stubs())
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        rows = [json.loads(x) for x in self.events.read_text().splitlines()]
        self.assertEqual([x['input'] for x in rows], [str(source),str(source)])
        self.assertEqual([x['frames'] for x in rows], [6,6])


    def test_legacy_tans_resume_resolves_root_before_creating_logs(self):
        self.assert_success(self.invoke())
        legacy = self.subject/'func/rest/PFM'
        (self.subject/'pfm').rename(legacy)
        self.events.unlink()
        result = self.invoke(START_FROM_MODULE='tans')
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        args = json.loads(self.capture.read_text())
        self.assertEqual(args[args.index('--pfm-dir')+1], str(legacy/'RidgePFM'))
        self.assertTrue((legacy/'QC/ModuleLogs').is_dir())


    def test_nsi_stop_is_honored_when_nsi_disabled(self):
        result = self.invoke(STOP_AFTER_MODULE='nsi', NSI_ENABLE=0, CRAWLING_SEED_FC_ENABLE=0)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertFalse(self.events.exists())
        self.assertTrue((self.subject/'pfm/dataset.json').is_file())


if __name__ == '__main__':
    unittest.main(verbosity=2)
