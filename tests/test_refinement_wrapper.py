"""Real wrapper/refinement orchestration; baseline coreg and BOLD stages are stubs."""
import json
import os
from pathlib import Path
import shutil
import sys
import unittest
import nibabel as nib
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
os.environ['MEFMRI_TEST_ROOT'] = str(ROOT)
sys.path.insert(0, str(ROOT / 'tests'))
import test_distortion_resume
sys.path.insert(0, str(ROOT / 'lib'))
from mefmri_refinement_store import digest
from mefmri_refinement_consumer import selection, mark_completed
from mefmri_refinement_stage import check_downstream


class WrapperTests(unittest.TestCase):
    def test_changed_bold_regenerates_denoising_and_blocks_late_resume(self):
        fixture = test_distortion_resume.DistortionResumeTests('runTest')
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        run = fixture.subject / 'func/rest/session_1/run_1'
        for name in ('target', 'warp'):
            (run / name).write_bytes(b'transform fixture')
        (run / 'IntermediateCoregTarget.txt').write_text(str(run / 'target'))
        (run / 'Intermediate2ACPCWarp.txt').write_text(str(run / 'warp'))
        bold = run / 'Rest_E1_acpc.nii.gz'
        bold.write_bytes(b'generation one; routing fixture only')
        selected = selection(run, 'independent')
        mark_completed(run, selected, 'Rest')
        fixture.stubs['denoise'].write_text(
            '#!/bin/bash\nset -e\nprintf "skip=%s\\n" "${MEICA_SKIP_TEDANA_IF_EXISTS:-unset}" >> "$EVENTS"\n')
        settings = dict(START_FROM_MODULE='denoise', STOP_AFTER_MODULE='denoise',
                        COREG_ALIGNMENT_STRATEGY='independent', MEICA_SKIP_TEDANA_IF_EXISTS=1)
        # All legacy runs need baseline pointers for the BOLD-generation guard.
        for other in (fixture.subject / 'func/rest').glob('session_*/run_*'):
            for pointer in ('IntermediateCoregTarget.txt', 'Intermediate2ACPCWarp.txt'):
                if other != run:
                    (other / pointer).write_text((run / pointer).read_text())
        first = fixture.run_pipeline(**settings)
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        self.assertEqual(fixture.events.read_text().splitlines(), ['skip=0'])
        second = fixture.run_pipeline(**settings)
        self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
        self.assertEqual(fixture.events.read_text().splitlines(), ['skip=0', 'skip=1'])
        bold.write_bytes(b'generation two; routing fixture only')
        mark_completed(run, selected, 'Rest')
        failed = fixture.run_pipeline(**dict(settings, START_FROM_MODULE='vol2surf', STOP_AFTER_MODULE='vol2surf'))
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn('Stale denoise', failed.stderr)
        self.assertEqual(fixture.events.read_text().splitlines(), ['skip=0', 'skip=1'])
        repaired = fixture.run_pipeline(**settings)
        self.assertEqual(repaired.returncode, 0, repaired.stdout + repaired.stderr)
        self.assertEqual(fixture.events.read_text().splitlines(), ['skip=0', 'skip=1', 'skip=0'])

    def test_rebuild_frozen_and_resume(self):
        fixture = test_distortion_resume.DistortionResumeTests('runTest')
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        shutil.rmtree(fixture.subject / 'func')  # Test-owned synthetic tree only.
        fixture.runs = [(1, 1), (2, 1)]
        fixture.make_runs()
        fixture.make_topup()
        x, y, z = np.indices((20, 22, 24), dtype=float)
        data = (50*np.exp(-((x-9)**2+(y-10)**2+(z-12)**2)/40)
                + 20*np.exp(-((x-6)**2+(y-14)**2+(z-8)**2)/9)).astype(np.float32)
        image = nib.Nifti1Image(data, np.diag([-2., 2., 2., 1.]))
        grid = fixture.base / 'grid.nii.gz'
        nib.save(image, grid)
        anatomy = fixture.subject / 'anat/T1w/T1w_acpc_dc_restore_brain.nii.gz'
        nib.save(image, anatomy)
        coreg = fixture.base / 'coreg_fixture.py'
        coreg.write_text('''import os,sys,shutil,subprocess
from pathlib import Path
import numpy as np
subject=Path(sys.argv[1]); grid=Path(sys.argv[2]); start=int(sys.argv[3])
xfms=subject/'func/xfms'/os.environ['FUNC_XFMS_DIRNAME']; xfms.mkdir(parents=True,exist_ok=True)
shutil.copyfile(grid,xfms/'T1w_acpc_brain_func.nii.gz')
np.savetxt(xfms/'identity.mat',np.eye(4))
warp=xfms/'baseline.nii.gz'
subprocess.run(['/usr/local/fsl/bin/convertwarp','--ref='+str(grid),'--premat='+str(xfms/'identity.mat'),'--out='+str(warp),'--relout'],check=True)
for run in (subject/'func/rest').glob('session_*/run_*'):
 if int(run.parent.name.split('_')[1])<start: continue
 native=run/'SBref.nii.gz';shutil.copyfile(grid,native)
 (run/'TR.txt').write_text('1.6\\n')
 (run/'PE.txt').write_text('-y\\n')
 (run/'IntermediateCoregTarget.txt').write_text(str(native)+'\\n')
 (run/'Intermediate2ACPCWarp.txt').write_text(str(warp)+'\\n')
''')
        fixture.stubs['coreg'].write_text('#!/bin/bash\nset -e\n' +
            '"'+sys.executable+'" "'+str(coreg)+'" "$3/$2" "$4" "$7"\n')
        head = fixture.base / 'head_fixture.py'
        head.write_text('''import os,sys
from pathlib import Path
sys.path.insert(0,sys.argv[1]+'/lib')
from mefmri_refinement_consumer import selection,mark_completed
subject=Path(sys.argv[2]); start=int(sys.argv[3])
for run in (subject/'func/rest').glob('session_*/run_*'):
 if int(run.parent.name.split('_')[1])<start: continue
 result=selection(run,os.environ['COREG_ALIGNMENT_STRATEGY'],os.environ['COREG_ALIGNMENT_BUNDLE'],os.environ['COREG_REFERENCE_VERSION'])
 (run/'Rest_E1_acpc.nii.gz').write_bytes(b'BOLD stage stub, routing test only')
 mark_completed(run,result,'Rest')
''')
        fixture.stubs['headmotion'].write_text('#!/bin/bash\nset -e\n' +
            '"'+sys.executable+'" "'+str(head)+'" "$1" "$3/$2" "$7"\n')
        settings = dict(DISTORTION_CORRECTION_MODE='topup', START_FROM_MODULE='coreg',
                        STOP_AFTER_MODULE='headmotion', AtlasTemplate=grid,
                        COREG_ALIGNMENT_STRATEGY='subject_refined', COREG_REFERENCE_POLICY='rebuild')
        result = fixture.run_pipeline(**settings)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        reference_root = fixture.subject / 'func/xfms/rest/SubjectRefinement'
        bundles = list(reference_root.glob('invocation-*/alignment/alignment.json'))
        self.assertEqual(len(bundles), 1)
        first = json.loads(bundles[0].read_text())
        old_paths = list((fixture.subject/'func/xfms'/first['xfms_namespace']).rglob('*'))
        old_paths += list((fixture.subject/'func/rest/session_1').rglob('*'))
        before = {str(p): digest(p) for p in old_paths if p.is_file()}
        fixture.runs = [(3, 1)]
        fixture.make_runs()
        settings.update(START_SESSION=3, COREG_REFERENCE_POLICY='frozen', COREG_REFERENCE_VERSION=first['version'])
        result = fixture.run_pipeline(**settings)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(before, {p: digest(p) for p in before})
        bundles = list(reference_root.glob('invocation-*/alignment/alignment.json'))
        second_path = next(p for p in bundles if json.loads(p.read_text())['policy'] == 'frozen')
        second = json.loads(second_path.read_text())
        self.assertEqual(first['version'], second['version'])
        self.assertNotEqual(first['xfms_namespace'], second['xfms_namespace'])
        # All-session downstream reuse accepts separately processed frozen
        # batches only when they share the same reference and valid bindings.
        check_downstream(fixture.subject, 'rest', 1, 'subject_refined', second_path.parent, second['version'])
        old_record = fixture.subject/'func/rest/session_1/run_1/BOLDAlignment.json'
        saved_text = old_record.read_text()
        changed = json.loads(saved_text)
        changed['alignment']['version'] = 'different-reference'
        old_record.write_text(json.dumps(changed))
        with self.assertRaisesRegex(ValueError, 'Mixed alignment'):
            check_downstream(fixture.subject, 'rest', 1, 'subject_refined', second_path.parent, second['version'])
        old_record.write_text(saved_text)
        settings.update(START_FROM_MODULE='headmotion', COREG_ALIGNMENT_BUNDLE=second_path.parent)
        result = fixture.run_pipeline(**settings)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == '__main__':
    unittest.main()
