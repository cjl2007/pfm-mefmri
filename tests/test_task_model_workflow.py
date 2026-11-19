"""Real wrapper ordering with preprocessing/PFM/TANS/model computations stubbed."""
import json
from pathlib import Path
import unittest
import test_pfm_task_workflow as fixtures

class TaskModelWorkflowTests(unittest.TestCase):
 def setUp(self):
  self.w=fixtures.TaskWorkflowTests('test_pool_requested_tasks_and_keep_reference_separate_from_output');self.w.setUp();self.addCleanup(self.w.doCleanups)
  self.model=self.w.base/'task model.sh'
  self.model.write_text('''#!/usr/bin/env bash
"$PIPELINE_PYTHON" - <<'PYCODE'
import json,os
with open(os.environ['EVENTS'],'a') as f:f.write(json.dumps({'stage':'task_model','tasks':json.loads(os.environ['TASK_MODEL_SELECTED_TASKS_JSON'])})+'\\n')
PYCODE
''')
 def invoke(self,**opts):
  settings=dict(TASK_MODEL_ENABLE=1,FUNC_TASK_MODEL_MODULE=self.model,STOP_AFTER_MODULE='task_model',PFM_INPUT_TASKS=['rest','effort']);settings.update(opts);return self.w.invoke(**settings)
 def test_once_after_combined_pfm_and_tans(self):
  r=self.invoke();self.assertEqual(r.returncode,0,r.stdout+r.stderr)
  events=[json.loads(l)for l in self.w.events.read_text().splitlines()];self.assertEqual([e['stage']for e in events],['pfm','tans','task_model']);self.assertEqual(events[-1]['tasks'],['rest','effort'])
 def test_model_only_resume_does_not_prepare_or_map_pfm(self):
  r=self.invoke(START_FROM_MODULE='task_model',PFM_ENABLE=0,TANS_ENABLE=0,NSI_ENABLE=0,CRAWLING_SEED_FC_ENABLE=0);self.assertEqual(r.returncode,0,r.stdout+r.stderr)
  self.assertEqual([json.loads(l)['stage']for l in self.w.events.read_text().splitlines()],['task_model'])
  self.assertFalse((self.w.subject/'pfm/dataset.json').exists())
 def test_stop_at_pfm_does_not_model(self):
  r=self.invoke(STOP_AFTER_MODULE='pfm',TANS_ENABLE=0);self.assertEqual(r.returncode,0,r.stdout+r.stderr)
  self.assertEqual([json.loads(l)['stage']for l in self.w.events.read_text().splitlines()],['pfm'])
 def test_disabled_model_boundary_errors_before_processing(self):
  r=self.invoke(TASK_MODEL_ENABLE=0);self.assertEqual(r.returncode,2,r.stdout+r.stderr);self.assertIn('requires TASK_MODEL_ENABLE=1',r.stderr);self.assertFalse(self.w.events.exists())
if __name__=='__main__':unittest.main()
