"""Synthetic task recovery, timing gates and immutable output reuse."""
from pathlib import Path
import importlib.util,json,tempfile,unittest
import nibabel as nib
import numpy as np
import pandas as pd
import yaml
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('task_model',ROOT/'lib/mefmri_task_model.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class TaskModelTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
  self.subject=Path(self.tmp.name)/'Subject with spaces';self.subject.mkdir()
  self.raw=self.subject/'func/unprocessed/motor/session_1/run_1';self.raw.mkdir(parents=True)
  self.processed=self.subject/'func/motor/session_1/run_1';self.processed.mkdir(parents=True)
  self.recipe=m.read_recipe(ROOT/'config/task_models/motor.yaml');self.recipe['render']['enabled']=False
  self.registry=Path(self.tmp.name)/'recipes';self.registry.mkdir();self.recipe_path=self.registry/'motor.yaml';self.recipe_path.write_text(yaml.safe_dump(self.recipe))
  self.events=self.raw/'Motor_S1_R1_events.csv'
  pd.DataFrame([dict(onset=10+i*25,duration=10,trial_type=['thumb','index'][i%2],phase='experimental',hand='right')for i in range(12)]).to_csv(self.events,index=False)
  self.timing=self.raw/'Motor_S1_R1_taskmodel.json';self.meta=dict(alignment_status='assumed_authorized',onset_offset_seconds=0,slice_reference='pipeline_slice_timing_provenance');self.timing.write_text(json.dumps(self.meta))
 def test_pending_alignment_refuses_fit(self):
  with self.assertRaisesRegex(ValueError,'alignment is pending'):m.load_events(self.events,self.recipe,dict(self.meta,alignment_status='pending'),220,1.6)
 def test_offset_and_removals_preserve_event_clock(self):
  ev,_=m.load_events(self.events,self.recipe,dict(self.meta,onset_offset_seconds=1.6),220,1.6,removed=2);self.assertAlmostEqual(ev.onset.iloc[0],8.4)
 def test_out_of_scan_events_fail(self):
  with self.assertRaisesRegex(ValueError,'beyond BOLD'):m.load_events(self.events,self.recipe,self.meta,100,1.6)
 def test_explicit_unsupported_task_fails(self):
  with self.assertRaisesRegex(ValueError,'No supported'):m.discover(self.subject,self.registry,['rest'],1)
  self.assertEqual(len(m.discover(self.subject,self.registry,['motor'],1)),1)
 def test_rank_group_reference_not_half_tr(self):
  times=np.repeat(np.arange(15)*.105,4);j=self.raw/'timing.json';j.write_text(json.dumps({'SliceTiming':times.tolist()}))
  p=self.processed/'SliceTiming.txt';np.savetxt(p,(8-np.repeat(np.arange(1,16),4))/15)
  reference,_=m.slice_reference(j,p,1.6);self.assertAlmostEqual(reference,.735);self.assertNotAlmostEqual(reference,.8)
 def test_synthetic_effect_and_hash_checked_resume(self):
  n=220;tr=1.6;rng=np.random.default_rng(4);events,_=m.load_events(self.events,self.recipe,self.meta,n,tr)
  motion=rng.normal(0,.00001,(n,6));x,fd,qc=m.build_design(n,tr,events,motion,rng.normal(size=(n,2)),self.recipe,.5)
  self.assertEqual(qc['spike_frames'],[]);beta=np.zeros((len(x.columns),12));beta[x.columns.get_loc('thumb'),:6]=3;beta[x.columns.get_loc('index'),6:]=3;y=x.values@beta+rng.normal(0,.08,(n,12))
  axis=nib.cifti2.BrainModelAxis.from_surface(np.arange(6),6,'CortexLeft')+nib.cifti2.BrainModelAxis.from_surface(np.arange(6),6,'CortexRight');header=nib.Cifti2Header.from_axes((nib.cifti2.SeriesAxis(0,tr,n),axis))
  seg=np.zeros((12,12,12),np.int16);seg[:6,:,:]=2;seg[6:,:,:]=4;anatomy=self.subject/'anat/T1w';anatomy.mkdir(parents=True);nib.save(nib.Nifti1Image(seg,np.eye(4)),anatomy/'aparc+aseg.nii.gz')
  raw=self.raw/'Motor_S1_R1_E1.nii.gz';image=nib.Nifti1Image(np.zeros((2,2,2,n),np.float32),np.eye(4));image.header.set_zooms((1,1,1,tr));nib.save(image,raw)
  (self.raw/'Motor_S1_R1_E1.json').write_text(json.dumps({'SliceTiming':[0,tr/3,2*tr/3]}));np.savetxt(self.processed/'SliceTiming.txt',[1/3,0,-1/3]);np.savetxt(self.processed/'MCF.par',motion)
  for tag in self.recipe['input_tags']:
   nib.save(nib.Cifti2Image(y.astype(np.float32),header),self.processed/('Motor_'+tag+'.dtseries.nii'))
   image=nib.Nifti1Image(rng.normal(100,.5,(12,12,12,n)).astype(np.float32),np.eye(4));image.header.set_zooms((1,1,1,tr));nib.save(image,self.processed/('Motor_'+tag+'.nii.gz'))
  job=m.discover(self.subject,self.registry,['motor'],1)[0];folder=m.fit_run(self.subject,job,jobs=1);data=np.asarray(nib.load(folder/'OCME/thumb_minus_index.dscalar.nii').dataobj)
  self.assertTrue((data[0,:6]>2.8).all());self.assertTrue((data[0,6:]< -2.8).all());self.assertTrue((data[5]<.05).all());self.assertEqual(m.fit_run(self.subject,job,jobs=1),folder)
  (folder/'OCME/design.tsv').write_text('tampered')
  with self.assertRaisesRegex(ValueError,'incomplete or changed'):m.fit_run(self.subject,job,jobs=1)
if __name__=='__main__':unittest.main()
