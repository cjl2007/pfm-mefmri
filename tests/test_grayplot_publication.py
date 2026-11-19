import importlib.util,json,os,subprocess,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
ROOT=Path(os.environ.get('MEFMRI_TEST_ROOT',Path(__file__).resolve().parents[1]));sys.path.insert(0,str(ROOT/'lib'))
import mefmri_pfm_grayplots as publication
from mefmri_qc_motion_grayplots import resolve_network_labels

class GrayplotPublicationTests(unittest.TestCase):
 def fixture(self):
  t=tempfile.TemporaryDirectory();self.addCleanup(t.cleanup);p=Path(t.name);subject=p/'Subject';root=subject/'pfm';label=root/'RidgePFM/RidgePFM_VTX.dlabel.nii';label.parent.mkdir(parents=True);label.write_bytes(b'labels')
  run=subject/'func/rest/session_1/run_1';run.mkdir(parents=True);(run/'stage.dtseries.nii').write_bytes(b'input')
  out=subject/'func/rest/QC/MotionGrayPlots';out.mkdir(parents=True);(out/'motion_grayplots_qc.html').write_text('old structural report')
  env={'PFM_DATASET_ROOT':str(root),'PFM_DATASET_TASKS_JSON':json.dumps([dict(task='rest',prefix='Rest')])}
  return subject,root,label,out,env
 def renderer(self,calls,fail=False,mutate=None):
  def render(cmd,**kwargs):
   calls.append(cmd)
   if fail:raise subprocess.CalledProcessError(1,cmd)
   out=Path(cmd[cmd.index('--out-dir')+1]);out.mkdir();(out/'motion_grayplots_qc.html').write_text('network report');(out/'network_order.json').write_text('{}')
   if mutate:mutate()
  return render
 def test_replaces_structural_report_at_stable_path_and_reuses(self):
  s,r,l,out,env=self.fixture();calls=[];run=self.renderer(calls)
  a=publication.run(ROOT,s,'1',env,run);b=publication.run(ROOT,s,'1',env,run)
  self.assertEqual(len(calls),1);self.assertEqual(a[0]['report'],str(out/'motion_grayplots_qc.html'));self.assertEqual(a[0]['report'],b[0]['report']);self.assertEqual((out/'motion_grayplots_qc.html').read_text(),'network report');self.assertFalse((r/'QC/NetworkGrayPlots').exists())
  l.write_bytes(b'new labels');publication.run(ROOT,s,'1',env,run);self.assertEqual(len(calls),2)
 def test_failed_render_preserves_current_report_and_legacy(self):
  s,r,l,out,env=self.fixture();legacy=r/'QC/NetworkGrayPlots/rest/aaaaaaaaaaaaaaaa';legacy.mkdir(parents=True);(legacy/'complete.json').write_text(json.dumps(dict(label_sha256='x',files=[],code={})));(legacy/'motion_grayplots_qc.html').write_text('legacy')
  with self.assertRaises(subprocess.CalledProcessError):publication.run(ROOT,s,'1',env,self.renderer([],fail=True))
  self.assertEqual((out/'motion_grayplots_qc.html').read_text(),'old structural report');self.assertTrue(legacy.exists())
  publication.run(ROOT,s,'1',env,self.renderer([]));self.assertFalse(legacy.exists())
 def test_changed_inputs_during_render_do_not_replace_report(self):
  s,r,l,out,env=self.fixture()
  with self.assertRaises(ValueError):publication.run(ROOT,s,'1',env,self.renderer([],mutate=lambda:l.write_bytes(b'changed')))
  self.assertEqual((out/'motion_grayplots_qc.html').read_text(),'old structural report')
 def test_qc_reuses_recorded_labels_and_refuses_silent_downgrade(self):
  s,r,l,out,env=self.fixture();(out/'network_order.json').write_text(json.dumps(dict(label_file=str(l))))
  with patch.dict(os.environ,{},clear=True):self.assertEqual(resolve_network_labels(s,'rest',out),l)
  l.unlink()
  with patch.dict(os.environ,{},clear=True):
   with self.assertRaises(ValueError):resolve_network_labels(s,'rest',out)
 def test_plain_qc_without_a_fitted_map_is_available(self):
  with tempfile.TemporaryDirectory() as t,patch.dict(os.environ,{},clear=True):self.assertIsNone(resolve_network_labels(Path(t),'rest',Path(t)/'QC'))

class NetworkDisplayTests(unittest.TestCase):
 def test_colors_order_and_display_do_not_change_signal_metrics(self):
  import numpy as np,nibabel as nb
  from matplotlib.axes import Axes
  import mefmri_qc_motion_grayplots as gray
  from mefmri_network_grayplot import network_order
  with tempfile.TemporaryDirectory() as t:
   t=Path(t);brain=nb.cifti2.BrainModelAxis.from_mask(np.ones(4,dtype=bool),name='CortexLeft')
   data=np.array([[0,1,4,2],[2,2,1,0],[3,0,2,5],[1,4,0,3]],dtype=float)
   image=nb.Cifti2Image(data,header=nb.cifti2.Cifti2Header.from_axes((nb.cifti2.SeriesAxis(0,1,4),brain)))
   axis=nb.cifti2.LabelAxis(['networks'],{0:('None',(0,0,0,0)),1:('Red',(1,0,0,1)),2:('Blue',(0,0,1,1))})
   labels=nb.Cifti2Image(np.array([[2,1,0,1]],dtype=np.int32),header=nb.cifti2.Cifti2Header.from_axes((axis,brain)));path=t/'labels.dlabel.nii';nb.save(labels,path)
   info=network_order(image,np.arange(4),path)
   np.testing.assert_array_equal(info['order'],[1,3,0,2]);np.testing.assert_array_equal(info['colors'],[[1,0,0,1],[1,0,0,1],[0,0,1,1],[.7,.7,.7,1]])
   z=gray.zscore_columns(data);original=z.copy();fd=np.array([0,.1,.4,.2]);dv=gray.dvars_from_z(z);coupling=gray.grayordinate_fd_corr(z,fd)
   shown=[];imshow=Axes.imshow
   def capture(ax,a,*args,**kwargs):
    shown.append(np.asarray(a).copy());return imshow(ax,a,*args,**kwargs)
   with patch.object(Axes,'imshow',capture):
    gray.plot_motion_grayplot_panel(t/'plot.png','Example',['OCME'],{'OCME':z},fd,None,None,{'OCME':dv},network_info=info)
   self.assertTrue(any(a.shape==(4,4) and np.array_equal(a,z[:,[1,3,0,2]].T) for a in shown))
   self.assertTrue(any(a.shape==(4,1,4) and np.array_equal(a[:,0,:],info['colors']) for a in shown))
   np.testing.assert_array_equal(z,original);np.testing.assert_array_equal(gray.dvars_from_z(z),dv);np.testing.assert_array_equal(gray.grayordinate_fd_corr(z,fd),coupling)

class CompartmentDisplayTests(unittest.TestCase):
 def test_mixed_networks_are_ordered_within_anatomical_blocks(self):
  import numpy as np,nibabel as nb
  from mefmri_network_grayplot import network_order,draw_network_strip
  import matplotlib.pyplot as plt
  from matplotlib.axes import Axes
  with tempfile.TemporaryDirectory() as t:
   cortex=nb.cifti2.BrainModelAxis.from_mask(np.ones(3,dtype=bool),name='CortexLeft')
   volume=nb.cifti2.BrainModelAxis.from_mask(np.ones((1,1,3),dtype=bool),name='ThalamusLeft',affine=np.eye(4));brain=cortex+volume
   image=nb.Cifti2Image(np.ones((4,6)),header=nb.cifti2.Cifti2Header.from_axes((nb.cifti2.SeriesAxis(0,1,4),brain)))
   table={0:('Unassigned',(0,0,0,0)),1:('Red',(1,0,0,1)),2:('Green',(0,1,0,1))};labels=nb.Cifti2Image(np.array([[2,1,0,1,0,2]],dtype=np.int32),header=nb.cifti2.Cifti2Header.from_axes((nb.cifti2.LabelAxis(['networks'],table),brain)));path=Path(t)/'labels.dlabel.nii';nb.save(labels,path)
   info=network_order(image,np.arange(6),path)
   np.testing.assert_array_equal(info['order'],[1,0,2,3,5,4]);self.assertEqual([s['assigned_rows']for s in info['compartments']],[2,2]);np.testing.assert_array_equal(info['colors'][4],[0,1,0,1]);self.assertEqual(info['legend'][0]['subcortical_rows'],1)
   fig,ax=plt.subplots();ax.imshow(np.ones((6,4)))
   draw_network_strip(ax,info);self.assertEqual(len(ax.lines),1);self.assertEqual(ax.lines[0].get_ydata()[0],2.5);plt.close(fig)
   # A cortical-only export must not borrow volume assignments from another map.
   data=np.asarray(labels.dataobj).copy();data[0,3:]=0;nb.save(nb.Cifti2Image(data,labels.header),path);info=network_order(image,np.arange(6),path)
   self.assertEqual(info['compartments'][1]['assigned_rows'],0);np.testing.assert_array_equal(info['colors'][3:],np.tile([.7,.7,.7,1],(3,1)))

if __name__=='__main__':unittest.main()
