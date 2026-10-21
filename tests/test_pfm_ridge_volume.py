"""Real Workbench ribbon export on a small synthetic cortical ribbon."""
import sys,unittest,tempfile,shutil,subprocess
from pathlib import Path
import nibabel as nib
import numpy as np
from nibabel.cifti2 import BrainModelAxis,ScalarAxis,LabelAxis
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from pfm_ridge_volume import export

@unittest.skipUnless(shutil.which('wb_command'),'Workbench required')
class RibbonExportTests(unittest.TestCase):
    def test_probability_scale_integer_identity_and_coverage(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder);surfaces={};vertices=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]],np.float32)
            faces=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]],np.int32)
            for h,cx in [('L',5),('R',14)]:
                surfaces[h]=[]
                for kind,radius in [('mid',2),('white',1.3),('pial',2.7)]:
                    coords=vertices*radius+np.array([cx,5,5],np.float32)
                    img=nib.gifti.GiftiImage(darrays=[nib.gifti.GiftiDataArray(coords,intent='NIFTI_INTENT_POINTSET'),nib.gifti.GiftiDataArray(faces,intent='NIFTI_INTENT_TRIANGLE')]);path=p/(h+kind+'.surf.gii');nib.save(img,path)
                    subprocess.run(['wb_command','-set-structure',str(path),'CORTEX_LEFT' if h=='L' else 'CORTEX_RIGHT'],check=True);surfaces[h].append(path)
            ref=p/'reference.nii.gz';nib.save(nib.Nifti1Image(np.zeros((20,11,11),np.uint8),np.eye(4)),ref)
            mask=np.zeros((20,11,11),bool);mask[10,5,5]=True
            axis=BrainModelAxis.from_surface(np.arange(4),4,'CortexLeft')+BrainModelAxis.from_surface(np.arange(4),4,'CortexRight')+BrainModelAxis.from_mask(mask,affine=np.eye(4),name='ThalamusLeft')
            scalar=p/'frequency.dscalar.nii';nib.save(nib.Cifti2Image(np.full((1,9),.6,np.float32),nib.Cifti2Header.from_axes((ScalarAxis(['frequency']),axis))),scalar)
            dest=p/'frequency.nii.gz';export(scalar,ref,surfaces,dest)
            values=nib.load(dest).get_fdata()[...,0];coverage=nib.load(p/'frequency_Coverage.nii.gz').get_fdata()
            self.assertTrue((coverage==1).any());np.testing.assert_allclose(values[coverage>0],.6,atol=1e-6);self.assertTrue(np.all(values[coverage==0]==0));self.assertEqual(coverage[10,5,5],2)
            # A nondefault zero label makes Workbench renumber its label table.
            table={0:('Unavailable',(0,0,0,0)),1:('A',(1,0,0,1)),2:('B',(0,0,1,1))}
            label=p/'wta.dlabel.nii';nib.save(nib.Cifti2Image(np.array([[1]*4+[2]*4+[2]],np.int32),nib.Cifti2Header.from_axes((LabelAxis(['WTA'],[table]),axis))),label)
            export(label,ref,surfaces,p/'wta.nii.gz');values=nib.load(p/'wta.nii.gz').get_fdata()[...,0]
            self.assertEqual(set(np.unique(values)),{0.,1.,2.});self.assertTrue((values[:10]==1).any());self.assertFalse((values[:10]==2).any());self.assertEqual(values[10,5,5],2)
            exported=p/'labels.txt'
            subprocess.run(['wb_command','-volume-label-export-table',str(p/'wta.nii.gz'),'1',str(exported)],check=True)
            self.assertIn('A',exported.read_text());self.assertIn('B',exported.read_text())
if __name__=='__main__':unittest.main()
