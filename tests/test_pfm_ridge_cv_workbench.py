"""A failed PNG capture must be retried, never treated as a cached view."""
import base64,io,subprocess,sys,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from PIL import Image,ImageDraw
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'lib'))
import pfm_ridge_cv_workbench as renderer

class WorkbenchCacheTests(unittest.TestCase):
 def test_png_failure_is_retried_then_cached(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td);resources=root/'resources';priors=resources/'PriorsPFM/Lynch2024/PriorsCortical.mat';priors.parent.mkdir(parents=True);priors.write_bytes(b'priors')
   template=resources/'CrawlingSeedFC/FlatMaps+Inflated.scene';template.parent.mkdir();template.write_text('scene')
   out=root/'cv';(out/'comparisons/p01').mkdir(parents=True)
   for rel in ['FullData_NetworkAssignments.dlabel.nii','comparisons/p01/A_NetworkAssignments.dlabel.nii','comparisons/p01/B_NetworkAssignments.dlabel.nii']:(out/rel).write_bytes(rel.encode())
   surfaces=root/'surfaces';surfaces.mkdir();args=SimpleNamespace(priors_mat=str(priors),surface_dir=str(surfaces),flat_dir=str(surfaces),subject='Example')
   failed=[False];full_calls=[]
   def command(cmd,**kwargs):
    if cmd[1]=='-scene-capture-image':
     scene=Path(cmd[2]);png=Path(cmd[4]);im=Image.new('RGB',(120,72),'white')
     if scene.stem.startswith('FullData'):
      full_calls.append(scene)
      if not failed[0]:
       failed[0]=True;im.save(png);raise subprocess.CalledProcessError(1,cmd)
     ImageDraw.Draw(im).rectangle((10,10,100,60),fill='#285fc3');im.save(png)
    return subprocess.CompletedProcess(cmd,0)
   with patch.object(renderer.features,'check',return_value=False),patch.object(renderer,'build_scene'),patch.object(renderer.subprocess,'run',side_effect=command):
    with self.assertRaises(subprocess.CalledProcessError):renderer.capture_grid(out,args,[dict(index=0,control=False)])
    result=renderer.capture_grid(out,args,[dict(index=0,control=False)])
    self.assertEqual(len(full_calls),2)
    for bank in ['FullData','A','B']:
     data=result['images'][bank]['0'];im=Image.open(io.BytesIO(base64.b64decode(data.split(',',1)[1])))
     self.assertEqual(im.getpixel((50,30)),(40,95,195))
    renderer.capture_grid(out,args,[dict(index=0,control=False)]);self.assertEqual(len(full_calls),2)

if __name__=='__main__':unittest.main()
