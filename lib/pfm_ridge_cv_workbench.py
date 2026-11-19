"""Cached Workbench captures for the CV report, using the README scene."""
from pathlib import Path
import base64, hashlib, io, json, os, subprocess
from concurrent.futures import ThreadPoolExecutor
from PIL import Image, features
from pfm_ridge_scene import build_scene

def capture_grid(out,args,rows):
    out=Path(out)
    template=Path(args.priors_mat).resolve().parents[2]/'CrawlingSeedFC/FlatMaps+Inflated.scene'
    if not template.is_file():raise FileNotFoundError(template)
    sources={'FullData':out/'FullData_NetworkAssignments.dlabel.nii',
             'A':out/'comparisons/p01/A_NetworkAssignments.dlabel.nii',
             'B':out/'comparisons/p01/B_NetworkAssignments.dlabel.nii'}
    env=os.environ.copy();env.update(QT_QPA_PLATFORM='offscreen',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    identity=dict(template=hashlib.sha256(template.read_bytes()).hexdigest(),
                  capture_code=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scene_code=hashlib.sha256(Path(__import__('pfm_ridge_scene').__file__).read_bytes()).hexdigest(),
                  subject=args.subject,width=1200,height=720,
                  sources={k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in sources.items()})
    identity['surfaces']=[(str(p.resolve()),p.stat().st_size,p.stat().st_mtime_ns)
                          for folder in [args.surface_dir,args.flat_dir]
                          for p in sorted(Path(folder).glob(args.subject+'.*.surf.gii'))]
    webp=features.check('webp');identity['image_format']='webp' if webp else 'png'
    key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:20]
    cache=out/'workbench_views'/key;cache.mkdir(parents=True,exist_ok=True)
    tasks=[(bank,row['index']) for bank in sources for row in rows
           if not row['control'] or bank=='FullData']
    def capture(task):
        bank,index=task;stem=cache/f'{bank}_{index:02d}';dest=stem.with_suffix('.webp' if webp else '.view.png')
        if not dest.exists():
            borders=[stem.with_suffix('.'+h+'.border') for h in ['L','R']]
            cmd=['wb_command','-cifti-label-to-border',str(sources[bank]),'-column',str(index+1)]
            for h,border in zip(['L','R'],borders):
                cmd+=['-border',str(Path(args.surface_dir)/f'{args.subject}.{h}.midthickness.32k_fs_LR.surf.gii'),str(border)]
            scene=stem.with_suffix('.scene');png=stem.with_suffix('.png')
            with stem.with_suffix('.log').open('w') as log:
                subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
                build_scene(template,sources[bank],args.surface_dir,args.flat_dir,args.subject,scene,index,border_paths=borders,compact=True)
                subprocess.run(['wb_command','-scene-capture-image',str(scene),'1',str(png),'-size-width-height','1200','720'],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
            im=Image.open(png).convert('RGB')
            if max(im.getextrema()[c][1]-im.getextrema()[c][0] for c in range(3))<20:raise ValueError('Blank Workbench capture: '+str(png))
            temp=dest.with_suffix('.tmp'+dest.suffix)
            if webp:im.save(temp,format='WEBP',lossless=True)
            else:im.save(temp,format='PNG')
            temp.replace(dest)
        return bank,index,('data:image/webp;base64,' if webp else 'data:image/png;base64,')+base64.b64encode(dest.read_bytes()).decode()
    images={bank:{} for bank in sources}
    with ThreadPoolExecutor(max_workers=3) as pool:
        for n,(bank,index,data) in enumerate(pool.map(capture,tasks),1):
            images[bank][str(index)]=data
            if n%12==0 or n==len(tasks):print(f'[ridge-cv] Workbench views {n}/{len(tasks)}',flush=True)
    (cache/'manifest.json').write_text(json.dumps(identity,indent=2)+'\n')
    return dict(renderer='Connectome Workbench',width=1200,height=720,images=images)
