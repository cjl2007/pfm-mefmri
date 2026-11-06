#!/usr/bin/env python3
"""Render taxonomy WTA, competing-label stripes, and session confidence to a GIF."""
from pathlib import Path
import argparse,json,os,subprocess,sys
import nibabel as nib
import numpy as np
from PIL import Image,ImageDraw,ImageFont
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from pfm_ridge_scene import build_scene

NETWORKS={'Lynch2024':['Default_Parietal','Salience','Somatomotor_Hand'],
          'Du2024':['DN-A','SAL/PMN','SMOT-A'],
          'Hermosillo2024':['DMN','Sal','SMd'],
          'Kong2019':['DefaultA','VentralAttentionA','SomatomotorA']}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results-dir',required=True);p.add_argument('--subject-dir',required=True)
    p.add_argument('--output',required=True);p.add_argument('--prefix',default='RidgePFM_VTX')
    p.add_argument('--taxonomies',nargs='+',default=list(NETWORKS));p.add_argument('--wta-only',action='store_true')
    p.add_argument('--scene-template',default=str(Path(__file__).resolve().parents[1]/'res0urces/CrawlingSeedFC/FlatMaps+Inflated.scene'))
    a=p.parse_args();results=Path(a.results_dir).resolve();subject=Path(a.subject_dir);output=Path(a.output).resolve();output.parent.mkdir(parents=True,exist_ok=True)
    renders=results/'movie';renders.mkdir(exist_ok=True);env=os.environ.copy();env.update(OMP_NUM_THREADS='4',QT_QPA_PLATFORM='offscreen')
    def font(size,bold=False):
        return ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans'+('-Bold' if bold else '')+'.ttf',size)
    frames=[];durations=[];records=[]
    for tax in a.taxonomies:
        folder=results/tax;wta=folder/(a.prefix+'.dlabel.nii');image=nib.load(wta);table=image.header.get_axis(0).label[0]
        stages=[('wta',wta,0,'Winner-take-all network assignments')]
        if not a.wta_only:
            conf=folder/'Confidence';freq=conf/(a.prefix+'_SessionAssignmentFrequency.dscalar.nii');fimg=nib.load(freq);names=list(fimg.header.get_axis(0).name)
            session_count=len(json.loads((conf/'confidence.json').read_text())['sessions'])
            colored=renders/(tax+'_frequency.dscalar.nii')
            subprocess.run(['wb_command','-cifti-palette',str(freq),'MODE_USER_SCALE',str(colored),'-pos-user','0','1','-palette-name','magma','-disp-pos','true','-disp-zero','true','-disp-neg','false','-interpolate','true'],check=True,env=env)
            stages.append(('stripes',conf/(a.prefix+'_NetworkStripes.dlabel.nii'),0,'Stripes: reproducible FC with competing labels'))
            stages += [('confidence_'+str(i),colored,names.index(n),'Session assignment frequency · '+n.replace('_',' ')) for i,n in enumerate(NETWORKS[tax])]
        for stage,source,index,title in stages:
            stem=renders/(tax+'_'+stage);scene=stem.with_suffix('.scene');png=stem.with_suffix('.png')
            build_scene(a.scene_template,source,subject/'anat/T1w/fsaverage_LR32k',subject/'anat/MNINonLinear/fsaverage_LR32k',subject.name,scene,index, border_paths=[folder/(a.prefix+"."+h+".border") for h in ["L","R"]], compact=True)
            with stem.with_suffix('.log').open('w') as log:
                subprocess.run(['wb_command','-scene-capture-image',str(scene),'1',str(png),'-size-width-height','1200','720'],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
            brain=Image.open(png).convert('RGB')
            if np.asarray(brain).std()<5:raise RuntimeError('Blank scene render: '+str(png))
            frame=Image.new('RGB',(1200,960),'white');frame.paste(brain,(0,100));d=ImageDraw.Draw(frame)
            d.text((28,15),'RidgePFM',font=font(28,True),fill='#142d45');d.text((28,58),title,font=font(22),fill='#334b61')
            box=d.textbbox((0,0),tax,font=font(24,True));x=1170-(box[2]-box[0]);d.text((x,20),tax,font=font(24,True),fill='#142d45')
            if stage.startswith('confidence'):
                import matplotlib
                cmap=matplotlib.colormaps['magma'];strip=(cmap(np.linspace(0,1,400))[:,:3]*255).astype(np.uint8)
                gradient=Image.fromarray(np.repeat(strip[None],18,axis=0));frame.paste(gradient,(400,835))
                d.text((385,860),'0',font=font(17),fill='#334b61');d.text((782,860),'1',font=font(17),fill='#334b61')
                d.text((28,900),f'Fraction of {session_count} independently mapped sessions · empirical consistency, separate from model fit',font=font(16),fill='#526476')
            else:
                labels=[(k,v) for k,v in table.items() if k>0]
                for j,(key,(name,color)) in enumerate(labels):
                    x=25+(j%4)*295;y=823+(j//4)*22;rgb=tuple(int(round(c*255)) for c in color[:3]);d.rectangle((x,y+3,x+12,y+15),fill=rgb)
                    text=name.replace('_',' ').replace('Default ','DMN ').replace('Somatomotor ','SM ').replace('Visual Dorsal/VentralStream','Visual stream')
                    d.text((x+18,y),text,font=font(12),fill='#334b61')
                if stage=='stripes':d.text((28,940),'Only reproducible competing assignments are striped; widths do not encode probabilities.',font=font(12),fill='#526476')
            annotated=stem.with_name(stem.name+'_annotated.png');frame.save(annotated);frames.append(frame);durations.append(3000 if stage in ['wta','stripes'] else 2400)
            records.append(dict(taxonomy=tax,stage=stage,source=str(source),map_index=index,scene=str(scene),image=str(annotated)))
            print(tax,stage,'rendered',flush=True)
    frames[0].save(output,save_all=True,append_images=frames[1:],duration=durations,loop=0,optimize=True,disposal=2)
    (renders/(output.stem+'.json')).write_text(json.dumps(dict(subject=subject.name,frames=records,durations_ms=durations,
        note='Taxonomy-specific labels; the selected networks are illustrative and are not asserted to be exact equivalents. Session frequencies are not calibrated biological membership probabilities.'),indent=2)+'\n')
    print(output,flush=True)
if __name__=='__main__':main()
