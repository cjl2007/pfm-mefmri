"""Brain-focused scan-reference movies and intensity-normalized RMS stability QC."""
from pathlib import Path
import json,os
os.environ.setdefault("MPLCONFIGDIR","/tmp/mefmri_stability_mpl")
import numpy as np
import nibabel as nib
from scipy.ndimage import binary_fill_holes,binary_erosion,distance_transform_edt
from nibabel.processing import resample_from_to
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

def rms_matrix(values):
    x=np.asarray(values,dtype=np.float64)
    if x.ndim!=2 or not x.shape[1] or not np.isfinite(x).all():raise ValueError('Expected finite scan-by-voxel values')
    d=np.mean(x*x,axis=1);sq=d[:,None]+d[None,:]-2*(x@x.T)/x.shape[1]
    out=100*np.sqrt(np.maximum(sq,0));out[np.abs(sq)<1e-12]=0;np.fill_diagonal(out,0)
    return out

def generate(subject_dir,func_name,xfms_name,runs,reportlets,resolve_reference,display_limits=None):
    reportlets=Path(reportlets);reportlets.mkdir(parents=True,exist_ok=True)
    t1path=Path(subject_dir)/'func/xfms'/xfms_name/'T1w_acpc_func.nii.gz'
    if not t1path.exists():t1path=t1path.with_name('T1w_acpc_brain_func.nii.gz')
    if not t1path.exists():return {}
    ref=nib.load(t1path);t1=np.asarray(ref.dataobj,dtype=np.float32)
    if t1.ndim!=3:raise ValueError('Anatomical reference must be 3D')
    brainpath=Path(subject_dir)/'anat/T1w/T1w_acpc_dc_restore_brain.nii.gz'
    if not brainpath.exists():raise ValueError('Missing brain-extracted T1w: '+str(brainpath))
    high=nib.load(brainpath);maskimg=nib.Nifti1Image((np.asarray(high.dataobj)>0).astype('uint8'),high.affine)
    mask=np.asarray(resample_from_to(maskimg,(ref.shape,ref.affine),order=0).dataobj)>0;mask=binary_fill_holes(mask)
    if mask.sum()<20:raise ValueError('Insufficient anatomical brain coverage')
    # A fixed physical margin is for display only; quantitative masks are undilated.
    zooms=nib.affines.voxel_sizes(ref.affine);display=distance_transform_edt(~mask,sampling=zooms)<=3
    edge=mask&~binary_erosion(mask);items=[];normalized=[]
    for run in sorted(runs,key=lambda r:(r.session,r.run)):
        path=resolve_reference(run,t1path,reportlets/'nifti')
        if path is None:continue
        im=nib.load(path)
        if im.shape[:3]!=ref.shape or not np.allclose(im.affine,ref.affine,atol=1e-4):raise ValueError('Mismatched reference geometry: '+str(path))
        data=np.asarray(im.dataobj,dtype=np.float32)
        if data.ndim==4 and data.shape[3]==1:data=data[...,0]
        if data.ndim!=3:raise ValueError('Reference must be a single 3D image: '+str(path))
        if not np.isfinite(data[mask]).all():raise ValueError('Nonfinite brain data: '+str(path))
        positive=data[mask&(data>0)]
        if not len(positive):raise ValueError('No positive brain intensities: '+str(path))
        scale=float(np.median(positive));normalized.append(data/scale);items.append((run,Path(path),scale))
    if not items:return {}
    x=np.stack([v[mask] for v in normalized]);b=np.stack([v[edge] for v in normalized]);distance=rms_matrix(x);boundary=rms_matrix(b);n=len(items)
    rows=[]
    for i,(run,path,scale) in enumerate(items):
        others=np.arange(n)!=i
        rows.append(dict(session=run.session,run=run.run,scan=i+1,source_reference=str(path),median_brain_intensity=scale,previous_scan=i if i else None,between_sessions=(run.session!=items[i-1][0].session) if i else None,adjacent_rms_percent=float(distance[i,i-1]) if i else None,adjacent_boundary_rms_percent=float(boundary[i,i-1]) if i else None,median_rms_to_other_scans_percent=float(np.median(distance[i,others])) if n>1 else None))
    scan=np.arange(1,n+1);cross=np.array([bool(r['between_sessions']) for r in rows]);fig,axs=plt.subplots(4,1,figsize=(11,12),constrained_layout=True)
    for key,label,color in [('adjacent_rms_percent','Whole brain','#3568aa'),('adjacent_boundary_rms_percent','Brain boundary','#b77925')]:axs[0].plot(scan,[r[key] for r in rows],'o-',ms=3,color=color,label=label)
    for p in scan[cross]:axs[0].axvline(p,color='gray',alpha=.3)
    axs[0].set(title='Neighboring scan-reference differences',ylabel='RMS difference (%)',xlabel='Scan number · vertical lines mark session transitions');axs[0].legend(frameon=False)
    axs[1].plot(scan,[r['median_rms_to_other_scans_percent'] for r in rows],'o-',ms=3,color='#168b80');axs[1].set(title='Each scan versus all other scans · median pairwise RMS',ylabel='Median RMS (%)',xlabel='Scan number')
    img=axs[2].imshow(distance,cmap='magma',origin='upper',interpolation='nearest');axs[2].set(title='Pairwise scan-reference RMS differences',xlabel='Scan index (zero-based)',ylabel='Scan index (zero-based)');fig.colorbar(img,ax=axs[2],label='RMS difference (%)')
    for j,(key,color) in enumerate([('adjacent_rms_percent','#3568aa'),('adjacent_boundary_rms_percent','#b77925')]):
        for k,between in enumerate([False,True]):
            vals=[r[key] for r in rows[1:] if r['between_sessions']==between];pos=j*3+k
            if vals:
                axs[3].scatter(np.linspace(pos-.12,pos+.12,len(vals)),vals,color=color,alpha=.65,s=18);axs[3].plot([pos-.2,pos+.2],[np.median(vals)]*2,color='black',lw=2)
    axs[3].set(xticks=[0,1,3,4],xticklabels=['Within session','Between sessions','Within session','Between sessions'],ylabel='Adjacent RMS (%)',title='Whole brain (blue) and boundary (ochre) · medians and individual transitions')
    plot=reportlets/'subj_desc-coreg_RMSStability.png';fig.savefig(plot,dpi=140);plt.close(fig)
    pairplot=reportlets/'subj_desc-coreg_RMSMatrix.png';fig,ax=plt.subplots(figsize=(7,6));m=ax.imshow(distance,cmap='magma');fig.colorbar(m,ax=ax,label='RMS difference (%)');ax.set(title='Pairwise RMS · normalized reference images',xlabel='Scan index (zero-based)',ylabel='Scan index (zero-based)');fig.tight_layout();fig.savefig(pairplot,dpi=140);plt.close(fig)
    np.save(reportlets/'subj_desc-coreg_pairwiseRMS.npy',distance)
    # Same slices, mask, intensity limits and contours for every scan; no fade/interpolated frames.
    coords=np.array(np.where(mask));center=np.median(coords,axis=1).astype(int);center[0]=min(ref.shape[0]-1,center[0]+int(round(8/zooms[0])));bounds=[(max(0,int(coords[d].min())-2),min(ref.shape[d],int(coords[d].max())+3)) for d in range(3)]
    limits=np.asarray(display_limits) if display_limits is not None else np.percentile(np.concatenate([v[mask][::8] for v in normalized]),[1,99.5]);tlim=np.percentile(t1[mask],[1,99.5]);frames=[]
    def pane(data,axis,lim):
        sl=np.rot90(np.take(data,center[axis],axis=axis));keep=np.rot90(np.take(display,center[axis],axis=axis));cont=np.rot90(np.take(mask,center[axis],axis=axis));other=[d for d in range(3) if d!=axis];a,c=bounds[other[1]];d,e=bounds[other[0]];# rotated array reverses its first dimension
        roi=(slice(ref.shape[other[1]]-c,ref.shape[other[1]]-a),slice(d,e));sl=sl[roi];keep=keep[roi];cont=cont[roi];g=np.uint8(np.clip((sl-lim[0])/(lim[1]-lim[0]),0,1)*255);rgb=np.repeat(g[:,:,None],3,axis=2);rgb[~keep]=0;rgb[cont&~binary_erosion(cont)]=(255,173,69);im=Image.fromarray(rgb);im.thumbnail((390,260));tile=Image.new('RGB',(400,270),'black');tile.paste(im,((400-im.width)//2,(270-im.height)//2));return tile
    from PIL import ImageDraw,ImageFont
    from matplotlib import font_manager
    font=ImageFont.truetype(font_manager.findfont("DejaVu Sans"),18)
    for i,(run,path,scale) in enumerate(items):
        canvas=Image.new('RGB',(1200,620),'black');draw=ImageDraw.Draw(canvas);draw.text((15,8),f'Scan {i+1} · session {run.session}, run {run.run} · functional reference',fill='white',font=font)
        for axis in range(3):canvas.paste(pane(normalized[i],axis,limits),(axis*400,32));canvas.paste(pane(t1,axis,tlim),(axis*400,330))
        draw.text((15,310),'Fixed T1 anatomy · fixed brain mask + 3 mm display margin',fill='white',font=font);frames.append(canvas)
    gif=reportlets/'subj_desc-coreg_unwarpedSBRefBBR_fade.gif'
    sampled=frames[::max(1,len(frames)//12)][:12];sheet=Image.new('RGB',(1600,660),'black')
    for j,im in enumerate(sampled):sheet.paste(im.resize((400,220)),((j%4)*400,(j//4)*220))
    palette=sheet.quantize(colors=256,method=Image.Quantize.MEDIANCUT,dither=Image.Dither.NONE)
    encoded=[im.quantize(palette=palette,dither=Image.Dither.NONE) for im in frames]
    encoded[0].save(gif,save_all=True,append_images=encoded[1:],duration=250,loop=0,disposal=2,optimize=False)
    frames[0].save(reportlets/'subj_desc-coreg_T1wAnchor.png')
    import csv
    transitions=reportlets/'subj_desc-coreg_RMSTransitions.tsv'
    with transitions.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
    summary={}
    for val,label in [(False,'within_session'),(True,'between_sessions')]:
        values=[r['adjacent_rms_percent'] for r in rows[1:] if r['between_sessions']==val];summary[label]=dict(transitions=len(values),median_rms_percent=float(np.median(values)) if values else None)
    (reportlets/'subj_desc-coreg_RMSMethods.json').write_text(json.dumps(dict(definition='100 * sqrt(mean((I_i/median_brain(I_i) - I_j/median_brain(I_j))**2))',mask_source=str(brainpath),boundary_definition='One reference-grid voxel inward from fixed brain-mask edge',display_margin_mm=3,intensity_limits_shared_across_scans=True,first_scan_adjacent_value=None,summary=summary,interpretation='Image consistency, affected by contrast, noise, and registration; not displacement in millimeters or proof of anatomical correctness'),indent=2)+'\n')
    return dict(display_limits=limits.tolist(),anchor_png='subj_desc-coreg_T1wAnchor.png',fade_gif=gif.name,frame_pngs=[],similarity_matrix_png=pairplot.name,stability_png=plot.name,metrics=rows,run_rms=distance.tolist(),matrix_kind='RMS difference; lower is more similar',session_summary=summary)
