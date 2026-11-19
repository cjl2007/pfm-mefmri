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
from mefmri_figure_style import configure_matplotlib,figure_font_path,report_css
configure_matplotlib()
from PIL import Image

def rms_matrix(values):
    x=np.asarray(values,dtype=np.float64)
    if x.ndim!=2 or not x.shape[1] or not np.isfinite(x).all():raise ValueError('Expected finite scan-by-voxel values')
    d=np.mean(x*x,axis=1);sq=d[:,None]+d[None,:]-2*(x@x.T)/x.shape[1]
    out=100*np.sqrt(np.maximum(sq,0));out[np.abs(sq)<1e-12]=0;np.fill_diagonal(out,0)
    return out

def render_alignment_frames(subject_dir, reference, normalized, items, mask, center, limits):
    """Sample fixed anatomical planes at display resolution, with fixed tissue contours.

    Sampling affects report pixels only. Quantitative QC continues to use the
    original functional reference grid and undilated brain mask.
    """
    from scipy.ndimage import map_coordinates, binary_dilation
    from PIL import ImageDraw, ImageFont
    from matplotlib import font_manager

    anatomy_dir = Path(subject_dir) / 'anat/T1w'
    brain_path = anatomy_dir / 'T1w_acpc_dc_restore_brain.nii.gz'
    anatomy = nib.as_closest_canonical(nib.load(brain_path))
    anatomy_data = np.asarray(anatomy.dataobj, dtype=np.float32)
    anatomical_mask = binary_fill_holes(anatomy_data > 0)
    spacing = nib.affines.voxel_sizes(anatomy.affine)
    anatomical_display = distance_transform_edt(~anatomical_mask, sampling=spacing) <= 3
    coords = np.array(np.where(anatomical_mask))
    bounds = [(max(0, int(coords[d].min()) - int(np.ceil(3 / spacing[d]))),
               min(anatomy.shape[d], int(coords[d].max()) + int(np.ceil(3 / spacing[d])) + 1))
              for d in range(3)]
    center_world = nib.affines.apply_affine(reference.affine, center)
    anatomical_center = nib.affines.apply_affine(np.linalg.inv(anatomy.affine), center_world)
    tissue_specs = []
    sources = {}
    ribbon_path = anatomy_dir / 'ribbon.nii.gz'
    aseg_path = anatomy_dir / 'aparc+aseg.nii.gz'
    if ribbon_path.exists():
        ribbon = nib.load(ribbon_path)
        labels = np.asarray(ribbon.dataobj)
        tissue_specs.extend([
            ('White matter', np.isin(labels, [2, 41]), ribbon.affine, (67, 221, 232)),
            ('Outer cortical gray matter', np.isin(labels, [2, 3, 41, 42]), ribbon.affine, (255, 173, 69)),
        ])
        sources['ribbon'] = str(ribbon_path)
    if aseg_path.exists():
        aseg = nib.load(aseg_path)
        labels = np.asarray(aseg.dataobj)
        tissue_specs.append(('Ventricular CSF', np.isin(labels, [4, 5, 14, 15, 43, 44, 72]),
                             aseg.affine, (193, 135, 255)))
        sources['ventricular_segmentation'] = str(aseg_path)
    if not tissue_specs:
        tissue_specs.append(('Brain boundary', anatomical_mask, anatomy.affine, (255, 173, 69)))
    tlim = np.percentile(anatomy_data[anatomical_mask], [1, 99.5])
    plane_geometry = []
    tile_width, tile_height = 600, 440
    for axis in range(3):
        other = [d for d in range(3) if d != axis]
        width_mm = (bounds[other[0]][1] - bounds[other[0]][0]) * spacing[other[0]]
        height_mm = (bounds[other[1]][1] - bounds[other[1]][0]) * spacing[other[1]]
        # One world-space location per output pixel, with isotropic physical scale.
        scale = min((tile_width - 24) / width_mm, (tile_height - 24) / height_mm)
        width = max(2, int(round(width_mm * scale)))
        height = max(2, int(round(height_mm * scale)))
        horizontal = np.linspace(bounds[other[0]][0], bounds[other[0]][1] - 1, width)
        vertical = np.linspace(bounds[other[1]][1] - 1, bounds[other[1]][0], height)
        yy, xx = np.meshgrid(vertical, horizontal, indexing='ij')
        ijk = np.empty((3, height * width))
        ijk[axis] = anatomical_center[axis]
        ijk[other[0]] = xx.ravel()
        ijk[other[1]] = yy.ravel()
        world = nib.affines.apply_affine(anatomy.affine, ijk.T)
        functional_ijk = nib.affines.apply_affine(np.linalg.inv(reference.affine), world).T
        keep = map_coordinates(anatomical_display.astype(np.float32), ijk, order=1,
                               mode='constant', cval=0).reshape(height, width) >= .5
        edges = []
        for name, tissue, affine, color in tissue_specs:
            tissue_ijk = nib.affines.apply_affine(np.linalg.inv(affine), world).T
            section = map_coordinates(tissue.astype(np.float32), tissue_ijk, order=1,
                                      mode='constant', cval=0).reshape(height, width) >= .5
            edge = binary_dilation(section & ~binary_erosion(section))
            edges.append((edge & keep, color))
        t1 = map_coordinates(anatomy_data, ijk, order=1, mode='constant', cval=0).reshape(height, width)
        plane_geometry.append((functional_ijk, keep, edges, t1, (width, height)))

    def panel(values, keep, edges, display_range):
        lo, hi = display_range
        gray = np.uint8(np.clip((values - lo) / max(float(hi - lo), 1e-8), 0, 1) * 255)
        rgb = np.repeat(gray[:, :, None], 3, axis=2)
        rgb[~keep] = 0
        for edge, color in edges:
            rgb[edge] = color
        image = Image.fromarray(rgb)
        tile = Image.new('RGB', (tile_width, tile_height), 'black')
        tile.paste(image, ((tile_width - image.width) // 2, (tile_height - image.height) // 2))
        return tile

    font_path = figure_font_path()
    font = ImageFont.truetype(font_path, 24)
    small_font = ImageFont.truetype(font_path, 21)
    anatomy_tiles = [panel(t1, keep, edges, tlim) for _, keep, edges, t1, _ in plane_geometry]
    frames = []
    for i, ((run, path, scale), data) in enumerate(zip(items, normalized)):
        canvas = Image.new('RGB', (1800, 1060), 'black')
        draw = ImageDraw.Draw(canvas)
        draw.text((20, 12), f'Scan {i + 1} · session {run.session}, run {run.run} · functional reference',
                  fill='white', font=font)
        for j, (name, tissue, affine, color) in enumerate(tissue_specs):
            draw.text((20 + j * 570, 48), name, fill=color, font=small_font)
        for axis, (ijk, keep, edges, t1, size) in enumerate(plane_geometry):
            values = map_coordinates(data, ijk, order=1, mode='constant', cval=0).reshape(size[1], size[0])
            canvas.paste(panel(values, keep, edges, limits), (axis * tile_width, 105))
            canvas.paste(anatomy_tiles[axis], (axis * tile_width, 580))
            draw.text((axis * tile_width + 20, 82), ['Sagittal', 'Coronal', 'Axial'][axis],
                      fill='#cbd6e5', font=small_font)
        draw.text((20, 550), 'Fixed T1 anatomy · identical slice coordinates and tissue contours', fill='white', font=font)
        draw.text((20, 1026), 'Display sampling only · fixed anatomical brain mask + 3 mm margin · source data unchanged',
                  fill='#cbd6e5', font=small_font)
        frames.append(canvas)
    return frames, dict(anatomy_source=str(brain_path), anatomy_voxel_size_mm=spacing.tolist(),
                        tissue_sources=sources, tissue_legend=[name for name, *_ in tissue_specs],
                        canvas_pixels=[1800, 1060], display_interpolation='linear',
                        fixed_world_slice_center_mm=center_world.tolist(),
                        quantitative_reference_grid_unchanged=True)


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
    frames,render_info=render_alignment_frames(subject_dir,ref,normalized,items,mask,center,limits)
    gif=reportlets/'subj_desc-coreg_unwarpedSBRefBBR_fade.gif'
    sampled=frames[::max(1,len(frames)//12)][:12];sheet=Image.new('RGB',(1600,660),'black')
    for j,im in enumerate(sampled):sheet.paste(im.resize((400,220)),((j%4)*400,(j//4)*220))
    palette=Image.new('P',(1,1));palette.putpalette([c for v in np.linspace(0,255,253).astype(int) for c in (int(v),)*3]+[67,221,232,255,173,69,193,135,255])
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
    (reportlets/'subj_desc-coreg_RMSMethods.json').write_text(json.dumps(dict(definition='100 * sqrt(mean((I_i/median_brain(I_i) - I_j/median_brain(I_j))**2))',mask_source=str(brainpath),boundary_definition='One reference-grid voxel inward from fixed brain-mask edge',display_margin_mm=3,display_renderer=render_info,intensity_limits_shared_across_scans=True,first_scan_adjacent_value=None,summary=summary,interpretation='Image consistency, affected by contrast, noise, and registration; not displacement in millimeters or proof of anatomical correctness'),indent=2)+'\n')
    return dict(display_limits=limits.tolist(),display_renderer=render_info,anchor_png='subj_desc-coreg_T1wAnchor.png',fade_gif=gif.name,frame_pngs=[],similarity_matrix_png=pairplot.name,stability_png=plot.name,metrics=rows,run_rms=distance.tolist(),matrix_kind='RMS difference; lower is more similar',session_summary=summary)
