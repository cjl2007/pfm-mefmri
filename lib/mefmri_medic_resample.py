"""Staged production MEDIC joint resampling interface (not yet installed).

Final intensities come exclusively from native, non-spatially-resampled echoes.
Motion matrices map unwarped native coordinates to the registration intermediate.
The base warp maps that intermediate to the requested reference (T1w or atlas).
"""
import argparse,hashlib,json,os,re,subprocess,tempfile,shutil
from pathlib import Path
from functools import lru_cache
import numpy as np
import nibabel as nib
from scipy.ndimage import spline_filter
from medic_geometry import AXES,scaled_voxel_affine,validate_grid,frame_indices
from joint_sampler import fsl_temporal_weights,sample_joint_frame

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def oriented_image(path,reorient):
    im=nib.load(path)
    q,qcode=im.get_qform(coded=True);s,scode=im.get_sform(coded=True)
    if qcode and scode and not np.allclose(q,s,atol=1e-3,rtol=0):
        raise ValueError(f'Conflicting qform/sform: {path}')
    if reorient:
        transform=nib.orientations.ornt_transform(nib.orientations.io_orientation(im.affine),nib.orientations.axcodes2ornt(('L','A','S')))
        # The timing file is indexed by the original slices. Don't silently
        # permute/reverse that axis without remapping the timing metadata.
        if not np.array_equal(transform[2],[2,1]):raise ValueError('Reorientation would change slice axis/order')
        im=im.as_reoriented(transform)
    if nib.aff2axcodes(im.affine) not in [('L','A','S'),('R','A','S')]:
        raise ValueError('Unsupported native orientation; explicit slice/field remapping required')
    return im

def preflight(args,require_motion=True):
    if args.axis not in AXES:raise ValueError('Invalid phase-encoding axis')
    if not np.isfinite(args.tr) or args.tr<=0:raise ValueError('Invalid TR')
    if not re.fullmatch(r'[A-Za-z0-9_.-]+',args.prefix):raise ValueError('Unsafe output prefix')
    inputs=[Path(p).resolve(strict=True) for p in args.inputs]
    if len(set(inputs))!=len(inputs):raise ValueError('Duplicate echo inputs')
    images=[oriented_image(p,args.reorient) for p in inputs];native=images[0]
    if native.ndim!=4:raise ValueError('Echo data must be 4D')
    for p,im in zip(inputs,images):
        if im.shape!=native.shape:raise ValueError(f'Echo dimensions differ: {p}')
        validate_grid(im,native)
        if abs(im.header.get_zooms()[3]-args.tr)>.001:raise ValueError(f'Echo TR disagrees: {p}')
    maps=nib.load(args.dmap);validate_grid(maps,native)
    if maps.ndim not in (3,4):raise ValueError('MEDIC map must be 3D or 4D')
    pairs=frame_indices(native.shape[3],1 if maps.ndim==3 else maps.shape[3],args.discard)
    offsets=np.loadtxt(args.timing,ndmin=1)
    if offsets.shape!=(native.shape[2],) or not np.isfinite(offsets).all() or np.any(abs(offsets)>=1):
        raise ValueError('Timing must provide one finite fractional-TR shift per slice')
    ref=nib.load(args.reference);mask=nib.load(args.mask)
    if ref.ndim!=3 or mask.shape!=ref.shape:raise ValueError('Reference/mask must be matching 3D images')
    validate_grid(mask,ref)
    if not np.any(mask.get_fdata()>0):raise ValueError('Empty output mask')
    Path(args.base_warp).resolve(strict=True)
    matrices=[]
    if require_motion:
        folder=Path(args.motion_dir)
        names=sorted(p for p in folder.iterdir() if re.fullmatch(r'MAT_[0-9]{4,}',p.name))
        expected=[folder/f'MAT_{t:04d}' for t in range(len(pairs))]
        if names!=expected:raise ValueError('Motion matrices must exactly index retained frames from zero')
        for p in names:
            m=np.loadtxt(p)
            if m.shape!=(4,4) or not np.isfinite(m).all() or not np.allclose(m[3],[0,0,0,1],atol=1e-6) or abs(np.linalg.det(m[:3,:3]))<1e-6:
                raise ValueError(f'Invalid motion matrix {p}')
            matrices.append(m)
    return inputs,images,maps,pairs,offsets,ref,mask,matrices

def run(args):
    checked=preflight(args,require_motion=not args.preflight_only)
    inputs,images,mapimage,pairs,offsets,ref,maskimage,matrices=checked
    if args.preflight_only:
        print(json.dumps(dict(preflight=True,echoes=len(images),retained_volumes=len(pairs))),flush=True);return
    out=Path(args.output_dir).resolve();out.mkdir(parents=True,exist_ok=True)
    lock=out/f'.{args.prefix}_medic.lock'
    try:fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600)
    except FileExistsError:raise RuntimeError(f'Concurrent run or stale lock; inspect before retrying: {lock}')
    os.close(fd)
    work=Path(tempfile.mkdtemp(prefix=f'.{args.prefix}_medic_',dir=out))
    try:
        source_paths=[*inputs,Path(args.dmap),Path(args.timing),Path(args.base_warp),Path(args.reference),Path(args.mask),*sorted(Path(args.motion_dir).glob('MAT_[0-9]*'))]
        hashes={str(p.resolve()):digest(p) for p in source_paths}
        # Keep the actual motion matrices after the shell module cleans MCF.
        # Transform provenance must not depend on optional scout retention.
        archive=work/'motion';archive.mkdir()
        for t,m in enumerate(matrices):np.savetxt(archive/f'MAT_{t:04d}',m,fmt='%.16g')
        shutil.copy2(args.timing,work/'SliceTiming.txt')
        base=work/'base_relative.nii.gz';env=os.environ.copy();env['FSLOUTPUTTYPE']='NIFTI_GZ'
        subprocess.run(['convertwarp',f'--ref={args.reference}',f'--warp1={args.base_warp}',f'--out={base}','--relout'],check=True,capture_output=True,env=env)
        warp=nib.load(base);validate_grid(warp,ref)
        if warp.shape!=(*ref.shape,3):raise ValueError('Base warp must be a dense three-component field')
        field=warp.get_fdata(dtype=np.float32)
        if not np.isfinite(field).all():raise ValueError('Nonfinite anatomical field')
        mask=maskimage.get_fdata()>0;indices=np.argwhere(mask)
        xyz=np.c_[indices,np.ones(len(indices))].T
        target=(scaled_voxel_affine(ref)@xyz)[:3]
        intermediate=target+field[tuple(indices.T)].T
        inverse_scaled=np.linalg.inv(scaled_voxel_affine(images[0]))
        pulls=[(inverse_scaled@np.linalg.inv(m)@np.vstack([intermediate,np.ones(len(indices))]))[:3] for m in matrices]
        maps=mapimage.get_fdata(dtype=np.float32)
        if not np.isfinite(maps).all():raise ValueError('Nonfinite MEDIC field')
        if maps.ndim==4 and maps.shape[3]>1:maps=maps[...,args.discard:]
        weights=fsl_temporal_weights(len(pairs),offsets,work/'temporal_weights')
        outputs=[]
        for e,im in enumerate(images,1):
            raw=im.get_fdata(dtype=np.float32)[...,args.discard:]
            if not np.isfinite(raw).all():raise ValueError(f'Nonfinite native echo {e}')
            @lru_cache(maxsize=raw.shape[2]*10)
            def coefficients(k,s):return spline_filter(raw[:,:,s,k],order=3,output=np.float32,mode='mirror')
            array=np.zeros((*ref.shape,len(pairs)),np.float32)
            for t,q in enumerate(pulls):
                values=sample_joint_frame(raw,maps,im.affine,args.axis,lambda k:pulls[k],weights,t,coefficients)
                if not np.isfinite(values).all():raise ValueError(f'Nonfinite resampled echo/frame: {e}/{t}')
                array[tuple(indices.T)+(np.full(len(indices),t),)]=values
                if t%25==0:print(f'MEDIC single-spatial-resampling echo={e} volume={t}/{len(pairs)}',flush=True)
            filename=f'{args.prefix}_E{e}_acpc.nii.gz';path=work/filename
            header=ref.header.copy();header.set_data_dtype(np.float32)
            result=nib.Nifti1Image(array,ref.affine,header);result.header.set_zooms((*ref.header.get_zooms()[:3],args.tr));nib.save(result,path)
            outputs.append(dict(filename=filename,sha256=digest(path),echo=e))
            coefficients.cache_clear()
        for path,h in hashes.items():
            if digest(path)!=h:raise ValueError(f'Input changed during resampling: {path}')
        code_hashes={p.name:digest(p) for p in [Path(__file__),Path(__file__).with_name('medic_geometry.py'),Path(__file__).with_name('joint_sampler.py')]}
        manifest=dict(schema=1,complete=True,method='joint_native_slice_temporal_spatial_v2',input_sha256=hashes,code_sha256=code_hashes,frame_pairs=pairs,outputs=outputs,tr=args.tr,axis=args.axis,base_relative_sha256=digest(base),motion_archive=str(archive),spatial_kernel='in-plane cubic B-spline; through-plane Keys cubic a=-0.5',temporal_kernel='installed FSL slicetimer impulse response, per contributing native slice',motion_semantics='Each contributing acquisition frame uses its own motion matrix and MEDIC field',work_dir=str(work),postprocessing='Output hashes describe pre-N4 intensity data; shell postprocessing may change intensities, never geometry.')
        # Publish only after every echo succeeds. Preserve replaced outputs and
        # record a completed manifest last; consumers must not trust partial runs.
        backup=work/'previous_outputs';backup.mkdir()
        published=[];saved=[]
        try:
            for row in outputs:
                destination=out/row['filename']
                if destination.exists():destination.replace(backup/destination.name);saved.append(destination.name)
                (work/row['filename']).replace(destination);published.append(destination.name)
            (work/'manifest.json').write_text(json.dumps(manifest,indent=2))
            (work/'manifest.json').replace(out/f'{args.prefix}_MEDIC_single_resampling.json')
        except BaseException:
            for name in published:(out/name).replace(work/name)
            for name in saved:(backup/name).replace(out/name)
            raise
        print('MEDIC single-spatial-resampling outputs committed',flush=True)
    finally:
        lock.unlink()

def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',nargs='+',required=True)
    for key in ['dmap','axis','timing','base-warp','reference','mask','output-dir']:
        p.add_argument('--'+key,required=True)
    p.add_argument('--motion-dir');p.add_argument('--discard',type=int,default=0)
    p.add_argument('--tr',type=float,required=True);p.add_argument('--prefix',default='Rest')
    p.add_argument('--reorient',action='store_true');p.add_argument('--preflight-only',action='store_true')
    return p
if __name__=='__main__':
    args=parser().parse_args()
    if not args.preflight_only and not args.motion_dir:raise SystemExit('--motion-dir is required for resampling')
    run(args)
