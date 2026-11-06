"""Prototype joint temporal/spatial sampling. Not installed in production.

Each temporal contributor is sampled once from its native intensity frame with
its own MEDIC field. Every native slice contributing to the spatial kernel keeps
its own temporal filter. In-plane cubic B-spline interpolation is separable from
four-tap interpolating cubic slice weights: there is no intermediate resampled
intensity volume, nor a spline prefilter mixing slices with different times.
This is not byte-equivalent to the historical unwarp-then-STC ordering.
"""
from pathlib import Path
import os,subprocess,shutil
import numpy as np
import nibabel as nib
from scipy.ndimage import map_coordinates,spline_filter
from medic_geometry import AXES

def fsl_temporal_weights(n,offsets,folder):
    offsets=np.asarray(offsets,dtype=float)
    if n<1 or offsets.ndim!=1 or not np.isfinite(offsets).all() or np.max(abs(offsets))>1:
        raise ValueError('Invalid temporal dimensions or offsets')
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    timing=folder/'offsets.txt';timing.write_text('\n'.join(map(str,offsets))+'\n')
    # Each x-position is a different temporal impulse. Measure the installed
    # FSL operator, including edge renormalization and lookup-table rounding.
    impulses=np.broadcast_to(np.eye(n,dtype=np.float32)[:,None,None,:],(n,1,len(offsets),n)).copy()
    src=folder/'impulses.nii.gz';out=folder/'impulses_stc.nii.gz'
    nib.save(nib.Nifti1Image(impulses,np.eye(4)),src)
    env=os.environ.copy();env['FSLOUTPUTTYPE']='NIFTI_GZ'
    executable=shutil.which('slicetimer')
    if executable is None:
        candidate=Path(os.environ.get('FSLDIR','/usr/local/fsl'))/'bin/slicetimer'
        if not candidate.is_file():raise FileNotFoundError('FSL slicetimer is required on PATH or under FSLDIR/bin')
        executable=str(candidate)
    subprocess.run([executable,'-i',str(src),'-o',str(out),'-r','1','--tcustom='+str(timing)],check=True,capture_output=True,env=env)
    # (slice, output-time, input-time)
    weights=nib.load(out).get_fdata(dtype=np.float32)[:,0,:,:].transpose(1,2,0)
    if not np.allclose(weights.sum(axis=2),1,atol=2e-6):raise ValueError('Temporal weights do not preserve constants')
    return weights

def sample_joint_frame(raw,maps,affine,axis,q,weights,t,coefficients=None):
    """q: callable(k)->3xN unwarped native coordinates for acquisition k,
    or a fixed 3xN array for tests with a time-invariant spatial transform.

    Motion/atlas inverse chain must already have been evaluated in q(k). The
    field for every contributing acquisition frame k is evaluated at q(k).
    """
    if axis not in AXES:raise ValueError('Invalid PE axis')
    if raw.ndim!=4 or maps.shape[:3]!=raw.shape[:3]:raise ValueError('Input grid mismatch')
    if maps.ndim not in (3,4) or (maps.ndim==4 and maps.shape[3] not in (1,raw.shape[3])):raise ValueError('Map frame mismatch')
    if weights.shape!=(raw.shape[2],raw.shape[3],raw.shape[3]):raise ValueError('Temporal filter shape mismatch')
    if not 0<=t<raw.shape[3]:raise ValueError('Output frame index out of bounds')
    get_coordinates=q if callable(q) else lambda k:q
    initial=np.asarray(get_coordinates(t),dtype=float)
    if initial.ndim!=2 or initial.shape[0]!=3 or not np.isfinite(initial).all():raise ValueError('Invalid pull coordinates')
    direction=np.linalg.inv(np.asarray(affine)[:3,:3])[:,AXES[axis]]
    result=np.zeros(initial.shape[1],dtype=np.float64)
    normalization=np.zeros(initial.shape[1],dtype=np.float64)
    contributors=np.flatnonzero(np.any(weights[:,t,:]!=0,axis=0))
    for k in contributors:
        coordinates=np.asarray(get_coordinates(k),dtype=float)
        if coordinates.shape!=initial.shape or not np.isfinite(coordinates).all():raise ValueError('Contributor coordinate mismatch')
        valid=np.all((coordinates>=0)&(coordinates<=np.array(raw.shape[:3])[:,None]-1),axis=0)
        field=maps if maps.ndim==3 else maps[...,0 if maps.shape[3]==1 else k]
        displacement=map_coordinates(field,coordinates,order=1,mode='nearest',prefilter=False)
        source=coordinates+direction[:,None]*displacement[None,:]
        center=np.floor(source[2]).astype(int)
        order=np.argsort(center);splits=np.flatnonzero(np.diff(center[order]))+1
        groups=np.split(order,splits)
        for group in groups:
            base=center[group[0]]
            for offset in [-1,0,1,2]:
                s=base+offset
                if not 0<=s<raw.shape[2]:continue
                temporal=weights[s,t,k]
                if temporal==0:continue
                distance=np.abs(source[2,group]-s)
                # Keys interpolating cubic, a=-0.5 (Catmull-Rom).
                zweight=np.where(distance<=1,1.5*distance**3-2.5*distance**2+1,
                                 np.where(distance<2,-.5*distance**3+2.5*distance**2-4*distance+2,0))
                weight=zweight*temporal
                c=coefficients(k,s) if coefficients else spline_filter(raw[:,:,s,k],order=3,mode='mirror')
                sampled=map_coordinates(c,source[:2,group],order=3,mode='constant',cval=0,prefilter=False)
                in_xy=np.all((source[:2,group]>=0)&(source[:2,group]<=np.array(raw.shape[:2])[:,None]-1),axis=0)
                result[group]+=weight*sampled*valid[group]
                normalization[group]+=weight*in_xy*valid[group]
    supported=normalization>.25
    result[supported]/=normalization[supported]
    result[~supported]=0
    return result.astype(np.float32)
