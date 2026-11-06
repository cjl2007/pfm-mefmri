"""MEDIC/warpkit physical-vector to FSL scaled-voxel conversion.

This reproduces the coordinate interpretation of the installed warpkit
resample_image function. It does not reinterpret a MEDIC map as a displacement
along an oblique voxel axis. Never infer a map's sign again from its PE suffix:
MEDIC's map already contains signed displacement values.
"""
import numpy as np
import nibabel as nib

AXES={'i':0,'i-':0,'x':0,'x-':0,'j':1,'j-':1,'y':1,'y-':1,'k':2,'k-':2,'z':2,'z-':2}

def scaled_voxel_affine(image):
    affine=np.asarray(image.affine,dtype=float)
    if affine.shape!=(4,4) or not np.isfinite(affine).all():
        raise ValueError('Invalid physical affine')
    if abs(np.linalg.det(affine[:3,:3]))<1e-8:
        raise ValueError('Singular physical affine')
    spacing=np.asarray(image.header.get_zooms()[:3],dtype=float)
    if np.any(spacing<=0) or not np.isfinite(spacing).all():
        raise ValueError('Invalid voxel spacing')
    # FSL scaled-mm convention flips x for neurological storage.
    out=np.diag([*spacing,1.])
    if np.linalg.det(affine[:3,:3])>0:
        out[0,0]*=-1
        out[0,3]=(image.shape[0]-1)*spacing[0]
    return out

def validate_grid(image,reference):
    if image.shape[:3]!=reference.shape[:3] or not np.allclose(image.affine,reference.affine,rtol=0,atol=1e-4):
        raise ValueError('MEDIC map and reoriented native data grids differ')
    # Warpkit decomposes rotations/scales. Do not silently accept shears, for
    # which that decomposition and FSL's full affine need not agree.
    a=image.affine[:3,:3];u=a/np.linalg.norm(a,axis=0)
    if not np.allclose(u.T@u,np.eye(3),atol=1e-4,rtol=0):
        raise ValueError('Sheared native grid is not supported')

def map_to_fsl(image,reference,axis,frame=0):
    if axis not in AXES:raise ValueError('Unsupported phase-encoding axis')
    validate_grid(image,reference)
    if image.ndim==3:
        if frame!=0:raise ValueError('Static displacement map has only frame 0')
        data=np.asarray(image.dataobj,dtype=np.float32)
    elif image.ndim==4:
        if not 0<=frame<image.shape[3]:raise ValueError('Displacement frame out of range')
        data=np.asarray(image.dataobj[...,frame],dtype=np.float32)
    else:raise ValueError('Displacement map must be 3D or 4D')
    if not np.isfinite(data).all():raise ValueError('Nonfinite MEDIC displacement')
    physical_to_scaled=scaled_voxel_affine(reference)[:3,:3]@np.linalg.inv(reference.affine[:3,:3])
    field=data[...,None]*physical_to_scaled[:,AXES[axis]]
    header=reference.header.copy();header.set_data_dtype(np.float32)
    result=nib.Nifti1Image(field.astype(np.float32),reference.affine,header)
    result.header.set_intent('vector')
    return result

def frame_indices(n_volumes,n_maps,discard=0):
    if n_volumes<1 or not 0<=discard<n_volumes:
        raise ValueError('Invalid volume count/discard offset')
    if n_maps not in (1,n_volumes):
        raise ValueError('MEDIC field count must be one or match original acquisition')
    return [(original,0 if n_maps==1 else original) for original in range(discard,n_volumes)]
