"""Independent pull-coordinate check for FSL displacement composition.

All supplied fields must be explicitly relative, dense FSL displacement fields.
This evaluates the chain with NumPy/SciPy; it does not compare two applywarp calls.
"""
import numpy as np
import nibabel as nib
from scipy.ndimage import map_coordinates
from medic_geometry import scaled_voxel_affine

def check_chain(reference,native,distortion,intermediate_to_ref,motion,combined,mask):
    for im in [intermediate_to_ref,combined]:
        if im.shape!=(*reference.shape[:3],3) or not np.allclose(im.affine,reference.affine,atol=1e-4):
            raise ValueError('Anatomical field does not match reference grid')
    if distortion.shape!=(*native.shape[:3],3):raise ValueError('Native field grid mismatch')
    r=scaled_voxel_affine(reference);s=scaled_voxel_affine(native)
    # Deterministic spatial sample, including all octants and brain boundaries.
    indices=np.argwhere(mask)[::13]
    xyz=np.c_[indices,np.ones(len(indices))].T
    target=(r@xyz)[:3]
    w=intermediate_to_ref.get_fdata(dtype=np.float32)
    intermediate=target+w[tuple(indices.T)].T
    q=(np.linalg.inv(motion)@np.vstack([intermediate,np.ones(len(indices))]))[:3]
    voxel=(np.linalg.inv(s)@np.vstack([q,np.ones(len(indices))]))[:3]
    valid=np.all((voxel>=1)&(voxel<np.array(native.shape[:3])[:,None]-2),axis=0)
    if valid.sum()<100:raise ValueError('Insufficient interior points for chain validation')
    field=distortion.get_fdata(dtype=np.float32)
    offset=np.array([map_coordinates(field[...,i],voxel[:,valid],order=1,mode='nearest',prefilter=False) for i in range(3)])
    expected=q[:,valid]+offset
    actual=target[:,valid]+combined.get_fdata(dtype=np.float32)[tuple(indices[valid].T)].T
    error=np.linalg.norm(expected-actual,axis=0)
    result=dict(points=int(valid.sum()),max_mm=float(error.max()),rms_mm=float(np.sqrt(np.mean(error**2))))
    if result['max_mm']>.01:raise ValueError('Composed pull-coordinate error: '+str(result))
    return result
