"""Isolated RidgePFM beta numerical routines; imports do not perform analysis."""
from pathlib import Path
import ctypes,subprocess,hashlib,shutil
import numpy as np
from scipy.sparse import csr_matrix

def normalize_blocks(data,groups):
    """Demean separately inside each selected run segment, then scale per vertex."""
    z=np.array(data,dtype=float,copy=True)
    for g in np.unique(groups):
        rows=groups==g;z[:,rows]-=z[:,rows].mean(1,keepdims=True)
    sd=z.std(1,keepdims=True)
    if np.any(sd<=1e-8) or not np.isfinite(z).all():raise ValueError('Invalid cortical time series')
    return np.ascontiguousarray(z/sd)

def solve_penalized(grams,projection,lam,multipliers):
    """Joint fit with penalty lam*diag(multipliers), including minimum-norm OLS."""
    if lam==0:
        e,q=np.linalg.eigh(grams);threshold=np.maximum(e.max(1,keepdims=True),1)*1e-10
        inv=np.divide(1,e,out=np.zeros_like(e),where=e>threshold)
        return np.einsum('bij,bj->bi',q,np.einsum('bij,bi->bj',q,projection)*inv)
    return np.linalg.solve(grams+np.diag(lam*np.asarray(multipliers)),projection[...,None])[...,0]

def eigensystem(grams,projection,multipliers):
    d=1/np.sqrt(multipliers)
    e,q=np.linalg.eigh(grams*d[None,:,None]*d[None,None,:])
    return e,q,np.einsum('bij,bi->bj',q,projection*d),d,grams,projection

def solve_eigen(system,lam):
    e,q,b,d,grams,projection=system
    if lam==0 and not np.array_equal(d,np.ones(len(d))):
        # Weighted pseudoinverses pick different solutions for singular designs.
        # All zero-penalty controls use the same Euclidean minimum-norm OLS.
        return solve_penalized(grams,projection,0,np.ones(len(d)))
    threshold=np.maximum(e.max(1,keepdims=True),1)*1e-10
    inv=np.divide(1,e+lam,out=np.zeros_like(e),where=(e>threshold) if lam==0 else np.ones(e.shape,bool))
    return np.einsum('bij,bj->bi',q,b*inv)*d

def kernel_library(directory):
    source=Path(__file__).with_name('pfm_ridge_cv_kernels.cpp');digest=hashlib.sha256(source.read_bytes()).hexdigest()[:16]
    path=Path(directory)/('cv_kernels_'+digest+'.so')
    if not path.exists():
        if not shutil.which('g++'):return None
        try:subprocess.run(['g++','-O3','-std=c++11','-fopenmp','-shared','-fPIC',str(source),'-o',str(path)],check=True,capture_output=True)
        except subprocess.CalledProcessError:return None
    try:lib=ctypes.CDLL(str(path))
    except OSError:return None
    d=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');i=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS');p=np.ctypeslib.ndpointer(dtype=np.int64,flags='C_CONTIGUOUS')
    lib.edge_correlations.argtypes=[d,ctypes.c_int,ctypes.c_int,p,i,d];lib.edge_correlations.restype=None
    lib.same_network_sums.argtypes=[i,ctypes.c_int,ctypes.c_int,ctypes.c_int,p,i,d,d,d,d];lib.same_network_sums.restype=None
    return lib

def local_correlations(lib,z,ptr,indices):
    out=np.zeros(len(indices),float)
    if lib is not None:lib.edge_correlations(z,len(z),z.shape[1],ptr,indices,out)
    else:
        for i in range(len(z)):out[ptr[i]:ptr[i+1]]=z[indices[ptr[i]:ptr[i+1]]]@z[i]/z.shape[1]
    return out

def evaluate_maps(lib,labels,z,area,ptr,indices,edge_corr,k):
    """Exact area-weighted Pearson pair averages, omitting <=10 mm edges.

    Each network's contrast is mean within-network FC minus mean FC to all
    other assigned cortical vertices, evaluated on frozen map labels.
    """
    labels=np.ascontiguousarray(labels,dtype=np.int32);nm,n=labels.shape;nt=z.shape[1]
    same=np.zeros((nm,k));samew=np.zeros((nm,k))
    if lib is not None:lib.same_network_sums(labels,nm,n,k,ptr,indices,edge_corr,area,same,samew)
    else:
        source=np.repeat(np.arange(n),np.diff(ptr));pairweight=area[source]*area[indices]
        for m,lab in enumerate(labels):
            keep=lab[source]==lab[indices];groups=lab[source[keep]]-1
            same[m]=np.bincount(groups,weights=edge_corr[keep]*pairweight[keep],minlength=k)
            samew[m]=np.bincount(groups,weights=pairweight[keep],minlength=k)
    adj=csr_matrix((np.ones(len(indices)),indices,ptr),shape=(n,n))
    edges=csr_matrix((edge_corr,indices,ptr),shape=(n,n))
    totalnear=area*(edges@area);totalnearw=area*(adj@area)
    totalts=area@z;totalarea=area.sum();out=[]
    for m,lab in enumerate(labels):
        indicator=csr_matrix((area,(lab-1,np.arange(n))),shape=(k,n))
        sums=indicator@z;weights=np.bincount(lab-1,weights=area,minlength=k)
        within_num=np.sum(sums*sums,axis=1)/nt-same[m]
        within_den=weights**2-samew[m]
        allnear=np.bincount(lab-1,weights=totalnear,minlength=k)
        allnearw=np.bincount(lab-1,weights=totalnearw,minlength=k)
        between_num=np.sum(sums*(totalts-sums),axis=1)/nt-(allnear-same[m])
        between_den=weights*(totalarea-weights)-(allnearw-samew[m])
        within=np.divide(within_num,within_den,out=np.zeros(k),where=within_den>1e-10)
        between=np.divide(between_num,between_den,out=np.zeros(k),where=between_den>1e-10)
        contrast=within-between
        # A missing/unevaluable class contributes zero, and is explicitly flagged.
        contrast[(within_den<=1e-10)|(between_den<=1e-10)]=0
        out.append(dict(macro_contrast=float(contrast.mean()),area_contrast=float(np.sum(weights*contrast)/totalarea),macro_within_r=float(within.mean()),area_within_r=float(np.sum(weights*within)/totalarea),within_r=within,between_r=between,contrast=contrast,area=weights,vertices=np.bincount(lab-1,minlength=k),unevaluable=(within_den<=1e-10)|(between_den<=1e-10)))
    return out
