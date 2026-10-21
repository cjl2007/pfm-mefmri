"""Optional compiled exact local-edge kernel; NumPy remains the fallback."""
import ctypes
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import numpy as np


def local_terms(z, design, excluded, cache):
    compiler=shutil.which('g++')
    if compiler is None:return None
    source=Path(__file__).with_suffix('.cpp')
    digest=hashlib.sha256(source.read_bytes()).hexdigest()[:20]
    cache=Path(cache);cache.mkdir(parents=True,exist_ok=True);library=cache/('local_'+digest+'.so')
    if not library.exists():
        fd,temporary=tempfile.mkstemp(suffix='.so',dir=cache);os.close(fd)
        try:
            subprocess.run([compiler,'-O3','-std=c++11','-fopenmp','-shared','-fPIC',str(source),'-o',temporary],check=True,capture_output=True)
            os.replace(temporary,library)
        except subprocess.CalledProcessError:
            Path(temporary).unlink(missing_ok=True)
            return None
    z=np.ascontiguousarray(z,dtype=np.float64);a=np.ascontiguousarray(design,dtype=np.float64)
    ptr=np.r_[0,np.cumsum([len(x) for x in excluded])].astype(np.int64)
    ids=np.concatenate(excluded).astype(np.int32)
    c=np.zeros((len(z),a.shape[1]),np.float64);g=np.zeros((len(z),a.shape[1],a.shape[1]),np.float64)
    d=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');i=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS');q=np.ctypeslib.ndpointer(dtype=np.int64,flags='C_CONTIGUOUS')
    try:
        lib=ctypes.CDLL(str(library))
    except OSError:
        return None
    lib.ridge_local.argtypes=[d,ctypes.c_int,ctypes.c_int,q,i,d,ctypes.c_int,d,d];lib.ridge_local.restype=None
    lib.ridge_local(z,len(z),z.shape[1],ptr,ids,a,a.shape[1],c,g)
    return c,g
