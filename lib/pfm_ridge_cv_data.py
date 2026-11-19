from pathlib import Path
import json,math
import numpy as np
import nibabel as nib
from nibabel.cifti2.cifti2_axes import LabelAxis

def make_splits(scan_info,fd,frames,tr,fd_threshold,gap_seconds,min_half_frames,min_partition_minutes):
    if not math.isfinite(tr) or tr<=0:raise ValueError('Input TR must be positive')
    if not math.isfinite(gap_seconds) or gap_seconds<0:raise ValueError('CV guard seconds must be finite and nonnegative')
    if min_half_frames<3 or not math.isfinite(min_partition_minutes) or min_partition_minutes<=0:raise ValueError('CV minimum half frames must be >=3 and minimum partition minutes >0')
    records=json.loads(Path(scan_info).read_text())
    if not isinstance(records,list) or len(records)<2:raise ValueError('Crossed-half CV requires at least two runs in ScanInfo.json')
    lengths=[int(r['n_timepoints']) for r in records]
    if any(n<1 for n in lengths):raise ValueError('ScanInfo run lengths must be positive')
    fd=np.asarray(fd,float).reshape(-1)
    if len(fd)!=sum(lengths):raise ValueError('FD rows do not match acquisition run lengths')
    retained=np.isfinite(fd)&(fd<fd_threshold)
    if frames==len(fd):input_raw=np.arange(len(fd));censored=False
    elif frames==int(retained.sum()):input_raw=np.flatnonzero(retained);censored=True
    else:raise ValueError('Input frame count matches neither uncensored nor FD-censored metadata; supply aligned ScanInfo/FD')
    # For uncensored input the input's actual frame set is used; no unrequested
    # censoring is silently applied. The raw timeline always defines the halves.
    present=np.zeros(len(fd),bool);present[input_raw]=True
    gap=int(math.ceil(gap_seconds/tr));membership=np.full(len(fd),-1,int);runid=np.full(len(fd),-1,int)
    runs=[];start=0
    for i,(n,record) in enumerate(zip(lengths,records)):
        ix=np.arange(start,start+n);local=np.arange(n);mid=n//2;first=local<mid-gap;second=local>=mid+gap;runid[ix]=i
        first_count=int(np.sum(first&present[ix]));second_count=int(np.sum(second&present[ix]));eligible=min(first_count,second_count)>=min_half_frames
        if eligible:
            membership[ix[first]]=i%2;membership[ix[second]]=1-i%2
        runs.append(dict(run_index=i,record=record,raw_start=start,raw_frames=n,midpoint=mid,first_retained=first_count,second_retained=second_count,eligible=eligible,reason=None if eligible else 'one_or_both_guarded_halves_too_short'))
        start+=n
    if sum(r['eligible'] for r in runs)<2:raise ValueError('CV needs at least two eligible runs after guard intervals and censoring')
    parts={}
    for label,tag in enumerate(('A','B')):
        rows=np.flatnonzero(membership[input_raw]==label);minutes=len(rows)*tr/60
        if minutes<min_partition_minutes:raise ValueError(f'CV partition {tag} has {minutes:.2f} retained minutes; needs {min_partition_minutes:g}. Reduce the declared minimum or acquire more data.')
        parts[tag]=dict(rows=rows,raw_rows=input_raw[rows],groups=runid[input_raw[rows]],frames=len(rows),minutes=minutes)
    assert not np.intersect1d(parts['A']['rows'],parts['B']['rows']).size
    return dict(partitions=parts,runs=runs,input_raw_rows=input_raw,raw_membership=membership,input_censored=censored,gap_frames_each_side=gap,gap_seconds_each_side=gap*tr,tr=tr)


def surface_areas(axis,left,right):
    values=np.zeros(len(axis),float)
    for structure,path in [('CIFTI_STRUCTURE_CORTEX_LEFT',left),('CIFTI_STRUCTURE_CORTEX_RIGHT',right)]:
        surf=nib.load(str(path));coords=np.asarray(surf.darrays[0].data,float);faces=np.asarray(surf.darrays[1].data,int)
        triangle=np.linalg.norm(np.cross(coords[faces[:,1]]-coords[faces[:,0]],coords[faces[:,2]]-coords[faces[:,0]]),axis=1)/2
        areas=np.zeros(len(coords))
        for j in range(3):np.add.at(areas,faces[:,j],triangle/3)
        mask=axis.name==structure;vertices=axis.vertex[mask]
        if len(vertices) and (vertices.min()<0 or vertices.max()>=len(areas)):raise ValueError('Surface vertices do not align with CIFTI')
        values[mask]=areas[vertices]
    if not np.isfinite(values).all() or np.any(values<=0):raise ValueError('Cortical surface areas must be finite and positive')
    return np.ascontiguousarray(values/values.mean())


def write_labels(path,labels,names,table,axis):
    hdr=nib.Cifti2Header.from_axes((LabelAxis(names,[table]*len(names)),axis));out=nib.Cifti2Image(np.asarray(labels,np.int32),hdr);out.set_data_dtype(np.int32);nib.save(out,path)



def plan_comparisons(split, minimum_minutes, maximum):
    """Disjoint run groups; alternate early/late WITHIN each comparison."""
    if minimum_minutes <= 0 or maximum < 1:
        raise ValueError('Minimum partition duration and maximum comparisons must be positive')
    parts=split['partitions']
    groups=sorted(set(parts['A']['groups']) & set(parts['B']['groups']))
    for n in range(min(int(maximum),len(groups)//2),0,-1):
        buckets=[groups[i::n] for i in range(n)]
        result=[]
        for bucket in buckets:
            pieces={side:[] for side in ['A','B']}
            for local,g in enumerate(bucket):
                # make_splits alternates by original run index. Rebalance the
                # polarity after grouping, including even-stride run groups.
                swap=(local % 2)!=(int(g) % 2)
                for side,other in [('A','B'),('B','A')]:
                    p=parts[other if swap else side]
                    keep=np.asarray(p['groups'])==g
                    pieces[side].extend(np.asarray(p['rows'])[keep].tolist())
            pair={}
            for side in ['A','B']:
                rows=np.asarray(sorted(pieces[side]),np.int64)
                raw=np.asarray(split['input_raw_rows'])[rows]
                run_groups=np.empty(len(raw),np.int32)
                for record in split['runs']:
                    take=(raw>=record['raw_start'])&(raw<record['raw_start']+record['raw_frames'])
                    run_groups[take]=record['run_index']
                pair[side]=dict(rows=rows,raw_rows=raw,groups=run_groups,
                                frames=len(rows),minutes=len(rows)*split['tr']/60)
            result.append(pair)
        if all(min(p['A']['minutes'],p['B']['minutes'])>=minimum_minutes for p in result):
            occupied=np.concatenate([p[s]['rows'] for p in result for s in ['A','B']])
            if len(np.unique(occupied))!=len(occupied):
                raise ValueError('Internal error: overlapping temporal comparisons')
            return result
    raise ValueError('Not enough retained data for a balanced comparison with two runs per A/B pair')


def normalized(data, rows, groups=None):
    """Run-centered split data; global centering for the ordinary all-data fit."""
    x=np.array(data[np.asarray(rows)].T,dtype=float,order='C')
    if groups is None:
        x-=x.mean(1,keepdims=True)
    else:
        groups=np.asarray(groups)
        for g in np.unique(groups):
            take=groups==g
            x[:,take]-=x[:,take].mean(1,keepdims=True)
    sd=x.std(1,keepdims=True)
    valid=np.isfinite(x).all(1)&np.isfinite(sd[:,0])&(sd[:,0]>1e-8)
    x=np.divide(x,sd,out=np.zeros_like(x),where=valid[:,None])
    return x,valid
