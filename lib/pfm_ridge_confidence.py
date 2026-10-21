#!/usr/bin/env python3
"""Session assignment frequencies and distinct pooled-split FC/label states.

Frequency is empirical session support, not a calibrated biological probability.
The four states use three disjoint session-half pairs (pairs reuse sessions).
"""
from __future__ import annotations
import argparse
import json
import subprocess
import sys
from pathlib import Path
import nibabel as nib
import numpy as np
from scipy.io import loadmat
from pfm_ridge_fusion import zscore_rows
from nibabel.cifti2 import ScalarAxis, LabelAxis


def save_maps(ref, values, names, path, table=None):
    first = ScalarAxis(names) if table is None else LabelAxis(names, [table]*len(names))
    header = nib.Cifti2Header.from_axes((first, ref.header.get_axis(1)))
    nib.save(nib.Cifti2Image(np.asarray(values), header), path)


def load_sessions(path, nframes):
    values = np.asarray(loadmat(path)['ScanIdx'])
    if values.shape != (nframes, 2) or not np.isfinite(values).all() or np.any(values != np.floor(values)) or np.any(values < 1):
        raise ValueError('ScanIdx must have exactly one [session, run] row per input frame; supply censored/aligned metadata')
    return values[:, 0].astype(int)


def session_splits(sessions, seed=20261021):
    ids = np.unique(sessions)
    if len(ids) < 4:
        raise ValueError('Four-state confidence requires at least four eligible sessions (two in each half)')
    rng = np.random.default_rng(seed)
    pairs = []
    for _ in range(3):
        order = rng.permutation(ids)
        pairs.append((order[:len(ids)//2], order[len(ids)//2:]))
    return pairs


def assignment_states(split_labels, correlations, threshold=.7):
    labels = np.asarray(split_labels)
    valid = (labels > 0).all(0) & np.isfinite(correlations).all(0)
    same = (labels == labels[0]).all(0)
    reproducible = np.min(correlations, axis=0) >= threshold
    states = np.zeros(labels.shape[1], np.int32)
    states[valid] = (1 + (~same[valid]).astype(int) + 2*(~reproducible[valid]).astype(int))
    return states


def profile_correlation(a, b, keep):
    a = np.where(keep, a, 0); b = np.where(keep, b, 0)
    n = keep.sum(0)
    cross = (a*b).sum(0) - a.sum(0)*b.sum(0)/np.maximum(n, 1)
    va = (a*a).sum(0)-a.sum(0)**2/np.maximum(n, 1)
    vb = (b*b).sum(0)-b.sum(0)**2/np.maximum(n, 1)
    den = np.sqrt(np.maximum(va*vb, 0))
    return np.divide(cross, den, out=np.full(a.shape[1], np.nan), where=(den > 1e-12)&(n > 2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--in-cifti', required=True); p.add_argument('--scan-idx', required=True)
    p.add_argument('--full-map-dlabel', required=True); p.add_argument('--distance-npy', required=True)
    p.add_argument('--priors-mat', required=True); p.add_argument('--left-surf', required=True); p.add_argument('--right-surf', required=True)
    p.add_argument('--outdir', required=True); p.add_argument('--outfile', default='RidgePFM_VTX')
    p.add_argument('--minimum-frames', type=int, default=30); p.add_argument('--fc-threshold', type=float, default=.7)
    p.add_argument('--seed', type=int, default=20261021)
    p.add_argument('--fc-features', type=int, default=0, help='0: exact all-cortex profiles; positive: seeded feature subset, recorded as sampled')
    p.add_argument('--local-exclusion-mm', type=float, default=10.0)
    p.add_argument('--fc-cache', help='Optional shared split-profile cache directory across taxonomies')
    p.add_argument('--split-plan', help='JSON list of three explicit [train sessions, reference sessions] pairs')
    args, model_args = p.parse_known_args()
    if not -1 <= args.fc_threshold <= 1 or args.fc_features < 0: p.error('Invalid FC threshold/features')
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    ref = nib.load(args.in_cifti); full = nib.load(args.full_map_dlabel)
    if full.header.get_axis(1) != ref.header.get_axis(1): raise ValueError('WTA geometry mismatch')
    sessions = load_sessions(args.scan_idx, ref.shape[0])
    eligible = [int(s) for s in np.unique(sessions) if np.sum(sessions == s) >= args.minimum_frames]
    splits = session_splits(np.asarray(eligible), args.seed)
    if args.split_plan:
        splits = [(np.asarray(a,int),np.asarray(b,int)) for a,b in json.loads(Path(args.split_plan).read_text())]
        if len(splits) != 3 or any(len(a)<2 or len(b)<2 or len(set(a)&set(b)) or set(a)|set(b)!=set(eligible) for a,b in splits):
            raise ValueError('Each of three split pairs must partition all eligible sessions with at least two per side')
    pri = loadmat(args.priors_mat, squeeze_me=True, struct_as_record=False)['Priors']
    names = [str(x) for x in np.asarray(pri.NetworkLabels).ravel()[:pri.FC.shape[1]]]
    table = full.header.get_axis(0).label[0]
    work = out/'resamples'; work.mkdir(exist_ok=True)
    base = [sys.executable, str(Path(__file__).with_name('pfm_ridge_fusion.py')), '--in-cifti', args.in_cifti,
            '--distance-npy', args.distance_npy, '--priors-mat', args.priors_mat,
            '--left-surf', args.left_surf, '--right-surf', args.right_surf,
            '--local-exclusion-mm', str(args.local_exclusion_mm), '--fit-r2-mode', ('exact' if any(flag in model_args and model_args[model_args.index(flag)+1] == '1' for flag in ['--fc-demean','--positive-fc-only']) else 'skip')] + model_args
    # Explicitly record command and frame identities. Never silently reuse stale maps.
    def fit(tag, frames):
        dest = work/tag; dest.mkdir(exist_ok=True)
        fp = dest/'frames.npy'; np.save(fp, frames)
        cmd = base + ['--frame-indices', str(fp), '--outdir', str(dest), '--outfile', 'map']
        subprocess.run(cmd, check=True)
        result = np.asarray(nib.load(dest/'map.dlabel.nii').dataobj)[0].astype(np.int32)
        return result
    session_labels = np.stack([fit('session_%04d'%s, np.flatnonzero(sessions == s)) for s in eligible])
    counts = np.stack([(session_labels == k).sum(0) for k in range(1,len(names)+1)])
    valid_count = (session_labels > 0).sum(0)
    frequency = np.divide(counts, valid_count, out=np.full(counts.shape,np.nan), where=valid_count > 0).astype(np.float32)
    save_maps(ref, frequency, names, out/(args.outfile+'_SessionAssignmentFrequency.dscalar.nii'))
    save_maps(ref, session_labels, ['Session %s'%s for s in eligible], out/(args.outfile+'_SessionAssignments.dlabel.nii'), table)
    split_labels = []
    for f, pair in enumerate(splits):
        for side, ids in enumerate(pair):
            split_labels.append(fit('split_%d_%d'%(f,side), np.flatnonzero(np.isin(sessions,ids))))
    split_labels = np.stack(split_labels)
    save_maps(ref, split_labels, ['Split %d half %d'%(f+1,s+1) for f in range(3) for s in range(2)], out/(args.outfile+'_SplitAssignments.dlabel.nii'), table)
    import hashlib
    source = Path(args.in_cifti).resolve(); st = source.stat()
    signature = dict(source=str(source), bytes=st.st_size, mtime_ns=st.st_mtime_ns,
                     sessions=hashlib.sha256(sessions.tobytes()).hexdigest(), splits=[[x.tolist() for x in pair] for pair in splits],
                     features=args.fc_features, seed=args.seed, exclusion=args.local_exclusion_mm,
                     distance=str(Path(args.distance_npy).resolve()), distance_mtime=Path(args.distance_npy).stat().st_mtime_ns)
    digest = hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()
    cache = Path(args.fc_cache) if args.fc_cache else out/'fc_profiles'; cache.mkdir(parents=True,exist_ok=True)
    cp = cache/(digest+'.npy')
    if cp.exists():
        correlations = np.load(cp)
    else:
        axis = ref.header.get_axis(1); cortex = np.flatnonzero(axis.surface_mask)
        features = cortex if not args.fc_features else np.sort(np.random.default_rng(args.seed).choice(cortex,min(len(cortex),args.fc_features),replace=False))
        data = np.asanyarray(ref.dataobj); dist = np.load(args.distance_npy,mmap_mode='r'); correlations=[]
        for f,pair in enumerate(splits):
            rows = [np.flatnonzero(np.isin(sessions,ids)) for ids in pair]
            z = [zscore_rows(data[r].T.astype(np.float32)) for r in rows]
            fc = np.full(ref.shape[1],np.nan,np.float32)
            for start in range(0,ref.shape[1],256):
                stop=min(start+256,ref.shape[1]); target=np.arange(start,stop)
                profiles=[v[features] @ v[start:stop].T / v.shape[1] for v in z]
                keep=np.ones(profiles[0].shape,bool)
                for j,g in enumerate(target):
                    if axis.surface_mask[g]:keep[:,j]=dist[g,features] > args.local_exclusion_mm
                fc[start:stop]=profile_correlation(profiles[0].astype(float),profiles[1].astype(float),keep)
            correlations.append(fc); print('FC profile split',f+1,'complete',flush=True)
        correlations=np.stack(correlations);np.save(cp,correlations)
    states = assignment_states(split_labels,correlations,args.fc_threshold)
    valid_session_counts = np.stack([
        np.sum(session_labels[np.isin(eligible, ids)] > 0, axis=0)
        for pair in splits for ids in pair
    ])
    minimum_valid_sessions = valid_session_counts.min(0)
    states[minimum_valid_sessions < 2] = 0
    save_maps(ref, minimum_valid_sessions[None].astype(np.float32), ['Minimum valid sessions per split side'],
              out/(args.outfile+'_MinimumValidSplitSessions.dscalar.nii'))
    save_maps(ref, correlations, ['FC profile r split %d'%(i+1) for i in range(3)],out/(args.outfile+'_FCReproducibility.dscalar.nii'))
    state_names=['Not estimable','Reproducible FC; same label','Reproducible FC; competing labels','Lower FC reproducibility; same label','Lower FC reproducibility; competing labels']
    colors=[(0,0,0,0),(.13,.42,.41,1),(.41,.72,.67,1),(.83,.65,.35,1),(.60,.36,.25,1)]
    save_maps(ref,states[None],['Assignment state'],out/(args.outfile+'_AssignmentStates.dlabel.nii'),dict(enumerate(zip(state_names,colors))))
    # Spatial stripes cycle through the complete split-label set, only in state 2.
    stripe=np.asarray(full.dataobj)[0].astype(np.int32).copy();axis=ref.header.get_axis(1)
    coords=np.zeros((len(axis),3))
    for name,sl,bm in axis.iter_structures():
        if name.endswith('CORTEX_LEFT'):coords[sl]=nib.load(args.left_surf).darrays[0].data[bm.vertex]
        elif name.endswith('CORTEX_RIGHT'):coords[sl]=nib.load(args.right_surf).darrays[0].data[bm.vertex]
        else:coords[sl]=nib.affines.apply_affine(axis.affine,bm.voxel)
    bands=np.floor((coords[:,0]+coords[:,1]+coords[:,2])/2).astype(int)
    membership=np.stack([(split_labels==k).any(0) for k in range(1,len(names)+1)])
    for g in np.flatnonzero(states==2):
        choices=np.flatnonzero(membership[:,g])+1;stripe[g]=choices[bands[g]%len(choices)]
    save_maps(ref,stripe[None],['Reproducible competing-label stripes'],out/(args.outfile+'_NetworkStripes.dlabel.nii'),table)
    save_maps(ref,membership.astype(np.float32),names,out/(args.outfile+'_SplitLabelMembership.dscalar.nii'))
    (out/'confidence.json').write_text(json.dumps(dict(signature=signature,sessions=eligible,excluded_sessions=[int(s) for s in np.unique(sessions) if s not in eligible],
        fc_threshold=args.fc_threshold,fc_features=args.fc_features,states=state_names,
        interpretation='Empirical frequency across independently mapped sessions; not calibrated probability of biological membership. Four states use three pooled-session splits; stripes only state 2. Stripe widths are not probabilities. Model fit remains separate.'),indent=2)+'\n')

if __name__ == '__main__': main()
