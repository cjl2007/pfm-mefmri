#!/usr/bin/env python3
"""Optional blocked CV of shared RidgePFM lambda/spatial parameters.

All outputs and signature-bound caches belong to --outdir. Input data and
resources are read-only. The workflow module applies the selected scalar pair
through its ordinary mapper after this command succeeds.
"""
from __future__ import annotations
import os
for _key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[_key]=os.environ.get('PFM_RF_CV_THREADS','4')
import argparse,csv,gc,hashlib,json,math,shutil,sys,time
from pathlib import Path
import numpy as np
import nibabel as nib
from scipy.io import loadmat
from scipy.sparse import csr_matrix
from pfm_ridge_cv_core import kernel_library,local_correlations
from pfm_ridge_fusion import clamp_prob_smooth,spatial_prior_term
from pfm_ridge_cv_data import make_splits,plan_comparisons,normalized,surface_areas,write_labels
from pfm_ridge_cv_metrics import prepare,evaluate,repeatability,aggregate_directions
from pfm_ridge_cv_selection import choose

SCHEMA='ridgepfm_fast_blocked_cv_v2'

def clean(value):
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    if isinstance(value,np.ndarray):return clean(value.tolist())
    if isinstance(value,np.generic):return clean(value.item())
    if isinstance(value,float) and not np.isfinite(value):return None
    return value

def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(clean(value),indent=2,allow_nan=False)+'\n');tmp.replace(path)

def save_npz(path,**values):
    path=Path(path);tmp=path.with_name(path.stem+'.tmp.npz');np.savez_compressed(tmp,**values);tmp.replace(path)

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024**2),b''):h.update(chunk)
    return h.hexdigest()

def positive_vector(text,name):
    try:values=sorted(set(float(v) for v in text.replace(',',' ').split()))
    except ValueError as e:raise ValueError(name+' must be a vector of positive numbers') from e
    if not values or any(not math.isfinite(v) or v<=0 for v in values):
        raise ValueError(name+' must contain finite positive numbers; zero spatial weight is excluded from CV')
    return values

def parser():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['in-cifti','distance-npy','priors-mat','scan-info-json','fd-txt','left-surf','right-surf','rotations-cifti','outdir']:
        p.add_argument('--'+name,required=True)
    p.add_argument('--subject',default='Subject');p.add_argument('--surface-dir',default='');p.add_argument('--flat-dir',default='')
    p.add_argument('--lambdas',default='0.01,0.03,0.1,0.3,1,3,10,30')
    p.add_argument('--spatial-weights',default='0.01,0.025,0.05,0.1,0.25,0.5')
    p.add_argument('--max-fc-loss',type=float,default=.1)
    p.add_argument('--minimum-partition-minutes',type=float,default=10)
    p.add_argument('--max-comparisons',type=int,default=10)
    p.add_argument('--n-rotations',type=int,default=100)
    p.add_argument('--guard-seconds',type=float,default=32)
    p.add_argument('--fd-threshold',type=float,default=.3);p.add_argument('--tr-seconds',type=float)
    p.add_argument('--fc-weight',type=float,default=1);p.add_argument('--fc-demean',type=int,choices=[0,1],default=0)
    p.add_argument('--positive-fc-only',type=int,choices=[0,1],default=0)
    p.add_argument('--spatial-prior-transform',choices=['logit','log'],default='logit')
    p.add_argument('--local-exclusion-mm',type=float,default=10);p.add_argument('--brain-structures-csv',default='')
    p.add_argument('--taxonomy',default='Lynch2024');p.add_argument('--render',type=int,choices=[0,1],default=1)
    p.add_argument('--keep-cache',type=int,choices=[0,1],default=0)
    p.add_argument('--preflight-only',action='store_true')
    return p

def preflight(args):
    args.lambda_values=positive_vector(args.lambdas,'lambda');args.weight_values=positive_vector(args.spatial_weights,'spatial weights')
    if not math.isfinite(args.max_fc_loss) or not 0<=args.max_fc_loss<=1:raise ValueError('FC loss must be between 0 and 1')
    if args.n_rotations<1 or args.max_comparisons<1:raise ValueError('Rotation and comparison counts must be positive')
    if not math.isfinite(args.fc_weight) or args.fc_weight<=0:raise ValueError('CV requires positive FC weight')
    if not math.isfinite(args.local_exclusion_mm) or args.local_exclusion_mm<0:raise ValueError('Local exclusion must be finite and nonnegative')
    if not math.isfinite(args.fd_threshold) or args.fd_threshold<=0:raise ValueError('FD threshold must be positive')
    if args.positive_fc_only:raise ValueError('CV beta supports signed FC; disable CV or set PFM_RF_POSITIVE_FC_ONLY=0')
    if args.taxonomy!='Lynch2024':raise ValueError('CV beta is validated for Lynch2024 parent cortical networks; other taxonomies remain available in fixed mode')
    requested={v.strip().upper().replace('CIFTI_STRUCTURE_','') for v in args.brain_structures_csv.split(',') if v.strip()}
    if requested and not {'CORTEX_LEFT','CORTEX_RIGHT'}.issubset(requested):raise ValueError('CV requires both cortical hemispheres in the requested structures')
    image=nib.load(args.in_cifti);axis=image.header.get_axis(1);ci=np.flatnonzero(axis.surface_mask)
    if not len(ci) or not np.array_equal(ci,np.arange(len(ci))):raise ValueError('CV beta requires the cortical models first in the CIFTI axis')
    records=json.loads(Path(args.scan_info_json).read_text())
    metadata_tr=[float(r['fd']['tr_seconds']) for r in records if r.get('fd',{}).get('tr_seconds') is not None]
    if metadata_tr and not np.allclose(metadata_tr,metadata_tr[0],atol=1e-8,rtol=0):raise ValueError('Mixed run TRs are not supported in this combined CIFTI')
    tr=args.tr_seconds if args.tr_seconds is not None else metadata_tr[0] if metadata_tr else float(image.header.get_axis(0).step)
    if not math.isfinite(tr) or tr<=0:raise ValueError('TR must be finite and positive')
    args.actual_tr=tr;args.tr_source='explicit' if args.tr_seconds is not None else 'ScanInfo' if metadata_tr else 'CIFTI header'
    split=make_splits(args.scan_info_json,np.loadtxt(args.fd_txt),image.shape[0],tr,args.fd_threshold,args.guard_seconds,30,args.minimum_partition_minutes)
    pairs=plan_comparisons(split,args.minimum_partition_minutes,args.max_comparisons)
    pri=loadmat(args.priors_mat,squeeze_me=True,struct_as_record=False)['Priors']
    if np.asarray(pri.FC).shape!=(len(ci),20):raise ValueError('CV beta expects the20 Lynch2024 cortical parent templates aligned to input cortex')
    rotation_image=nib.load(args.rotations_cifti)
    if rotation_image.header.get_axis(1)!=axis[ci]:raise ValueError('Rotation cortical axis does not match input cortex')
    if rotation_image.shape[0]<args.n_rotations:raise ValueError('Requested more rotations than the supplied resource contains')
    rotations=np.rint(np.asarray(rotation_image.dataobj)[:args.n_rotations]).astype(np.int64)
    if np.any((rotations<0)|(rotations>len(ci))):raise ValueError('Invalid rotation indices')
    for rot in rotations:
        valid=rot>0
        if not np.array_equal(axis[ci].name[valid],axis[ci].name[rot[valid]-1]):raise ValueError('Rotation crosses cortical hemispheres')
    if args.render:
        for h in ['L','R']:
            for folder,kind in [(args.flat_dir,'flat')]+[(args.surface_dir,k) for k in ['midthickness','inflated','very_inflated','pial','white']]:
                if not folder or not (Path(folder)/f'{args.subject}.{h}.{kind}.32k_fs_LR.surf.gii').is_file():
                    raise ValueError('Missing Workbench scene surface: '+str(Path(folder)/f'{args.subject}.{h}.{kind}.32k_fs_LR.surf.gii'))
    return image,ci,split,pairs,pri,rotations

def binding(args):
    inputs=[args.in_cifti,args.distance_npy,args.priors_mat,args.scan_info_json,args.fd_txt,args.left_surf,args.right_surf,args.rotations_cifti]
    if args.render:
        inputs += [str(Path(folder)/f'{args.subject}.{h}.{kind}.32k_fs_LR.surf.gii') for h in ['L','R'] for folder,kind in [(args.flat_dir,'flat'),(args.surface_dir,'inflated')]]
    if args.render:inputs.append(str(Path(args.priors_mat).resolve().parents[2]/'CrawlingSeedFC/FlatMaps+Inflated.scene'))
    records=[]
    for value in inputs:
        p=Path(value).resolve();s=p.stat();r=dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
        if s.st_size<100000000:r['sha256']=digest(p)
        records.append(r)
    names=['pfm_ridge_cv.py','pfm_ridge_cv_data.py','pfm_ridge_cv_metrics.py','pfm_ridge_cv_selection.py','pfm_ridge_cv_core.py','pfm_ridge_cv_kernels.cpp','pfm_ridge_cv_report.py','pfm_ridge_cv_report.html','pfm_ridge_cv_workbench.py','pfm_ridge_scene.py','crawling_seed_fc_movie.py','pfm_ridge_fusion.py']
    config={k:v for k,v in vars(args).items() if k not in ['outdir','preflight_only','keep_cache']}
    value=dict(schema=SCHEMA,settings=config,inputs=records,code={n:digest(Path(__file__).with_name(n)) for n in names})
    return hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest()[:20],value

class Fitter:
    def __init__(self,args,pri,nc,ngray,cache):
        self.args=args;self.nc=nc;self.cache=cache
        self.a=np.asarray(pri.FC,np.float32).astype(float);self.a-=self.a.mean(0);self.a/=np.maximum(np.linalg.norm(self.a,axis=0),1e-8)
        self.term=spatial_prior_term(clamp_prob_smooth(np.asarray(pri.Spatial,np.float32)[:nc,:20]),args.spatial_prior_transform)
        self.atlas=(np.argmax(self.term,axis=1)+1).astype(np.uint8)
        self.settings=[dict(key=f'l{l:g}_w{w:g}',lam=l,spatial=w,control=False) for l in args.lambda_values for w in args.weight_values]+[dict(key='prior_only',lam=None,spatial=None,control=True)]
        geometry=cache/'geometry.npz'
        if geometry.exists():
            with np.load(geometry) as g:self.ptr=g['ptr'];self.idx=g['indices'];self.e=g['eigen'];self.q=g['vectors']
        else:
            dist=np.load(args.distance_npy,mmap_mode='r')
            if dist.shape!=(ngray,ngray):raise ValueError('Distance matrix does not match the input grayordinate axis')
            near=[np.flatnonzero(dist[i,:nc]<=args.local_exclusion_mm).astype(np.int32) for i in range(nc)]
            self.ptr=np.r_[0,np.cumsum([len(v) for v in near])].astype(np.int64);self.idx=np.concatenate(near);del near,dist
            adj=csr_matrix((np.ones(len(self.idx)),self.idx,self.ptr),shape=(nc,nc))
            products=(self.a[:,:,None]*self.a[:,None,:]).reshape(nc,400)
            gram=(self.a.T@self.a)[None]-(adj@products).reshape(nc,20,20);del products
            self.e,self.q=np.linalg.eigh(gram);self.e=np.maximum(self.e,0)
            save_npz(geometry,ptr=self.ptr,indices=self.idx,eigen=self.e,vectors=self.q)
        self.core=kernel_library(cache)

    def maps(self,z):
        z=np.ascontiguousarray(z[:self.nc]);corr=local_correlations(self.core,z,self.ptr,self.idx)
        edges=csr_matrix((corr,self.idx,self.ptr),shape=(self.nc,self.nc))
        projection=z@(z.T@self.a)/z.shape[1]-edges@self.a
        if self.args.fc_demean:
            adj=csr_matrix((np.ones(len(self.idx)),self.idx,self.ptr),shape=edges.shape)
            count=self.nc-np.diff(self.ptr)
            total=z@z.sum(0)/z.shape[1]-edges@np.ones(self.nc)
            mean=np.divide(total,count,out=np.zeros(self.nc),where=count>0)
            projection-=mean[:,None]*(self.a.sum(0)[None]-adj@self.a)
        b=np.einsum('vij,vi->vj',self.q,projection);result=[]
        for lam in self.args.lambda_values:
            beta=np.einsum('vij,vj->vi',self.q,b/(self.e+lam));fc=(beta-beta.mean(1,keepdims=True))/np.maximum(beta.std(1,keepdims=True),1e-8)
            for weight in self.args.weight_values:result.append((np.argmax(self.args.fc_weight*fc+weight*self.term,axis=1)+1).astype(np.uint8))
        return np.asarray(result+[self.atlas])

def run(args):
    start=time.time();root=Path(args.outdir).resolve();root.mkdir(parents=True,exist_ok=True)
    try:image,ci,split,pairs,pri,rotations=preflight(args)
    except (ValueError,OSError,KeyError) as e:
        write(root/'selection_failure.json',dict(state='preflight_failed',reason=str(e)))
        from pfm_ridge_cv_report import failure_report
        failure_report(root,str(e));raise
    if args.preflight_only:
        write(root/'preflight.json',dict(state='passed',comparisons=len(pairs),tr=args.actual_tr,tr_source=args.tr_source,partitions=pairs));print('[ridge-cv] preflight passed');return
    ident,bound=binding(args);out=root/'runs'/ident;out.mkdir(parents=True,exist_ok=True);cache=out/'cache';cache.mkdir(exist_ok=True)
    write(out/'binding.json',bound)
    result_path=out/'results.json'
    def progress(stage,**extra):
        record=dict(stage=stage,elapsed_minutes=(time.time()-start)/60,signature=ident,**extra);write(out/'status.json',record);write(root/'status.json',record);print('[ridge-cv]',json.dumps(record),flush=True)
    if result_path.exists() and (out/'completion.json').exists():
        complete=json.loads((out/'completion.json').read_text())
        for relative,sha in complete['outputs'].items():
            if digest(out/relative)!=sha:raise ValueError('Completed CV output changed: '+relative)
        payload=json.loads(result_path.read_text());progress('reuse_completed_grid')
    else:
        nc=len(ci);caxis=image.header.get_axis(1)[ci];progress('geometry',comparisons=len(pairs),grid_cells=len(args.lambda_values)*len(args.weight_values))
        fitter=Fitter(args,pri,nc,image.shape[1],cache);area=surface_areas(caxis,args.left_surf,args.right_surf)
        names=[str(v).strip() for v in np.asarray(pri.NetworkLabels).ravel()[:20]];colors=np.asarray(pri.NetworkColors)[:20,:3]
        table={0:('???',(1,1,1,0))};table.update({i+1:(n,tuple(colors[i])+(1,)) for i,n in enumerate(names)})
        write(out/'partitions.json',dict(tr=args.actual_tr,tr_source=args.tr_source,guard_seconds_each_side=split['gap_seconds_each_side'],input_raw_rows=split['input_raw_rows'],runs=split['runs'],comparisons=pairs))
        data=np.asarray(image.dataobj);hs=[];ds=[];network_dice=[];network_effect=[];supports=[];pair_details=[]
        for pi,pair in enumerate(pairs):
            directory=out/'comparisons'/f'p{pi+1:02d}';directory.mkdir(parents=True,exist_ok=True);scorepath=directory/'scores.npz'
            if scorepath.exists():
                with np.load(scorepath) as r:h=r['fc'];d=r['dice'];nd=r['network_dice'];ne=r['network_fc'];support=r['support'];targets=int(r['targets'])
            else:
                progress('fit_comparison',comparison=pi+1);z={};valid={};labels={}
                for side in ['A','B']:
                    p=pair[side];z[side],valid[side]=normalized(data,p['rows'],p['groups'])
                    if not valid[side][:nc].all():raise ValueError('Invalid cortical time series in '+side)
                    labpath=directory/(side+'_labels.npz')
                    if labpath.exists():labels[side]=np.load(labpath)['labels']
                    else:
                        labels[side]=fitter.maps(z[side]);save_npz(labpath,labels=labels[side])
                        write_labels(directory/(side+'_NetworkAssignments.dlabel.nii'),labels[side],[s['key'] for s in fitter.settings],table,caxis)
                for side in ['A','B']:
                    if not (directory/(side+'_NetworkAssignments.dlabel.nii')).exists():write_labels(directory/(side+'_NetworkAssignments.dlabel.nii'),labels[side],[s['key'] for s in fitter.settings],table,caxis)
                mask=valid['A']&valid['B'];targets=int(mask.sum());directions=[]
                for train,test in [('A','B'),('B','A')]:
                    path=directory/(train+'to'+test+'_profiles.npz')
                    if path.exists():raw=dict(np.load(path))
                    else:
                        progress('prepare_profiles',comparison=pi+1,test=test);metric=prepare(z[test][:nc],z[test][mask])
                        raw=evaluate(metric,labels[train],rotations,lambda n,total:progress('profile_rotations',comparison=pi+1,train=train,test=test,rotation=n,total=total));save_npz(path,**raw);del metric
                    directions.append(raw)
                h,ne,support=aggregate_directions(*directions)
                nd=repeatability(labels['A'],labels['B'],area);d=nd.mean(1)
                save_npz(scorepath,fc=h,dice=d,network_dice=nd,network_fc=ne,support=support,targets=targets)
                del z,valid,labels,directions;gc.collect()
            hs.append(h);ds.append(d);network_dice.append(nd);network_effect.append(ne);supports.append(support);pair_details.append(dict(comparison=pi+1,minutes_A=pair['A']['minutes'],minutes_B=pair['B']['minutes'],whole_brain_targets=targets))
        h=np.mean(hs,axis=0);d=np.mean(ds,axis=0);nd=np.mean(network_dice,axis=0);ne=np.mean(network_effect,axis=0)
        try:selection=choose(h,d,fitter.settings,args.max_fc_loss)
        except ValueError as e:
            failure=dict(state='no_supported_selection',reason=str(e),signature=ident,fc_by_pair=hs,dice_by_pair=ds,rotation_support=supports)
            write(root/'selection_failure.json',failure);write(out/'selection_failure.json',failure)
            from pfm_ridge_cv_report import failure_report
            failure_report(root,str(e));raise
        rows=[]
        for i,s in enumerate(fitter.settings):
            rows.append(dict(s,index=i,fc=h[i],dice=d[i],network_dice=nd[i],network_fc=ne[i],eligible=bool(selection['eligible'][i]),pareto=bool(selection['pareto'][i]),valid_rotations=int(np.min(np.asarray(supports)[:,i])),exclusion_reason='' if selection['eligible'][i] else 'control' if s['control'] else 'no valid paired rotation or insufficient 20-network Dice support in at least one comparison'))
        fullpath=out/'FullData_labels.npz'
        if fullpath.exists():full=np.load(fullpath)['labels']
        else:
            progress('all_data_grid_refit');z,valid=normalized(data,np.arange(len(data)))
            if not valid[:nc].all():raise ValueError('Invalid cortical signal in all-data fit')
            full=fitter.maps(z);save_npz(fullpath,labels=full)
            write_labels(out/'FullData_NetworkAssignments.dlabel.nii',full,[s['key'] for s in fitter.settings],table,caxis);del z
        if not (out/'FullData_NetworkAssignments.dlabel.nii').exists():write_labels(out/'FullData_NetworkAssignments.dlabel.nii',full,[s['key'] for s in fitter.settings],table,caxis)
        del data
        selected=selection['index'];write_labels(out/'RidgePFM_CV_Selected_Cortex.dlabel.nii',full[selected:selected+1],['CV-selected cortex'],table,caxis)
        payload=dict(schema=SCHEMA,subject=args.subject,state='complete',signature=ident,rows=rows,selection=selection,lambda_values=args.lambda_values,spatial_weights=args.weight_values,network_names=names,colors=colors,comparisons=len(pairs),max_comparisons=args.max_comparisons,min_partition_minutes=args.minimum_partition_minutes,rotations_requested=args.n_rotations,rotations_common=args.n_rotations,pairs=pair_details,partitions={s:dict(frames=pairs[0][s]['frames'],minutes=pairs[0][s]['minutes']) for s in ['A','B']},guard_seconds_each_side=split['gap_seconds_each_side'],retained_minutes=image.shape[0]*args.actual_tr/60,TR_seconds=args.actual_tr,TR_source=args.tr_source,source_input=str(Path(args.in_cifti).resolve()),cortex_vertices=nc,targets=pair_details[0]['whole_brain_targets'],render_enabled=bool(args.render),eligible_candidates=int(selection['eligible'].sum()),scope='Cortical CV selection; standard mapper applies the saved pair to requested structures',preprocessing_scope='Fixed run-level preprocessing; conditional temporal CV, not nested upstream preprocessing',fc_by_pair=hs,dice_by_pair=ds,rotation_support=supports,source_settings=dict(lam=None,spatial=None),source_unchanged=True,elapsed_minutes=(time.time()-start)/60)
        for item in bound['inputs']:
            stat=Path(item['path']).stat()
            if [stat.st_size,stat.st_mtime_ns]!=[item['bytes'],item['mtime_ns']]:raise ValueError('Input changed during CV: '+item['path'])
        payload=clean(payload)
        write(result_path,payload)
        with (out/'grid_scores.csv').open('w') as f:
            columns=['key','lam','spatial','control','fc','dice','eligible','pareto','valid_rotations','exclusion_reason'];w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows([{k:clean(r[k]) for k in columns} for r in rows])
        with (out/'network_scores.csv').open('w') as f:
            w=csv.writer(f);w.writerow(['key','network','mean_Dice','mean_heldout_profile_advantage'])
            for row in rows:
                for j,name in enumerate(names):w.writerow([row['key'],name,clean(nd[row['index'],j]),clean(ne[row['index'],j])])
        output_files=['results.json','partitions.json','grid_scores.csv','network_scores.csv','FullData_labels.npz','FullData_NetworkAssignments.dlabel.nii','RidgePFM_CV_Selected_Cortex.dlabel.nii']
        for pi in range(len(pairs)):
            output_files += [f'comparisons/p{pi+1:02d}/{s}_{name}' for s in ['A','B'] for name in ['labels.npz','NetworkAssignments.dlabel.nii']]
        write(out/'completion.json',dict(state='complete',outputs={p:digest(out/p) for p in output_files}))
    from pfm_ridge_cv_report import render_report
    render_report(out,args,payload)
    ident=payload['signature'];choice=payload['selection'];winner=payload['rows'][choice['index']]
    selected=dict(schema=SCHEMA,lambda_value=winner['lam'],spatial_weight=winner['spatial'],max_fc_loss=args.max_fc_loss,FC_advantage_retained=winner['fc']/choice['maximum_fc'],network_Dice=winner['dice'],signature=ident,run_directory=str(out),report=str(root/'report.html'),criterion='Maximize mean network Dice on Pareto frontier retaining declared fraction of maximum positive held-out FC-profile advantage',comparisons=payload['comparisons'])
    write(out/'selected_parameters.json',selected);write(root/'selected_parameters.json',selected)
    page=(out/'report.html').read_text().replace('<head>','<head><base href="runs/'+ident+'/">',1);(root/'report.html').write_text(page)
    if (root/'selection_failure.json').exists():(root/'selection_failure.json').unlink()
    if not args.keep_cache and cache.exists():shutil.rmtree(cache)
    progress('complete',lambda_value=winner['lam'],spatial_weight=winner['spatial'])

if __name__=='__main__':
    try:run(parser().parse_args())
    except (ValueError,OSError,KeyError) as e:
        print('[ridge-cv] ERROR:',str(e),file=sys.stderr);sys.exit(2)
