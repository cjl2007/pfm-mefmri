#!/usr/bin/env python3
"""Reuse a verified numerical grid when only the FC-loss preference changes.

The numerical engine and all input/resource/code bindings are unchanged. Each
policy receives its own immutable run directory and selected cortical map.
"""
from pathlib import Path
import copy,json,os,shutil,sys,time
import pfm_ridge_cv as engine  # sets thread limits before NumPy is imported
import numpy as np
import nibabel as nib
from pfm_ridge_cv_selection import choose
from pfm_ridge_cv_data import write_labels

def equivalent_grid(a,b):
    a=copy.deepcopy(a);b=copy.deepcopy(b)
    for value in [a,b]:
        value['settings'].pop('max_fc_loss',None)
        # Canonical numeric vectors remain part of the binding.
        value['settings'].pop('lambdas',None);value['settings'].pop('spatial_weights',None)
    return a==b

def clone_grid(source,destination,bound,signature,loss):
    source=Path(source);destination=Path(destination)
    if destination.exists():return False
    record=json.loads((source/'completion.json').read_text())
    for relative,sha in record['outputs'].items():
        if engine.digest(source/relative)!=sha:raise ValueError('Cached numerical grid changed: '+str(source/relative))
    payload=json.loads((source/'results.json').read_text());rows=payload['rows']
    h=[np.nan if r['fc'] is None else r['fc'] for r in rows];d=[np.nan if r['dice'] is None else r['dice'] for r in rows]
    selection=choose(h,d,rows,loss)
    temp=destination.with_name(destination.name+'.reuse.tmp')
    if temp.exists():shutil.rmtree(temp)
    temp.mkdir(parents=True)
    try:
        # Raw immutable binary grids are shared; text/provenance is copied.
        skip={'results.json','completion.json','binding.json','selected_parameters.json','RidgePFM_CV_Selected_Cortex.dlabel.nii','report.html','status.json','cache','cached_grid_source.json'}
        for p in source.rglob('*'):
            relative=p.relative_to(source)
            if not p.is_file() or relative.parts[0] in skip:continue
            dest=temp/relative;dest.parent.mkdir(parents=True,exist_ok=True)
            if p.suffix=='.npz' or p.name.endswith('.nii'):
                try:os.link(p,dest)
                except OSError:shutil.copy2(p,dest)
            else:shutil.copy2(p,dest)
        payload.update(signature=signature,selection=engine.clean(selection),grid_reused_from_signature=payload['signature'])
        full=nib.load(temp/'FullData_NetworkAssignments.dlabel.nii');labels=np.load(temp/'FullData_labels.npz')['labels']
        write_labels(temp/'RidgePFM_CV_Selected_Cortex.dlabel.nii',labels[selection['index']:selection['index']+1],['CV-selected cortex'],full.header.get_axis(0).label[0],full.header.get_axis(1))
        engine.write(temp/'results.json',payload);engine.write(temp/'binding.json',bound)
        engine.write(temp/'cached_grid_source.json',dict(source_run=str(source.resolve()),source_signature=payload['grid_reused_from_signature'],requested_max_fc_loss=loss,code_inputs_and_numerical_settings_match=True,dispatcher_sha256=engine.digest(__file__)))
        engine.write(temp/'completion.json',dict(state='complete',outputs={p:engine.digest(temp/p) for p in record['outputs']}))
        temp.replace(destination)
    except Exception:
        shutil.rmtree(temp);raise
    return True

def run(args):
    # Delegate preflight failures to the engine's ordinary failure reporting.
    try:engine.preflight(args);signature,bound=engine.binding(args)
    except (ValueError,OSError,KeyError):return engine.run(args)
    root=Path(args.outdir).resolve();destination=root/'runs'/signature
    reused=False;begin=time.time()
    if not args.preflight_only and not destination.exists() and (root/'runs').exists():
        current=json.loads((root/'selected_parameters.json').read_text()).get('run_directory') if (root/'selected_parameters.json').exists() else None
        sources=sorted((root/'runs').iterdir(),key=lambda p:(str(p.resolve())!=current,p.name))
        for source in sources:
            if len(source.name)!=20 or any(c not in '0123456789abcdef' for c in source.name):continue
            if not (source/'binding.json').exists() or not (source/'completion.json').exists():continue
            old=json.loads((source/'binding.json').read_text())
            if equivalent_grid(old,bound):
                reused=clone_grid(source,destination,bound,signature,args.max_fc_loss)
                if reused:print('[ridge-cv] reuse verified grid; selecting a new FC-loss preference without FC/null recomputation',flush=True);break
    engine.run(args)
    selected_path=root/'selected_parameters.json'
    if args.preflight_only or not selected_path.exists():return
    selected=json.loads(selected_path.read_text());out=Path(selected['run_directory'])
    if (out/'cached_grid_source.json').exists():
        lineage=json.loads((out/'cached_grid_source.json').read_text());selected['grid_reused_from_signature']=lineage['source_signature']
        selected['grid_reselection_seconds']=time.time()-begin
        engine.write(selected_path,selected);engine.write(out/'selected_parameters.json',selected)

if __name__=='__main__':
    try:run(engine.parser().parse_args())
    except (ValueError,OSError,KeyError) as e:
        print('[ridge-cv] ERROR:',str(e),file=sys.stderr);sys.exit(2)
