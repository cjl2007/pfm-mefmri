#!/usr/bin/env python3
"""Reuse unchanged PFM preprocessing for CV, preserving its input identity."""
from pathlib import Path
import argparse,hashlib,json,shutil

def record(path):
    p=Path(path).resolve();s=p.stat();return dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns)

def context(args):
    paths=[args.source,args.left_surf,args.right_surf]
    if args.scan_idx and Path(args.scan_idx).is_file():paths.append(args.scan_idx)
    wb=shutil.which('wb_command')
    if wb:paths.append(wb)
    return dict(inputs=[record(p) for p in paths],regress_enable=args.regress_enable,distance_mm=args.distance_mm,sigma_mm=args.sigma_mm,regression_code=hashlib.sha256(Path(args.regression_code).read_bytes()).hexdigest(),cache_code=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),pipeline_code=hashlib.sha256(Path(args.pipeline_code).read_bytes()).hexdigest())

def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)

def check(args):
    cache=Path(args.cache);requested=context(args)
    if cache.exists():
        try:
            saved=json.loads(cache.read_text());output=Path(saved['output']['path'])
            if saved['context']==requested and output.is_file() and record(output)==saved['output']:return str(output)
        except (OSError,ValueError,KeyError):pass
    write(str(cache)+'.pending',requested);return ''

def save(args):
    cache=Path(args.cache);planned=json.loads(Path(str(cache)+'.pending').read_text())
    if planned!=context(args):raise ValueError('PFM input or preprocessing context changed during preparation')
    write(cache,dict(context=planned,output=record(args.prepared_output)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['check','record'],required=True)
    for key in ['cache','source','left-surf','right-surf','regression-code','pipeline-code']:p.add_argument('--'+key,required=True)
    p.add_argument('--scan-idx',default='');p.add_argument('--regress-enable',type=int,required=True);p.add_argument('--distance-mm',type=float,required=True);p.add_argument('--sigma-mm',type=float,required=True);p.add_argument('--prepared-output')
    a=p.parse_args()
    if a.mode=='check':print(check(a))
    else:save(a)
