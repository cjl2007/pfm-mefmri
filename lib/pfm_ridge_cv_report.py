"""Offline full-grid report with flat and inflated native cortical views."""
from pathlib import Path
import base64,json,html
import numpy as np

def enc(x):return base64.b64encode(np.ascontiguousarray(x).tobytes()).decode()

def failure_report(root,reason):
    root=Path(root);(root/'report.html').write_text('<!doctype html><html><meta charset="utf-8"><title>RidgePFM CV could not select parameters</title><body style="font:16px system-ui;margin:40px"><h1>No supported CV selection</h1><p>'+html.escape(reason)+'</p><p>No fixed-parameter fallback was substituted. See selection_failure.json for details.</p></body></html>')

def render_report(out,args,payload):
    out=Path(out);d=json.loads(json.dumps(payload,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item()))
    if args.render:
        from pfm_ridge_cv_workbench import capture_grid
        d['workbench']=capture_grid(out,args,d['rows'])
    else:
        d['workbench']=None
    a=np.load(out/'comparisons/p01/A_labels.npz')['labels'];b=np.load(out/'comparisons/p01/B_labels.npz')['labels']
    for row in d['rows']:
        row['vertices_A']=np.bincount(a[row['index']],minlength=21)[1:21].tolist();row['vertices_B']=np.bincount(b[row['index']],minlength=21)[1:21].tolist()
    d['axis_limits']={}
    for name in ['fc','dice']:
        vals=[r[name] for r in d['rows'] if not r['control'] and r[name] is not None]
        lo=min(vals);hi=max(vals);pad=max((hi-lo)*.09,.015);d['axis_limits'][name]=[max(-1,lo-pad),min(1,hi+pad)]
    text=Path(__file__).with_name('pfm_ridge_cv_report.html').read_text()
    text=text.replace('__SUBJECT__',html.escape(args.subject)).replace('__DATA__',json.dumps(d,separators=(',',':'),allow_nan=False).replace('<','\\u003c'))
    (out/'report.html').write_text(text)
