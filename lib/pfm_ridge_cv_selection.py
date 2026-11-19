"""Explicit FC-loss budget, with no data-dependent metric rescaling."""
import numpy as np

def choose(fc,dice,settings,loss=.1):
    if not 0<=loss<=1: raise ValueError('FC loss must be between 0 and 1')
    h=np.asarray(fc,float);r=np.asarray(dice,float)
    eligible=np.isfinite(h)&np.isfinite(r)&np.array([not s.get('control',False) for s in settings])
    ix=np.flatnonzero(eligible)
    if not len(ix):raise ValueError('No evaluable grid cells')
    best=int(ix[np.argmax(h[ix])]);maximum=h[best]
    if maximum<=0:raise ValueError('No positive held-out FC advantage; no automatic CV selection')
    threshold=(1-loss)*maximum
    allowed=eligible&(h>=threshold-1e-12)
    pareto=eligible.copy()
    for i in ix:pareto[i]=not np.any(eligible&(h>=h[i])&(r>=r[i])&((h>h[i])|(r>r[i])))
    pool=np.flatnonzero(allowed&pareto)
    selected=int(sorted(pool,key=lambda i:(-r[i],-h[i],i))[0])
    return dict(index=selected,fc_best=best,threshold=float(threshold),maximum_fc=float(maximum),loss_budget=loss,eligible=eligible,allowed=allowed,pareto=pareto)
