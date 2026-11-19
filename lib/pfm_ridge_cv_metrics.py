"""Exact signed mean whole-brain FC-profile correlation. No PCA/SVD/eigen solve."""
import numpy as np
from scipy.sparse import csr_matrix

def prepare(seeds,targets):
    x=np.asarray(seeds,float);y=np.asarray(targets,float)
    mean=y.mean(0);metric=y.T@y/len(y)-np.outer(mean,mean)
    metric=(metric+metric.T)/2
    norm2=np.einsum('ij,ij->i',x@metric,x)
    if not np.isfinite(norm2).all() or np.any(norm2<=1e-12):raise ValueError('Unusable FC-profile spatial variance')
    return x/np.sqrt(norm2[:,None]),metric

def mean_grid(prepared,labels):
    w,metric=prepared;labels=np.asarray(labels);m,n=labels.shape;k=20
    if n!=len(w):raise ValueError('Cortical axis mismatch')
    valid=labels>0
    rows=(np.arange(m)[:,None]*k+labels.astype(np.int64)-1)[valid]
    cols=np.broadcast_to(np.arange(n),labels.shape)[valid]
    membership=csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(m*k,n))
    counts=np.asarray(membership.sum(1)).ravel()
    sums=membership@w
    squared=np.einsum('ij,ij->i',sums@metric,sums)
    average=np.divide(squared-counts,counts*(counts-1),out=np.full(len(counts),np.nan),where=counts>1).reshape(m,k)
    if np.nanmax(abs(average))>1+1e-7:raise ValueError('Profile correlation outside valid range')
    counts=counts.reshape(m,k)
    weighted=np.sum(np.where(counts==0,0,average*counts),axis=1)/counts.sum(1)
    return weighted,average,counts

def evaluate(prepared,labels,rotations,callback=None):
    n=labels.shape[1];real=[];null=[];net_real=[];net_null=[]
    for i,rotation in enumerate(rotations):
        valid=(rotation>0)&(rotation<=n)
        original=labels.copy();original[:,~valid]=0
        moved=np.zeros_like(labels);moved[:,valid]=labels[:,rotation[valid]-1]
        a=mean_grid(prepared,original);b=mean_grid(prepared,moved)
        real.append(a[0]);null.append(b[0]);net_real.append(a[1]);net_null.append(b[1])
        if callback and (i==0 or (i+1)%8==0):callback(i+1,len(rotations))
    real=np.asarray(real);null=np.asarray(null)
    return dict(observed=real,null=null,effect=real-null,network_observed=np.asarray(net_real),network_null=np.asarray(net_null))

def aggregate_directions(forward,reverse):
    """Average supported rotations, pairing A→B and B→A for each draw.

    A positive-size singleton makes that draw undefined. Drop the draw for
    that cell in both directions; do not replace it with zero or discard
    the remaining draws. Each temporal comparison still has equal weight.
    """
    effect=(forward['effect']+reverse['effect'])/2
    valid=np.isfinite(effect)
    support=valid.sum(axis=0)
    fc=np.divide(np.where(valid,effect,0).sum(axis=0),support,
                 out=np.full(effect.shape[1],np.nan),where=support>0)
    network=((forward['network_observed']-forward['network_null'])+
             (reverse['network_observed']-reverse['network_null']))/2
    network_valid=valid[:,:,None]&np.isfinite(network)
    counts=network_valid.sum(axis=0)
    network_fc=np.divide(np.where(network_valid,network,0).sum(axis=0),counts,
                         out=np.full(network.shape[1:],np.nan),where=counts>0)
    return fc,network_fc,support

def repeatability(a,b,area):
    values=[]
    for aa,bb in zip(a,b):
        ca=np.bincount(aa,weights=area,minlength=21)[1:21]
        cb=np.bincount(bb,weights=area,minlength=21)[1:21]
        same=aa==bb;inter=np.bincount(aa[same],weights=area[same],minlength=21)[1:21]
        values.append(np.divide(2*inter,ca+cb,out=np.full(20,np.nan),where=ca+cb>0))
    return np.asarray(values)

def select(h,d,settings,tolerance=.005):
    h=np.atleast_2d(h);d=np.atleast_2d(d);mh=h.mean(0);md=d.mean(0)
    eligible=np.isfinite(h).all(0)&np.isfinite(d).all(0)&np.asarray([not s['control'] for s in settings])
    indexes=np.flatnonzero(eligible)
    if not len(indexes):raise ValueError('No eligible parameter candidates')
    best=int(indexes[np.argmax(mh[indexes])]);se=np.std(h-h[:,best,None],axis=0,ddof=1)/np.sqrt(len(h)) if len(h)>1 else np.zeros(len(mh))
    allowed=eligible&((mh[best]-mh)<=np.maximum(tolerance,se)+1e-12)
    pareto=eligible.copy()
    for i in indexes:
        pareto[i]=not np.any(eligible&(mh>=mh[i])&(md>=md[i])&((mh>mh[i])|(md>md[i])))
    pool=np.flatnonzero(allowed&pareto)
    chosen=int(sorted(pool,key=lambda i:(-md[i],-mh[i],i))[0])
    return dict(chosen=chosen,fc_best=best,eligible=eligible,pareto=pareto,allowed=allowed,
                mean_fc=mh,mean_dice=md,paired_se_to_fc_best=se,tolerance=tolerance,
                uncertainty_rule='paired SE across disjoint A/B comparisons when >=2; practical floor; heuristic, not a formal equivalence test')
