"""Fixed-graph, prior-relative posterior excursions. No physical/fractal claims."""
import numpy as np
from scipy.special import logsumexp

def posterior(costs, prior):
    q=np.asarray(costs,dtype=float); pi=np.asarray(prior,dtype=float)
    if q.ndim!=1 or pi.shape!=q.shape or not len(q): raise ValueError('shape')
    if not np.isfinite(q).all() or not np.isfinite(pi).all() or np.any(pi<=0): raise ValueError('finite costs and positive prior required')
    pi=pi/pi.sum()
    # Subtract a common offset before normalizing; common modes carry no contrast information.
    lq=-.5*(q-q.min()); lr=lq-logsumexp(np.log(pi)+lq)
    return np.exp(np.log(pi)+lr),lr

def graph_components(mask, edges):
    mask=np.asarray(mask,dtype=bool)
    if mask.ndim!=1: raise ValueError('mask shape')
    adj=[[] for _ in mask]
    for i,j in edges:
        if int(i)!=i or int(j)!=j or min(i,j)<0 or max(i,j)>=len(mask) or i==j: raise ValueError('invalid edge')
        i,j=int(i),int(j);adj[i].append(j);adj[j].append(i)
    seen=set();groups=[]
    for i in np.flatnonzero(mask):
        if int(i) in seen:continue
        stack=[int(i)];seen.add(int(i));g=[]
        while stack:
            j=stack.pop();g.append(j)
            for k in adj[j]:
                if mask[k] and k not in seen:seen.add(k);stack.append(k)
        groups.append(sorted(g))
    return groups

def excursion(costs,prior,edges,threshold=1.):
    if not np.isfinite(threshold) or threshold<=0:raise ValueError('positive threshold required')
    p,lr=posterior(costs,prior);mask=lr>np.log(threshold)
    groups=graph_components(mask,edges)
    return dict(components=len(groups),node_groups=groups,
                component_masses=[float(p[g].sum()) for g in groups],
                selected_mass=float(p[mask].sum()),posterior=p,log_ratio=lr)

def path_summary(values, threshold, connected=None, mass=None):
    """Vectorized component counts, followed by explicit component-mass sums.
    `connected[j]` declares whether nodes j and j+1 are adjacent.
    Missing edges never bridge gaps. Empty excursions return zero components.
    """
    v=np.atleast_2d(values);m=v>threshold;n=v.shape[1]
    c=np.ones(n-1,dtype=bool) if connected is None else np.asarray(connected,dtype=bool)
    if c.shape!=(n-1,):raise ValueError('edge shape')
    start=m & np.c_[np.ones(len(v),dtype=bool),~(m[:,:-1]&c)]
    end=m & np.c_[~(m[:,1:]&c),np.ones(len(v),dtype=bool)]
    count=start.sum(axis=1);longest=np.zeros(len(v),dtype=int);largest=np.zeros(len(v))
    w=np.broadcast_to(np.ones(n)/n if mass is None else mass,v.shape)
    for k in range(len(v)):
        a=np.flatnonzero(start[k]);b=np.flatnonzero(end[k])+1
        if len(a):
            longest[k]=max(b-a);cw=np.r_[0,np.cumsum(w[k])];largest[k]=max(cw[b]-cw[a])
    return count,longest,largest

def path_edges(n):return [(i,i+1) for i in range(n-1)]

def cell_weights(x):
    x=np.asarray(x,dtype=float)
    if x.ndim!=1 or len(x)<2 or not np.isfinite(x).all() or np.any(np.diff(x)<=0):raise ValueError('strictly increasing grid required')
    w=np.empty_like(x);w[0]=(x[1]-x[0])/2;w[-1]=(x[-1]-x[-2])/2;w[1:-1]=(x[2:]-x[:-2])/2
    return w
