"""Independent BFS/SciPy and exact invariance checks for the scientific operator."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
R=Path(__file__).resolve().parents[2];sys.path.insert(0,str(R))
from research.fractalgo.core import *

def main():
    rng=np.random.default_rng(2026111901);checks=[]
    for n in [1,2,7,19]:
        edges=path_edges(n);q=rng.normal(size=n);pi=np.exp(rng.normal(size=n));a=excursion(q,pi,edges,.8)
        b=excursion(q+100,pi*7,edges,.8)
        assert a['node_groups']==b['node_groups'] and np.allclose(a['posterior'],b['posterior'])
        perm=rng.permutation(n);inv=np.argsort(perm)
        c=excursion(q[perm],pi[perm],[(int(inv[i]),int(inv[j])) for i,j in edges],.8)
        assert sorted(map(len,a['node_groups']))==sorted(map(len,c['node_groups']))
        assert np.isclose(sum(a['component_masses']),a['selected_mass'])
        checks.append(f'common_shift_prior_scaling_relabel_mass_n{n}')
    # Independent graph implementation for arbitrary, disconnected graphs.
    for trial in range(40):
        n=15;mask=rng.random(n)>.5;adj=np.triu(rng.random((n,n))<.15,1);adj=adj|adj.T
        edges=list(zip(*np.where(np.triu(adj,1))));g=graph_components(mask,edges)
        sub=csr_matrix(adj[np.ix_(mask,mask)])
        assert len(g)==connected_components(sub,directed=False,return_labels=False)
    checks.append('40_independent_scipy_graph_checks')
    for trial in range(40):
        n=20;connected=rng.random(n-1)>.2;v=rng.normal(size=n);mass=rng.random(n);mass/=mass.sum()
        out=path_summary(v,.4,connected,mass);g=graph_components(v>.4,[(i,i+1) for i in range(n-1) if connected[i]])
        assert out[0][0]==len(g) and out[1][0]==max(map(len,g),default=0)
        assert np.isclose(out[2][0],max([mass[z].sum() for z in g],default=0))
    checks.append('40_path_vs_bfs_gap_and_mass_checks')
    n=80;va=np.r_[np.ones(40),-np.ones(40)];vb=np.tile([1.,-1.],40)
    for z in [-5,-1,.1,1,5]:
        a=excursion(z*va,np.ones(n),path_edges(n));b=excursion(z*vb,np.ones(n),path_edges(n))
        assert a['components']==1 and b['components']==40
        assert np.allclose(np.sort(a['posterior']),np.sort(b['posterior']))
    checks.append('exact_gaussian_witness_5_realizations')
    assert excursion(np.zeros(3),[1,2,3],path_edges(3))['components']==0
    for prior in [[0,1],[1,-1],[1,float('nan')]]:
        try:posterior([0,1],prior)
        except ValueError:pass
        else:raise AssertionError('invalid prior accepted')
    checks.append('empty_excursion_and_invalid_prior_rejection')
    out=R/'research/fractalgo/results';out.mkdir(exist_ok=True)
    (out/'tests.json').write_text(json.dumps({'status':'passed','checks':checks,'seed':2026111901},indent=2)+'\n')
    print('PASS',len(checks),'check groups, including 80 independent graph comparisons')
if __name__=='__main__':main()
