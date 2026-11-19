"""Same GP marginal bands, different joint laws: quantile uncertainty and mesh.

Controlled experiment, NOT the undisclosed implementation or data of B-INPOP.
Mean is set to the exact conditional profile to isolate propagation of uncertainty.
"""
from pathlib import Path
import json,csv
import numpy as np
from scipy.linalg import cho_factor,cho_solve,eigh
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq
R=Path(__file__).resolve().parents[1];OUT=R/'results'
SEED=20261117;NSIM=4096;ELL=.5;SIGMA_Q=3.

def quantiles_and_mix(x,Q):
    logf=-.5*(Q-Q.min(axis=1)[:,None]);f=np.exp(logf);d=np.diff(logf,axis=1)
    # Integrate exp(linearly interpolated Q) exactly on every mesh interval.
    # Trapezoidal integration of exp(Q) is not equivalent for independent draws.
    small=abs(d)<1e-6
    rel=np.ones_like(d)
    rel[~small]=np.expm1(d[~small])/d[~small]
    rel[small]=1+d[small]/2+d[small]**2/6+d[small]**3/24
    seg=np.diff(x)[None,:]*f[:,:-1]*rel
    cum=np.cumsum(seg,axis=1)
    norm=cum[:,-1];cdf=np.c_[np.zeros(len(Q)),cum/norm[:,None]]
    idx=np.argmax(cdf>=.997,axis=1);ar=np.arange(len(Q))
    low=cdf[ar,idx-1];hi=cdf[ar,idx]
    fraction=(.997-low)/(hi-low);di=d[ar,idx-1]
    t=fraction.copy();ok=abs(di)>1e-8
    t[ok]=np.log1p(fraction[ok]*np.expm1(di[ok]))/di[ok]
    q=x[idx-1]+t*(x[idx]-x[idx-1])
    mix=(f/norm[:,None]).mean(axis=0);mixcdf=cdf.mean(axis=0)
    k=np.searchsorted(mixcdf,.997)-1;di=d[:,k];ok=abs(di)>1e-8
    def target(t):
        fraction=np.full(len(di),t)
        fraction[ok]=np.expm1(di[ok]*t)/np.expm1(di[ok])
        return np.mean(cdf[:,k]+seg[:,k]/norm*fraction)-.997
    tm=brentq(target,0,1)
    return q,float(x[k]+tm*(x[k+1]-x[k])),mix

def main():
    emulator=json.loads((OUT/'gp_emulator.json').read_text());zhat=emulator['z_hat']
    x=np.linspace(0,3,801);nodes=np.array([0.,.6,1.5,3.])
    kernel=lambda a,b:SIGMA_Q**2*np.exp(-.5*((a[:,None]-b[None,:])/ELL)**2)
    cross=kernel(x,nodes);chol=cho_factor(kernel(nodes,nodes),lower=True)
    C=kernel(x,x)-cross@cho_solve(chol,cross.T);C=(C+C.T)/2
    eig,U=eigh(C);assert eig[0]>-1e-10
    keep=eig>1e-10;L=U[:,keep]*np.sqrt(eig[keep])
    assert np.max(abs(L@L.T-C))<1e-10
    rng=np.random.default_rng(SEED);coherent=rng.standard_normal((NSIM,L.shape[1]))@L.T
    marginal=np.sqrt(np.maximum(np.diag(C),0));independent=rng.standard_normal((NSIM,len(x)))*marginal
    mean=x**4-2*zhat*x*x;rows=[];curves={};bands={}
    for count in [51,101,201,401,801]:
        ids=np.linspace(0,800,count,dtype=int);xx=x[ids]
        for law,draw in [('joint_GP',coherent),('independent_marginals',independent)]:
            q,mix,density=quantiles_and_mix(xx,mean[ids]+draw[:,ids])
            row=dict(law=law,grid_points=count,quantile_997_mean=float(q.mean()),
                quantile_997_std=float(q.std(ddof=1)),q997_distribution_025=float(np.quantile(q,.025)),
                q997_distribution_975=float(np.quantile(q,.975)),pooled_posterior_quantile_997=mix,
                std_first_half=float(q[:NSIM//2].std(ddof=1)),std_second_half=float(q[NSIM//2:].std(ddof=1)))
            rows.append(row)
            if count==801:
                curves[law+'__pooled_density']=density;curves[law+'__q997_draws']=q
    # Marginalize unnormalized likelihood, an explicitly different target from
    # averaging normalized posteriors. Both choices must be named, not conflated.
    avgL=np.exp(-mean/2+np.diag(C)/8);cdf=cumulative_trapezoid(avgL,x,initial=0);cdf/=cdf[-1]
    avgLq=float(np.interp(.997,cdf,x))
    # The increment law is the test that one-point distributions cannot specify.
    k=np.argmin(abs(x-.3));increment=[]
    for steps in [1,2,4,8]:
        i=k;j=k+steps
        joint=float(C[i,i]+C[j,j]-2*C[i,j]);ind=float(C[i,i]+C[j,j])
        increment.append(dict(delta_x=float(x[j]-x[i]),variance_increment_joint=joint,
            variance_increment_independent=ind,variance_ratio=ind/joint))
    result=dict(assertion_tested='Joint conditional GP versus independent draws with identical pointwise variances',
        scope='Controlled Gaussian-function uncertainty law; exact mean profile; not a replication of published INPOP data',
        integration='Exact integration of exp(-Q/2) on linearly interpolated Q segments; exact per-curve quantile inversion; pooled quantile by scalar root',
        seed=SEED,simulations=NSIM,ell=ELL,sigma_delta_Q=SIGMA_Q,conditioning_nodes=nodes.tolist(),
        max_variance_delta_Q=float(np.diag(C).max()),eigen_rank=int(keep.sum()),
        exact_mean_quantile_997=emulator['results'][0]['exact_quantile_997'],
        quantile_997_of_averaged_unnormalized_likelihood=avgLq,
        increments=increment,results=rows)
    (OUT/'gp_realizations.json').write_text(json.dumps(result,indent=2)+'\n')
    with (OUT/'gp_realizations.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    np.savez_compressed(OUT/'gp_realizations_curves.npz',x=x,mean=mean,pointwise_sd=marginal,
        joint_draw_examples=coherent[:5],independent_draw_examples=independent[:5],**curves)
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
