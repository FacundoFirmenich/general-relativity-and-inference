"""Proper hierarchical marginal likelihood on the inherited Cassini design.

This is a conditional optical-amplitude audit, NOT a new PPN or graviton bound.
The empirical weights and optical template are conditioned upon. Priors are
explicit sensitivity scenarios, not certified instrument calibrations.
"""
from pathlib import Path
import json, sys, time
import numpy as np
from scipy.linalg import eigh, cho_factor, cho_solve
from scipy.special import logsumexp, ndtri, ndtr
from scipy.optimize import minimize, brentq
from scipy.stats import norm, halfnorm, qmc, rankdata
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analyses/gp/code'))
from gp_covarianza_cassini import intkernel
BASE=1e-13; AUNIT=1e-3

class ConditionalModel:
    def __init__(self):
        p=ROOT/'analyses/gp/sources/cassini'
        z=np.load(p/'dss25_diseno_erfa_constant.npz')
        t=np.load(p/'contraste_robusto_pesos.npz')
        self.pid=z['pass_id']; self.w=z['weights']; X=z['design']
        # This version explicitly uses the hybrid optical residual from the base audit.
        yy=self.w*(z['residual']+t['constant__delta_SPI'])
        xx=X*self.w[:,None]
        Q,s,_=np.linalg.svd(xx/np.linalg.norm(xx,axis=0),full_matrices=False)
        assert s[-1]>1e-8
        tt=self.w*t['constant__template_per_A']*AUNIT
        self.blocks=[]; rows=[]; eigen=[]
        for pno in np.unique(self.pid):
            ids=np.flatnonzero(self.pid==pno)
            dt=z['epoch'][ids,None]-z['epoch'][None,ids]
            K=intkernel(dt,7200.)*(BASE*self.w[ids,None])*(BASE*self.w[None,ids])
            e,U=eigh(K); assert e[0]>-1e-7*max(1,e[-1])
            eigen.extend(np.maximum(e,0)); rows.append(U.T@np.c_[yy[ids],tt[ids],Q[ids]])
            self.blocks.append((ids,U))
        self.ev=np.array(eigen); Z=np.vstack(rows)
        self.y,self.t,self.Q=Z[:,0],Z[:,1],Z[:,2:]
        self.n,self.p=self.Q.shape
        self.observed=yy

    def stats(self,x,B=100.):
        u,v=x
        if not np.isfinite(u+v) or abs(u)>30 or abs(v)>30: return (-np.inf,0.,0.)
        d=np.exp(2*u)+np.exp(2*v)*self.ev; inv=1/d
        DQ=self.Q*inv[:,None]
        H=np.eye(self.p)/B**2+self.Q.T@DQ
        cf=cho_factor(H,lower=True,check_finite=False)
        b=self.Q.T@(inv[:,None]*np.c_[self.y,self.t])
        sol=cho_solve(cf,b,check_finite=False)
        vals=np.c_[self.y,self.t].T@(inv[:,None]*np.c_[self.y,self.t])-b.T@sol
        ld=np.log(d).sum()+2*self.p*np.log(B)+2*np.log(np.diag(cf[0])).sum()
        return -.5*(self.n*np.log(2*np.pi)+ld+vals[0,0]),max(vals[1,1],0.),vals[0,1]

    def loglike(self,x,sa=3.,B=100.):
        L,j,h=self.stats(x,B)
        return L if sa==0 else L-.5*np.log1p(sa*sa*j)+.5*sa*sa*h*h/(1+sa*sa*j)

    def inverse(self,x,z,sa=3.,B=100.):
        d=np.exp(2*x[0])+np.exp(2*x[1])*self.ev
        F=np.c_[B*self.Q,sa*self.t] if sa else B*self.Q
        H=np.eye(F.shape[1])+F.T@(F/d[:,None])
        return z/d-(F/d[:,None])@np.linalg.solve(H,F.T@(z/d))

    def unrotate(self,z):
        out=np.empty_like(z); start=0
        for ids,U in self.blocks:
            out[ids]=U@z[start:start+len(ids)];start+=len(ids)
        return out

def logprior(x,ts=1.):
    u,v=x
    if abs(u)>30 or abs(v)>30:return -np.inf
    return -.5*(u/.35)**2-np.log(.35*np.sqrt(2*np.pi))+.5*np.log(2/np.pi)+v-np.log(ts)-.5*(np.exp(v)/ts)**2

def transform(q,ts=1.):
    q=np.clip(q,1e-14,1-1e-14)
    return np.array([.35*ndtri(q[0]),np.log(ts*ndtri((q[1]+1)/2))])

def tensor_reference(m,n=141,sa=3.,ts=1.,B=100.,ubounds=(-1.8,2.5),vbounds=(-11.,3.)):
    # Bounds are deliberately broader than the posterior; convergence is recorded.
    u=np.linspace(*ubounds,n);v=np.linspace(*vbounds,2*n-1)
    x=np.array(np.meshgrid(u,v,indexing='ij')).reshape(2,-1).T
    st=np.array([m.stats(p,B) for p in x]); l0,j,h=st.T
    ll=l0 if sa==0 else l0-.5*np.log1p(sa*sa*j)+.5*sa*sa*h*h/(1+sa*sa*j)
    lp=np.array([logprior(p,ts) for p in x])
    weights=np.ones((len(u),len(v))); weights[[0,-1],:]*=.5;weights[:,[0,-1]]*=.5
    lw=ll+lp+np.log(weights.ravel()*(u[1]-u[0])*(v[1]-v[0]))
    z=logsumexp(lw);pw=np.exp(lw-z)
    return dict(logz=float(z),x=x,weights=pw,stats=st,
        mean=sa*sa*h/(1+sa*sa*j),variance=sa*sa/(1+sa*sa*j))

def summary(ref):
    w=ref['weights']; means=ref['mean']; s=np.sqrt(ref['variance'])
    quant=[]
    for p in [.025,.5,.975]:
        quant.append(AUNIT*brentq(lambda a: w@ndtr((a-means)/s)-p,-100,100))
    return dict(logz=ref['logz'],A_quantiles_025_50_975=quant,
        s_mean=float(w@np.exp(ref['x'][:,0])),tau_mean=float(BASE*w@np.exp(ref['x'][:,1])))

def rank_rhat(chains):
    n=chains.shape[1]//2
    c=np.concatenate([chains[:,:n],chains[:,-n:]],axis=0)
    out=[]
    for j in range(c.shape[-1]):
        z=ndtri((rankdata(c[:,:,j].ravel())-.375)/(c.shape[0]*n+.25)).reshape(c.shape[:2])
        def r(z):
            W=np.var(z,axis=1,ddof=1).mean();B=n*np.var(z.mean(axis=1),ddof=1)
            return np.sqrt(((n-1)*W/n+B/n)/W)
        folded=abs(c[:,:,j]-np.median(c[:,:,j]))
        f=ndtri((rankdata(folded.ravel())-.375)/(folded.size+.25)).reshape(folded.shape)
        out.append(float(max(r(z),r(f))))
    return out

def main():
    start=time.time();out=ROOT/'results';out.mkdir(exist_ok=True)
    m=ConditionalModel(); f=lambda x: -(m.loglike(x)+logprior(x))
    opt=minimize(f,[0.,0.],method='Nelder-Mead',options={'xatol':1e-9})
    assert opt.success
    # Independent dense Cholesky check against Woodbury on the actual design.
    x=opt.x;F=np.c_[100*m.Q,3*m.t]
    C=np.diag(np.exp(2*x[0])+np.exp(2*x[1])*m.ev)+F@F.T
    cf=cho_factor(C,lower=True);dense=-.5*(m.n*np.log(2*np.pi)+2*np.log(np.diag(cf[0])).sum()+m.y@cho_solve(cf,m.y))
    difference=abs(dense-m.loglike(x));assert difference<1e-5
    r1=tensor_reference(m,101);r2=tensor_reference(m,201)
    err=abs(r1['logz']-r2['logz'])
    if err>.005:
        r1=r2;r2=tensor_reference(m,301);err=abs(r1['logz']-r2['logz'])
    assert err<.005,('quadrature not converged',err)
    main_summary=summary(r2)
    # Reuse sufficient statistics for prior sensitivity at fixed B; no posterior-as-prior reuse.
    xgrid=r2['x'];st=r2['stats'];l0,j,h=st.T
    prior0=np.array([logprior(x) for x in xgrid]);weight0=r2['weights']
    llmain=l0-.5*np.log1p(9*j)+4.5*h*h/(1+9*j)
    sensitivity=[]
    for ts in [.3,1.,3.]:
        lp=np.array([logprior(x,ts) for x in xgrid]);base=lp-prior0
        z0=r2['logz']+logsumexp(np.log(np.maximum(weight0,1e-320))+base+l0-llmain)
        for sa in [1.,3.,10.]:
            ll=l0-.5*np.log1p(sa*sa*j)+.5*sa*sa*h*h/(1+sa*sa*j)
            lw=np.log(np.maximum(weight0,1e-320))+base+ll-llmain
            normz=logsumexp(lw)
            rr=dict(logz=r2['logz']+normz,x=xgrid,weights=np.exp(lw-normz),mean=sa*sa*h/(1+sa*sa*j),variance=sa*sa/(1+sa*sa*j))
            sensitivity.append(dict(tau_prior_scale=ts*BASE,A_prior_sd=sa*AUNIT,log_B10=rr['logz']-z0,**summary(rr)))
    # Nested sampling directly integrates the likelihood over the declared proper priors.
    import dynesty
    nested=[]
    for seed in [20261118,20261119]:
        sampler=dynesty.NestedSampler(m.loglike,transform,ndim=2,nlive=250,rstate=np.random.default_rng(seed),bound='multi',sample='rwalk')
        sampler.run_nested(dlogz=.05,print_progress=False)
        res=sampler.results
        nested.append(dict(seed=seed,nlive=250,logz=float(res.logz[-1]),logz_error=float(res.logzerr[-1]),calls=int(np.sum(res.ncall)),difference_from_quadrature=float(res.logz[-1]-r2['logz'])))
    # Four genuinely independent Metropolis chains. Proposal freezes before production.
    eps=.001;H=np.empty((2,2))
    for i in range(2):
        for k in range(2):
            ei=np.eye(2)[i]*eps;ek=np.eye(2)[k]*eps
            H[i,k]=(f(opt.x+ei+ek)-f(opt.x+ei-ek)-f(opt.x-ei+ek)+f(opt.x-ei-ek))/(4*eps**2)
    prop=np.linalg.cholesky(np.linalg.inv(H))*2.38/np.sqrt(2)
    chains=[];accept=[]
    for seed in range(4):
        rng=np.random.default_rng(20261200+seed);pos=opt.x+prop@rng.normal(size=2)*3;val=-f(pos);trace=[];nacc=0
        for it in range(14000):
            trial=pos+prop@rng.normal(size=2);new=-f(trial)
            if np.log(rng.uniform())<new-val:pos,val=trial,new;nacc+=int(it>=2000)
            if it>=2000:trace.append(pos.copy())
        chains.append(trace);accept.append(nacc/12000)
    chains=np.asarray(chains);rh=rank_rhat(chains)
    assert max(rh)<1.02,('MCMC mixing',rh)
    import emcee
    taus=[]
    for chain in chains:
        taus.append(emcee.autocorr.integrated_time(chain,has_walkers=False,quiet=True).tolist())
    # Scrambled Sobol digital nets: repeated scrambles estimate integration variation.
    qz=[]
    for seed in range(8):
        points=qmc.Sobol(2,scramble=True,seed=20261300+seed).random_base2(14)
        ll=np.array([m.loglike(transform(q)) for q in points])
        qz.append(float(logsumexp(ll)-np.log(len(ll))))
    # Posterior predictive: Matheron's conditional simulation, including all latent terms.
    rng=np.random.default_rng(20261400);ids=rng.choice(len(xgrid),256,p=weight0)
    obs_stats=[];rep_stats=[]
    def discrepancy(r):
        pairs=[];pm=[]
        for pid in np.unique(m.pid):
            z=r[m.pid==pid];pairs.extend(z[:-1]*z[1:]);pm.append(np.mean(z))
        return [np.mean(pairs),np.max(abs(np.array(pm)))]
    for ix in ids:
        x=xgrid[ix];s=np.exp(x[0]);tau=np.exp(x[1])
        latent=m.Q@(100*rng.normal(size=m.p))+tau*np.sqrt(m.ev)*rng.normal(size=m.n)+m.t*3*rng.normal()
        noise=s*rng.normal(size=m.n);y0=latent+noise
        latent_post=m.y-noise-s*s*m.inverse(x,m.y-y0)
        noise_rep=s*rng.normal(size=m.n)
        obs_stats.append(discrepancy(m.unrotate(m.y-latent_post)/s))
        rep_stats.append(discrepancy(m.unrotate(noise_rep)/s))
    obs_stats=np.array(obs_stats);rep_stats=np.array(rep_stats)
    # Conjugate simulation-based calibration: tests the analytic Gaussian update only.
    a=rng.normal(0,3,4000);y=a+rng.normal(size=4000);means=.9*y;sd=np.sqrt(.9)
    pit=ndtr((a-means)/sd);coverage=float(np.mean((pit>.025)&(pit<.975)))
    assert .93<coverage<.97
    result=dict(scope='Conditional optical amplitude A; fixed lambda=2.81884565 solar radii; not gamma, f_R0 or graviton mass',
        observations=m.n,nuisance_dimension=m.p,kernel_ell_seconds=7200,noise_prior='log(s) ~ Normal(0, 0.35^2); tau ~ HalfNormal(1e-13)',
        amplitude_prior='A ~ Normal(0,0.003^2)',nuisance_prior='beta in fixed whitened orthonormal design ~ Normal(0,100^2 I)',
        prior_status='Explicit sensitivity scenarios; no independent instrument calibration used',
        posterior=main_summary,sensitivity=sensitivity,nested=nested,
        numerical_checks=dict(dense_woodbury_loglike_difference=difference,quadrature_logz_change=err,
            finite_grid_bounds=[[-1.8,2.5],[-11,3]],rank_normalized_split_rhat=rh,mcmc_acceptance=accept,autocorrelation_times=taus,
            MCMC_hyper_mean=np.mean(chains,axis=(0,1)).tolist(),quadrature_hyper_mean=(weight0@xgrid).tolist(),
            sobol_logz=qz,sobol_logz_mean=float(np.mean(qz)),sobol_logz_sd=float(np.std(qz,ddof=1)),
            conjugate_SBC_95_coverage=coverage,SBC_replications=4000,SBC_scope='Conjugate Gaussian submodel only'),
        posterior_predictive=dict(draws=256,statistics=['within-pass adjacent residual product','largest absolute pass mean'],
            tail_probabilities=np.mean(rep_stats>=obs_stats,axis=0).tolist(),scope='Posterior predictive checks are diagnostics, not independent model validation'),
        runtime_seconds=time.time()-start)
    (out/'hierarchical.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(out/'hierarchical_arrays.npz',x=xgrid,weights=weight0,mean=r2['mean'],variance=r2['variance'],chains=chains,ppc_obs=obs_stats,ppc_rep=rep_stats)
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
