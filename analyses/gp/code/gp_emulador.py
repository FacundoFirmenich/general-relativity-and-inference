"""Controlled GP emulation of a Cassini-derived conditional optical profile.

x=sqrt(A/sigma_A) >= 0 is dimensionless, NOT graviton mass. N=1696.
The exact profile is known, making interpolation checks independent of a GP band.
"""
from pathlib import Path
import json,csv
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scipy.optimize import minimize_scalar
from scipy.integrate import cumulative_trapezoid
R=Path(__file__).resolve().parents[1];OUT=R/'results'
N=1696

def cdf_quantile(x,cost,p=.997):
    f=np.exp(-.5*(cost-np.min(cost)));cdf=cumulative_trapezoid(f,x,initial=0);cdf/=cdf[-1]
    return float(np.interp(p,cdf,x)),cdf

def gp_fit(x,y,coordinate):
    z=x/3 if coordinate=='x' else (x/3)**2
    dist=z[:,None]-z[None,:];one=np.ones(len(x));jitter=1e-12
    def fit(logell,details=False):
        ell=np.exp(logell);K=np.exp(-.5*(dist/ell)**2)+jitter*np.eye(len(x))
        chol=cho_factor(K,lower=True,check_finite=True)
        inv1=cho_solve(chol,one);invy=cho_solve(chol,y)
        mu=(one@invy)/(one@inv1);alpha=cho_solve(chol,y-mu)
        var=(y-mu)@alpha/len(x)
        if var<=0:return np.inf
        obj=.5*len(x)*np.log(var)+np.log(np.diag(chol[0])).sum()
        return (ell,mu,var,alpha,chol,inv1) if details else obj
    # Deterministic broad scan plus local refinements, avoiding one starting point.
    scan=np.linspace(-3,3,121);vals=np.array([fit(v) for v in scan]);candidates=[]
    for idx in np.argsort(vals)[:5]:
        lo=scan[max(0,idx-1)];hi=scan[min(len(scan)-1,idx+1)]
        if hi>lo:candidates.append(minimize_scalar(fit,bounds=(lo,hi),method='bounded',options={'xatol':1e-10}))
    opt=min(candidates,key=lambda r:r.fun)
    ell,mu,var,alpha,chol,inv1=fit(opt.x,True)
    def pred(xx):
        zz=xx/3 if coordinate=='x' else (xx/3)**2
        Kstar=np.exp(-.5*((zz[:,None]-z[None,:])/ell)**2)
        mean=mu+Kstar@alpha
        v=var*(1-np.einsum('ij,ji->i',Kstar,cho_solve(chol,Kstar.T))+(1-Kstar@inv1)**2/(one@inv1))
        return mean,np.sqrt(np.maximum(v,0))
    # GP derivative at the physical boundary x=0.
    if coordinate=='x':
        rprime=z*np.exp(-.5*(z/ell)**2)/(3*ell**2)
        dm=float(rprime@alpha)
        dv=var*(1/(9*ell**2)-rprime@cho_solve(chol,rprime)+(rprime@inv1)**2/(one@inv1))
    else:dm=0.;dv=0.
    return pred,dict(ell_coordinate=ell,process_variance=var,constant_mean=mu,jitter=jitter,
        derivative_at_zero=dm,derivative_standard_deviation_at_zero=float(np.sqrt(max(dv,0))),neg_loglik_profile=float(opt.fun))

def main():
    src=json.loads((R/'sources/cassini/contraste_robusto.json').read_text())
    row=next(v for v in src['results'] if v['model']=='constant' and v['estimator']=='unprotected')
    zhat=row['A_hybrid_SPI']/row['nominal_sigma_A'];x=np.linspace(0,3,12001)
    exact=x**4-2*zhat*x*x
    qtrue,cdf=cdf_quantile(x,exact);rows=[];curves={'x':x,'delta_Q_exact':exact}
    for n in [4,7,13]:
        xt=np.linspace(0,3,n);yt=xt**4-2*zhat*xt**2
        for coordinate in ['x','x_squared']:
            pred,meta=gp_fit(xt,yt,coordinate);mean,sd=pred(x);estimate,pg=cdf_quantile(x,mean)
            em=mean-exact
            # Differentiate log likelihood relative to boundary; constants cancel.
            logerr=-.5*(em-em[0]);region=x<=qtrue
            training_error=np.max(abs(pred(xt)[0]-yt))
            rd=dict(n_training=n,coordinate=coordinate,exact_quantile_997=qtrue,gp_quantile_997=estimate,
                relative_quantile_error=estimate/qtrue-1,
                max_abs_log_odds_error_to_exact_quantile=float(np.max(abs(logerr[region]))),
                max_abs_Q_error_full_domain=float(np.max(abs(em))),
                max_abs_Q_error_training=float(training_error),
                max_gp_sd_Q=float(sd.max()),**meta)
            rows.append(rd);curves[f'mean__{n}__{coordinate}']=mean;curves[f'sd__{n}__{coordinate}']=sd
    # Mesh refinement of the exact reference quantile.
    fine=np.linspace(0,3,24001);refine=cdf_quantile(fine,fine**4-2*zhat*fine*fine)[0]
    assert abs(refine-qtrue)<1e-7
    np.savez_compressed(OUT/'gp_emulator_curves.npz',**curves)
    obj=dict(assertion_tested='Boundary and tail distortions of a fitted likelihood emulator on a known exact conditional profile',
        physical_scope='x=sqrt(A/sigma_A) is a dimensionless diagnostic coordinate on nonnegative optical amplitudes; no conversion to graviton mass',
        profile_definition='DeltaQ(x)=x^4-2*z_hat*x^2; posterior proportional to exp(-DeltaQ/2), uniform prior in x on [0,3]',
        z_hat=zhat,sigma_A=row['nominal_sigma_A'],prior='Uniform in x for BOTH kernels, so kernel-coordinate changes are not prior changes',
        exact_quantile_mesh_difference=abs(refine-qtrue),N=N,results=rows)
    (OUT/'gp_emulator.json').write_text(json.dumps(obj,indent=2)+'\n')
    with (OUT/'gp_emulator.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(obj,indent=2))
if __name__=='__main__':main()
