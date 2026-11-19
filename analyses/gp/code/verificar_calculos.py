"""Independent numerical checks of the load-bearing finite assertions."""
from pathlib import Path
import sys,json
import numpy as np
import mpmath as mp
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.linalg import cho_factor,cho_solve
R=Path(__file__).resolve().parents[1];OUT=R/'results';sys.path.insert(0,str(R/'code'))
from gp_realizaciones import quantiles_and_mix
from gp_covarianza_cassini import intkernel

checks={}
# Exact piecewise-exponential integration versus adaptive quadrature of a
# separately constructed linearly interpolated Q, including quantile inversion.
x=np.array([0.,.2,.9,1.7,3.]);Q=np.array([[0.,1.,-.4,8.,35.],[1.,-.2,2.,5.,29.]])
qs,qmix,density=quantiles_and_mix(x,Q);refs=[]
for row in Q:
    f=lambda t:np.exp(-np.interp(t,x,row)/2)
    Z=quad(f,0,3,points=x[1:-1],epsabs=1e-12,epsrel=1e-12)[0]
    def CDF(t):return quad(f,0,t,points=x[(x>0)&(x<t)],epsabs=1e-12,epsrel=1e-12)[0]/Z
    refs.append(brentq(lambda t:CDF(t)-.997,0,3,xtol=1e-13))
checks['max_quantile_difference_independent_quadrature']=float(np.max(abs(qs-refs)))
assert checks['max_quantile_difference_independent_quadrature']<1e-10

# Recompute sparse GP predictions in 60 decimal places at fixed fitted ell.
# This tests the reported numerical effect, not uniqueness of hyperparameter MLE.
mp.mp.dps=60;j=json.loads((OUT/'gp_emulator.json').read_text());zhat=mp.mpf(str(j['z_hat']))
refcurve=np.load(OUT/'gp_emulator_curves.npz');errors=[]
for row in j['results']:
    if row['n_training']!=4:continue
    ell=mp.mpf(str(row['ell_coordinate']));coord=row['coordinate'];n=4
    xt=[mp.mpf(i) for i in range(4)]
    z=[t/3 if coord=='x' else (t/3)**2 for t in xt]
    K=mp.matrix([[mp.exp(-((a-b)/ell)**2/2)+(mp.mpf('1e-12') if i==k else 0) for k,b in enumerate(z)] for i,a in enumerate(z)])
    y=mp.matrix([t**4-2*zhat*t*t for t in xt]);one=mp.matrix([1]*n)
    iy=mp.lu_solve(K,y);io=mp.lu_solve(K,one);mu=(one.T*iy)[0]/(one.T*io)[0];alpha=mp.lu_solve(K,y-mu*one)
    maxerr=0.
    for idx in np.linspace(0,len(refcurve['x'])-1,31,dtype=int):
        t=mp.mpf(str(refcurve['x'][idx]));zz=t/3 if coord=='x' else (t/3)**2
        mean=mu+sum(mp.exp(-((zz-z[k])/ell)**2/2)*alpha[k] for k in range(n))
        maxerr=max(maxerr,abs(float(mean)-refcurve[f'mean__4__{coord}'][idx]))
    errors.append(dict(coordinate=coord,max_absolute_Q_difference=maxerr))
assert max(v['max_absolute_Q_difference'] for v in errors)<1e-5
checks['sparse_gp_60_digit_checks']=errors

# Cassini: full dense Cholesky/Schur calculation versus block spectral whitening.
z=np.load(R/'sources/cassini/dss25_diseno_erfa_constant.npz');tt=np.load(R/'sources/cassini/contraste_robusto_pesos.npz')
epoch,X,w,y,pid=(z[k] for k in ['epoch','design','weights','residual','pass_id'])
WX=X*w[:,None];U,s,V=np.linalg.svd(WX/np.linalg.norm(WX,axis=0),full_matrices=False)
tw=tt['constant__template_per_A']*w*1e-3
C=np.eye(len(y));sigma=1e-13
for p in np.unique(pid):
    ids=np.flatnonzero(pid==p);K=intkernel(epoch[ids,None]-epoch[None,ids],7200.)
    C[np.ix_(ids,ids)]+=K*(sigma*w[ids,None])*(sigma*w[None,ids])
cf=cho_factor(C,lower=True);CX=cho_solve(cf,U);Ct=cho_solve(cf,tw)
Mt=Ct-CX@np.linalg.solve(U.T@CX,U.T@Ct);I=tw@Mt
raw=w*(y+tt['constant__delta_SPI']);estimate=1e-3*Mt@raw/I;sd=1e-3/np.sqrt(I)
published=json.loads((OUT/'gp_cassini_refits.json').read_text())
row=next(r for r in published['results'] if r['model']=='constant' and r['ell_seconds']==7200 and r['sigma_GP_count']==1e-13)
checks['dense_vs_block_A_difference']=float(estimate-row['A_hybrid'])
checks['dense_vs_block_sigma_difference']=float(sd-row['sigma_A_conditional'])
assert abs(checks['dense_vs_block_A_difference'])<1e-10
assert abs(checks['dense_vs_block_sigma_difference'])<1e-10
checks['max_joint_injection_deviation']=max(abs(r['joint_injection_recovery']-1) for r in published['results'])
assert checks['max_joint_injection_deviation']<1e-8

# Within pointwise one-sigma bands does not imply substantial joint density.
rho=.999;K=np.array([[1.,rho],[rho,1.]]);v=np.array([1.,-1.])
distance=float(v@np.linalg.solve(K,v));formula=2/(1-rho)
assert np.linalg.eigvalsh(K).min()>0
assert abs(distance/formula-1)<1e-10
checks['pointwise_band_counterexample']={'rho':rho,'pointwise_standardized_offsets':v.tolist(),
    'mahalanobis_squared':distance,'log_density_ratio_to_mean':-.5*distance,
    'scope':'Counterexample to a general implication, not a measurement of the INPOP kernel'}

(OUT/'independent_checks.json').write_text(json.dumps(dict(status='passed',review='self-review with separate numerical formulations',checks=checks),indent=2)+'\n')
print(json.dumps(checks,indent=2))
