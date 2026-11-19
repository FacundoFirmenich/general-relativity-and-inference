"""Actual Cassini design: GP covariance refits versus null-trained GP subtraction.
Kernels and amplitudes are sensitivity hypotheses, not inferred instrument noise.
"""
from pathlib import Path
import json,csv,hashlib
import numpy as np
from scipy.linalg import eigh
from scipy.special import erf,roots_legendre
R=Path(__file__).resolve().parents[1];SRC=R/'sources/cassini';OUT=R/'results'
AUNIT=1e-3;T=300.;BASE=1e-13

def intkernel(dt,ell):
    def F(x):return ell**2*np.exp(-.5*(x/ell)**2)+ell*x*np.sqrt(np.pi/2)*erf(x/(np.sqrt(2)*ell))
    norm=2*(F(T)-F(0))
    return (F(dt+T)-2*F(dt)+F(dt-T))/norm

def main():
    template=np.load(SRC/'contraste_robusto_pesos.npz');rows=[];checks=[];saved={}
    for model in ['constant','four_batches']:
        z=np.load(SRC/f'dss25_diseno_erfa_{model}.npz')
        epoch,X,w,y,pid=(z[k] for k in ['epoch','design','weights','residual','pass_id'])
        t=template[model+'__template_per_A'];delta=template[model+'__delta_SPI']
        WX=X*w[:,None];Q,s,V=np.linalg.svd(WX/np.linalg.norm(WX,axis=0),full_matrices=False)
        assert s[-1]>1e-8
        tw=w*t*AUNIT;p=tw-Q@(Q.T@tw);data=np.c_[w*y,w*(y+delta)]
        z0=data-Q@(Q.T@data)
        cases=[(0.,0.)]+[(ell,sig) for ell in [300.,1800.,7200.,86400.] for sig in [3e-14,1e-13,3e-13]]
        cache={}
        for ell,sigma in cases:
            groups=[]
            if ell and ell not in cache:
                for passno in np.unique(pid):
                    ids=np.flatnonzero(pid==passno);dt=epoch[ids,None]-epoch[None,ids]
                    K=intkernel(dt,ell);eig,U=eigh(K*(BASE*w[ids,None])*(BASE*w[None,ids]))
                    assert eig[0]>-1e-6*max(1.,eig[-1]),(ell,eig[0])
                    groups.append((ids,np.maximum(eig,0),U))
                cache[ell]=groups
            groups=cache.get(ell,[])
            def power(v,exponent):
                ans=v.copy()
                for ids,eig,U in groups:
                    f=(1+(sigma/BASE)**2*eig)**exponent
                    ans[ids]=(U*f)@(U.T@v[ids])
                return ans
            BQ=power(Q,-.5);Bt=power(tw,-.5);By=power(data,-.5)
            U,ss,V=np.linalg.svd(BQ,full_matrices=False)
            pp=Bt-U@(U.T@Bt);info=pp@pp
            h=AUNIT*power(pp,-.5)/info
            fit=h@data;se=AUNIT/np.sqrt(info)
            naive=AUNIT*(p@power(z0,-1))/(p@p)
            attenuation=float(p@power(p,-1)/(p@p))
            injection=3e-3
            naive_shift=AUNIT*p@(power(z0[:,0]+injection/AUNIT*p,-1)-power(z0[:,0],-1))/(p@p)
            joint_shift=float(h@(data[:,0]+injection*w*t)-h@data[:,0])
            assert abs(joint_shift/injection-1)<1e-8
            assert abs(naive_shift/injection-attenuation)<1e-8
            assert -.000001<=attenuation<=1.000001
            row=dict(model=model,ell_seconds=ell,sigma_GP_count=sigma,n=len(epoch),
                A_ODF=float(fit[0]),A_hybrid=float(fit[1]),sigma_A_conditional=se,
                A_naive_subtract_ODF=float(naive[0]),A_naive_subtract_hybrid=float(naive[1]),
                injected_amplitude=injection,naive_injection_recovery=attenuation,
                joint_injection_recovery=joint_shift/injection,
                information_A=float(info/AUNIT**2),
                nuisance_constraint_max=float(np.max(abs(h@Q))),template_constraint_error=abs(float(h@(w*t))-1))
            rows.append(row)
            if ell==7200 and sigma==1e-13:saved[model+'__joint_raw_weights']=h*w
        # The integrated continuous kernel is checked with separate quadrature.
        x,g=roots_legendre(32);tt=T*(x+1)/2;gg=g/2
        for ell in [300.,1800.,7200.,86400.]:
            def integ(d):return np.sum(gg[:,None]*gg[None,:]*np.exp(-.5*((d+tt[:,None]-tt[None,:])/ell)**2))
            for dt in [0.,300.,1200.,7200.]:
                ref=integ(dt)/integ(0);got=float(intkernel(np.array(dt),ell))
                checks.append(abs(got-ref));assert abs(got-ref)<1e-8
    with (OUT/'gp_cassini_refits.csv').open('w') as f:
        wri=csv.DictWriter(f,fieldnames=list(rows[0]));wri.writeheader();wri.writerows(rows)
    info=dict(assertion_tested='Conditional GLS refits and exact injection response of a null-trained GP subtraction',
        scope='Fixed observational design and empirical baseline weights; no observed calibration covariance or graviton bound',
        kernel='Squared-exponential continuous process integrated over 300s counts, normalized to constant count variance, independent by pass',
        max_integrated_kernel_absolute_error=max(checks),results=rows)
    (OUT/'gp_cassini_refits.json').write_text(json.dumps(info,indent=2)+'\n')
    np.savez_compressed(OUT/'gp_cassini_weights.npz',**saved)
    for row in rows:print(row['model'],row['ell_seconds'],row['sigma_GP_count'],row['A_hybrid'],row['sigma_A_conditional'],row['naive_injection_recovery'])
if __name__=='__main__':main()
