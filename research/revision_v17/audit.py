"""Finite Gaussian contrast audit and exact analytic controls; no new observations."""
from pathlib import Path
import sys,json,hashlib,time
import numpy as np
from scipy.linalg import eigh
from scipy.special import ndtr,ndtri
from scipy.integrate import quad
from scipy.stats import norm
R=Path(__file__).resolve().parents[2]

def gaussian_distance(m0,S0,m1,S1,rtol=1e-10):
    """KL(P0||P1), Hellinger squared = 1-affinity. Singular support explicit.

    Numerical ranks are conditional on rtol; jitter is never silently inserted.
    Exact support inequality gives mutual singularity. Near equality is tolerance
    classified, so rank/tolerance are part of the audit, not a physical theorem.
    """
    m0,m1,S0,S1=map(np.asarray,(m0,m1,S0,S1))
    bases=[];ranks=[]
    for S in [S0,S1]:
        if not np.allclose(S,S.T,rtol=rtol,atol=rtol*np.linalg.norm(S,2)):raise ValueError('asymmetric covariance')
        v,U=eigh(S);tol=rtol*np.max(abs(v))
        if v.min() < -tol:raise ValueError('negative covariance eigenvalue')
        bases.append(U[:,v>tol]);ranks.append(int(np.sum(v>tol)))
    U,V=bases;dm=m1-m0
    common=(ranks[0]==ranks[1] and np.linalg.norm(U@U.T-V@V.T)<rtol*10 and np.linalg.norm(dm-U@(U.T@dm))<rtol*10*max(np.linalg.norm(m0),np.linalg.norm(m1),np.sqrt(np.linalg.norm(S0,2)),np.sqrt(np.linalg.norm(S1,2)),np.finfo(float).tiny))
    if not common:return dict(KL=None,KL_status='INFINITE',H2=1.,rank=ranks,support='different',rank_rtol=rtol)
    if not ranks[0]:return dict(KL=0.,KL_status='FINITE',H2=0.,rank=ranks,support='common point',rank_rtol=rtol)
    A=U.T@S0@U;B=U.T@S1@U;d=U.T@dm;M=(A+B)/2
    la=np.linalg.slogdet(A)[1];lb=np.linalg.slogdet(B)[1];lm=np.linalg.slogdet(M)[1]
    kl=.5*(np.trace(np.linalg.solve(B,A))+d@np.linalg.solve(B,d)-len(d)+lb-la)
    logaff=.25*(la+lb)-.5*lm-.125*d@np.linalg.solve(M,d)
    return dict(KL=float(max(0,kl)),KL_status='FINITE',H2=float(max(0,-np.expm1(logaff))),rank=ranks,support='common',rank_rtol=rtol)

def contrast_audit(b0,K0,b1,K1,N=1,weights=None,rtol=1e-10):
    n=len(b0)
    if n<2 or N<=0 or not np.isfinite(N):raise ValueError('n >= 2 and N > 0 required')
    weights=np.ones(n) if weights is None else np.asarray(weights)
    if weights.shape!=(n,) or np.any(weights<=0) or not np.isfinite(weights).all():raise ValueError('fixed positive finite weights required')
    D=np.c_[np.eye(n-1),-np.ones(n-1)]
    b0,b1,K0,K1=map(np.asarray,(b0,b1,K0,K1))
    if b0.shape!=(n,) or b1.shape!=(n,) or K0.shape!=(n,n) or K1.shape!=(n,n):raise ValueError('incompatible shapes')
    if not all(np.isfinite(a).all() for a in [b0,b1,K0,K1]):raise ValueError('nonfinite input')
    for K in [K0,K1]:
        if not np.allclose(K,K.T,rtol=rtol,atol=rtol*np.linalg.norm(K,2)) or np.linalg.eigvalsh(K).min() < -rtol*np.linalg.norm(K,2):raise ValueError('invalid full cost covariance')
    mu0=-N/2*(D@b0);mu1=-N/2*(D@b1)
    S0=(N/2)**2*(D@K0@D.T);S1=(N/2)**2*(D@K1@D.T)
    return dict(contrast_mean_reference=mu0.tolist(),contrast_covariance_reference=S0.tolist(),contrast_mean_alternative=mu1.tolist(),contrast_covariance_alternative=S1.tolist(),**gaussian_distance(mu0,S0,mu1,S1,rtol))

def main():
    cases={};n=4;b=np.arange(n)*.2;K=np.eye(n)
    cases['common_mode']=contrast_audit(b,K,b+3,K+2*np.ones((n,n)))
    cases['covariance_loss']=contrast_audit(np.zeros(3),np.full((3,3),.9)+.1*np.eye(3),np.zeros(3),np.eye(3))
    cases['singular_common_support']=gaussian_distance(np.zeros(2),np.diag([1,0]),np.array([1,0]),np.diag([2,0]))
    cases['singular_different_support']=gaussian_distance(np.zeros(2),np.diag([1,0]),np.array([0,1]),np.diag([1,0]))
    assert cases['common_mode']['H2']<1e-12 and cases['common_mode']['KL']<1e-12
    assert cases['singular_different_support']['KL_status']=='INFINITE'
    c=gaussian_distance(np.array([0.]),np.array([[1.]]),np.array([.7]),np.array([[1.]]))
    assert abs(c['KL']-.7**2/2)<1e-12 and abs(c['H2']-(1-np.exp(-.7**2/8)))<1e-12
    # Bijection invariance checked using a different contrast reference basis.
    D=np.c_[np.eye(2),-np.ones(2)];M=np.array([[1.,2.],[-1,1.]])
    C=D@(np.full((3,3),.9)+.1*np.eye(3))@D.T;B=D@np.eye(3)@D.T
    cc=gaussian_distance(np.zeros(2),M@C@M.T,np.zeros(2),M@B@M.T)
    assert abs(cc['H2']-cases['covariance_loss']['H2'])<1e-12
    scaled=gaussian_distance(np.zeros(2),C*1e-24,np.zeros(2),B*1e-24)
    assert scaled['rank']==[2,2] and abs(scaled['H2']-cc['H2'])<1e-12
    alpha0=.04;alpha_star=.05;z=ndtri(1-alpha0);rows=[]
    for shift in [.02,.05,.15]:
        tv=2*ndtr(abs(shift)/2)-1
        numeric=quad(lambda x:abs(norm.pdf(x)-norm.pdf(x-shift))/2,-12,12,epsabs=1e-12,points=[shift/2])[0]
        assert abs(tv-numeric)<1e-10
        rows.append(dict(bias=shift,sd=1.,TV=tv,risk=norm.sf(z-shift),uniform_upper=alpha0+tv,certified=bool(alpha0+tv<=alpha_star)))
    memory=dict(fock_probability=(1-np.exp(-1))/2,thermal_probability=(1+np.exp(-3))/2,scope='Exact characteristic-function counterexample to moment closure, not GP regression or instrument data')
    result=dict(scope='Analytic model controls and reusable Gaussian contrast-law audit',cases=cases,reference_rejection=alpha0,target_rejection=alpha_star,threshold=z,mean_error_budget=2*ndtri((1+alpha_star-alpha0)/2),transport=rows,memory=memory,checks=['Common-mode invariance','Singular support rejection','One-dimensional analytic KL/Hellinger','Invertible contrast-coordinate invariance','Covariance unit-rescaling invariance','TV independently integrated for three shifts'])
    out=R/'research/revision_v17/results';out.mkdir(exist_ok=True)
    (out/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    figures(result)

def figures(result):
    alpha0=result['reference_rejection'];z=result['threshold']
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    cmap=LinearSegmentedColormap.from_list('allowed',['#18A7D1','#7653C7','#9A5498'])
    bb,ss=np.meshgrid(np.linspace(-.10,.17,60),np.linspace(.8,1.2,60));risk=norm.sf((z-bb)/ss)
    for lang in ['en','es']:
        en=lang=='en';fig=plt.figure(figsize=(10,4.4),layout='constrained');ax=fig.add_subplot(121,projection='3d')
        ax.plot_surface(bb,ss,risk,cmap=cmap,linewidth=0,alpha=.94,rasterized=True)
        ax.contour(bb,ss,risk,levels=[.05],colors=['#17365D'],linewidths=2)
        ax.set(xlabel='Mean shift b' if en else 'Sesgo b',ylabel='Noise scale σ' if en else 'Escala de ruido σ',zlabel='Rejection probability' if en else 'Probabilidad de rechazo',title='A fixed rule under calibration change' if en else 'Regla fija ante cambios de calibración');ax.view_init(25,-125)
        ax=fig.add_subplot(122);xs=np.linspace(0,.17,200)
        ax.plot(xs,alpha0+2*ndtr(xs/2)-1,color='#7653C7',label='Uniform TV bound' if en else 'Cota uniforme TV',linestyle='--')
        ax.plot(xs,norm.sf(z-xs),color='#18A7D1',label='Exact rule risk' if en else 'Riesgo exacto de la regla')
        ax.axhline(.05,color='#9A5498',label='Required level' if en else 'Nivel exigido',linestyle=':')
        ax.set(xlabel='Mean shift b' if en else 'Sesgo b',ylabel='Rejection probability' if en else 'Probabilidad de rechazo',title='Certificate and rule-specific calculation' if en else 'Certificado y cálculo específico');ax.legend(frameon=False,fontsize=8);ax.spines[['top','right']].set_visible(False)
        fig.savefig(R/f'figures/governance_transport_{lang}.pdf');plt.close(fig)
if __name__=='__main__':main()
