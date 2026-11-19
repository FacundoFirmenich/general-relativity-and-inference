"""Scoped information audit: exact model controls, not new gravity observations.

Run with Python + NumPy/SciPy/Matplotlib; no downloads or mutable inputs.
Units: natural logarithms (nats), standardized Gaussian noise unless stated.
"""
from pathlib import Path
import json, sys
import numpy as np
from scipy.integrate import quad
from scipy.special import logsumexp
from scipy.stats import norm, chi2, ncx2

R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(R))
from research.revision_v17.audit import gaussian_distance


def renyi_gaussian(m0, S0, m1, S1, order, rtol=1e-10):
    """D_order(P0||P1), including degenerate Gaussian affine supports.

    For distinct affine supports Gaussian laws are mutually singular. Return
    JSON null + INFINITE, never an undocumented finite regularization.
    Near-rank decisions use the declared relative tolerance; these are numerical
    classifications. For order > 1, finite overlap also requires a PD precision.
    """
    if not np.isfinite(order) or order <= 0:
        raise ValueError('positive finite Renyi order required')
    if order != 1 and abs(order-1) < 1e-7:
        raise ValueError('near-one cancellation: use exact order 1 or higher precision')
    m0, S0, m1, S1 = map(lambda x: np.asarray(x, dtype=float), (m0,S0,m1,S1))
    n = len(m0)
    if m0.shape != (n,) or m1.shape != (n,) or S0.shape != (n,n) or S1.shape != (n,n):
        raise ValueError('incompatible shapes')
    if not all(np.isfinite(x).all() for x in [m0,m1,S0,S1]):
        raise ValueError('nonfinite Gaussian input')
    check = gaussian_distance(m0,S0,m1,S1,rtol)
    out = dict(order=float(order), value=None, status='INFINITE', rank=check['rank'], rank_rtol=rtol)
    if check['support'] == 'different':
        return dict(out, reason='mutually singular affine supports')
    if order == 1:
        return dict(out, value=check['KL'], status='FINITE', reason='KL limit')
    if check['rank'][0] == 0:
        return dict(out, value=0., status='FINITE', reason='identical point laws')
    ev,U = np.linalg.eigh(S0)
    U = U[:,ev > rtol*np.max(abs(ev))]
    A, B, d = U.T@S0@U, U.T@S1@U, U.T@(m1-m0)
    Ai, Bi = np.linalg.inv(A), np.linalg.inv(B)
    a = order
    precision = a*Ai + (1-a)*Bi
    ve = np.linalg.eigvalsh(precision)
    # The tolerance at the finite/infinite boundary is returned explicitly.
    if ve.min() <= rtol*max(np.linalg.norm(a*Ai,2),np.linalg.norm((1-a)*Bi,2)):
        return dict(out, reason='overlap precision not numerically positive definite')
    h = (1-a)*Bi@d
    ld = lambda S: np.linalg.slogdet(S)[1]
    log_overlap = (-a*ld(A)-(1-a)*ld(B)-ld(precision)
                   -(1-a)*d@Bi@d+h@np.linalg.solve(precision,h))/2
    val = log_overlap/(a-1)
    if val < -1e-8:
        raise ArithmeticError('negative divergence beyond roundoff')
    return dict(out, value=float(max(0.,val)), status='FINITE', reason='common-support Gaussian overlap')


def noise_mixture():
    """Unknown label M in {1,3}, equal weights; Y=theta+N(0,M^2).

    Bayes posterior is a two-normal mixture for a normal theta prior.
    Quad integrates predictive posterior variance; independent computation
    uses Var(theta)-E[predictive posterior mean squared].
    """
    sd = np.array([1.,3.]); w = np.array([.5,.5])
    def fi(x):
        logs = np.log(w)+norm.logpdf(x,scale=sd)
        l = logsumexp(logs)
        if l < -740: return 0.
        score = np.sum(np.exp(logs-l)*(-x/sd**2))
        return float(np.exp(l)*score**2)
    J, err = quad(fi,-np.inf,np.inf,epsabs=2e-11,epsrel=2e-11,limit=150)
    oracle_J = float(np.sum(w/sd**2)); rows=[]
    for v in [.25,1.,4.]:
        vj = v*sd**2/(v+sd**2)
        def stats(y):
            logs=np.log(w)+norm.logpdf(y,scale=np.sqrt(v+sd**2))
            lp=logsumexp(logs); z=np.exp(logs-lp)
            mj=v/(v+sd**2)*y
            mean=np.sum(z*mj)
            var=np.sum(z*(vj+(mj-mean)**2))
            return float(np.exp(lp)),float(mean),float(var)
        risk,e = quad(lambda y: stats(y)[0]*stats(y)[2],-np.inf,np.inf,epsabs=2e-11,epsrel=2e-11)
        em2,emerr = quad(lambda y: stats(y)[0]*stats(y)[1]**2,-np.inf,np.inf,epsabs=2e-11,epsrel=2e-11)
        bound=1/(J+1/v); oracle_risk=float(w@vj)
        assert abs(risk-(v-em2)) < 3e-10
        assert bound <= risk+1e-10 and oracle_risk <= risk+1e-10
        rows.append(dict(prior_variance=v,Fisher_mixture=J,prior_Fisher=1/v,
                         Bayes_risk=risk,van_Trees_lower_bound=bound,oracle_Bayes_risk=oracle_risk,
                         retained_K=1-risk/v,oracle_K=1-oracle_risk/v,
                         risk_quad_error=e,variance_identity_error=abs(risk-v+em2)))
    return dict(noise_sd=sd.tolist(),prior_label_weights=w.tolist(),Fisher_mixture=J,
                oracle_Fisher=oracle_J,Fisher_quad_error=err,scenarios=rows)


def rms_control():
    """Exact rejection probabilities, no Monte Carlo approximations."""
    level=.05; h=3.; z=norm.isf(level/2)
    rows=[]
    for n in [16,64,256,1024,4096,16384]:
        theta=h/np.sqrt(n)
        matched=float(norm.sf(z-h)+norm.cdf(-z-h))
        energy=float(ncx2.sf(chi2.isf(level,n),n,n*theta**2))
        rows.append(dict(n=n,theta=theta,matched_power=matched,RMS_power=energy,
                         raw_Fisher_at_null=n,RMS_Fisher_at_null=0.))
    # Verify the noncentral chi-square Poisson-mixture representation using an
    # independent summation of central chi-square tails.
    from scipy.stats import poisson
    n=256; cut=chi2.isf(level,n);k=np.arange(100)
    alternative=float(np.sum(poisson.pmf(k,h*h/2)*chi2.sf(cut,n+2*k)))
    assert abs(alternative-rows[2]['RMS_power']) < 1e-12
    assert all(rows[i+1]['RMS_power']<rows[i]['RMS_power'] for i in range(len(rows)-1))
    return dict(level=level,local_amplitude_h=h,rows=rows,
                noncentral_tail_independent_error=abs(alternative-rows[2]['RMS_power']),
                scope='Yi iid N(theta,1), known variance; two-sided matched test and upper-tail energy test',
                exact_theorem='For Z~N(theta*u,I), J_norm(theta=0)=0; J_raw=||u||^2. A symmetric theta prior has E(theta|norm Z)=0 and quadratic K=0.')


def mantissa_control():
    constants={'ln10':float(np.log(10)),'ln9':float(np.log(9)),
               'ln20_95pct_saturation':float(np.log(20)),
               'ln19_5_95_window':float(np.log(19)),
               'third_decade':float(10**(1/3)),'two_pi_over_e':float(2*np.pi/np.e),
               'retain_90pct_exponential_time_over_tau':float(-np.log(.9)),
               'logistic_10_90_width_over_tau':float(2*np.log(9))}
    qs={str(p):float(norm.ppf(p)) for p in [.97,.98,.985,.99]}
    fractions=[.5,.9,.95,.99,.999]
    threshold=[dict(fraction=f,beta=b,time_over_tau=float((-np.log1p(-f))**(1/b)))
               for b in [.5,1.,2.] for f in fractions]
    # Counterexample to universality: logistic and stretched exponentials.
    assert constants['logistic_10_90_width_over_tau'] > 4
    benford=dict(band_probability=float(np.log10(2.3/2.1)),
                 lower_equal_width_probability=float(np.log10(1.2)),
                 equal_mass_upper_endpoint=float(8*2.3/2.1),
                 density_derivative='-1/(m^2*ln(10)) < 0 on [1,10)',
                 scope='Conditional on Benford significands; not a test on selected reported constants')
    fits=[dict(solve_for=9,estimate=2.1,standard_uncertainty=2.3),
          dict(solve_for=12,estimate=1.35,standard_uncertainty=2.47),
          dict(solve_for=8,estimate=.21,standard_uncertainty=2.43)]
    for r in fits:
        r['unit']='1e-5 in gamma-1'
        r['standardized_estimate']=r['estimate']/r['standard_uncertainty']
    scales=[1.,np.log(10),1/2.3,10.]
    conversion=[dict(scale=float(s),estimate=2.1*s,standard_uncertainty=2.3*s,
                     z=float(2.1*s/(2.3*s))) for s in scales]
    assert max(abs(r['z']-2.1/2.3) for r in conversion)<1e-14
    return dict(constants=constants,gaussian_quantiles=qs,stretched_exponential=threshold,
                benford=benford,published_Cassini_fits=fits,reparameterization=conversion,
                physical_mantissa_artifact='NOT_ESTABLISHED',
                source='Bertotti, Iess, Tortora (2003), doi:10.1038/nature01997; Fig.3 and Method, pp.375-376; author-uploaded full text read',
                original_main_points=1094,original_main_solve_for=9,reconstructed_epochs=1696,
                boundary='No published observation refit; no derivation of sigma_gamma from ln(10).')


def calibration_control():
    """Full prior-predictive calibration of this two-component toy hierarchy.
    Independent exact posterior CDF; not SBC of the Cassini hierarchy.
    """
    seed=2026111818; n=4000; v=1.; rng=np.random.default_rng(seed)
    theta=rng.normal(0,np.sqrt(v),n)
    sd=rng.choice([1.,3.],n); y=theta+rng.normal(size=n)*sd
    components=np.array([1.,3.])
    logs=norm.logpdf(y[:,None],scale=np.sqrt(v+components**2))
    weights=np.exp(logs-logsumexp(logs,axis=1)[:,None])
    means=y[:,None]*v/(v+components**2)
    post_sd=np.sqrt(v*components**2/(v+components**2))
    pit=np.sum(weights*norm.cdf((theta[:,None]-means)/post_sd),axis=1)
    bins=np.histogram(pit,bins=np.linspace(0,1,11))[0]
    cover={str(level):dict(fraction=float(np.mean((pit>=(1-level)/2)&(pit<=(1+level)/2))),
                           binomial_mc_se_under_nominal=float(np.sqrt(level*(1-level)/n))) for level in [.9,.95]}
    # Fail loudly on a gross implementation-calibration discrepancy; these are
    # probabilistic diagnostics at a fixed seed, not proofs of correctness.
    for level,row in cover.items():
        assert abs(row['fraction']-float(level))<6*row['binomial_mc_se_under_nominal']
    for k in [1,2,5,20]:
        C=np.ones((k,k));one=np.ones(k)
        assert abs(one@np.linalg.pinv(C)@one-1)<1e-12
    return dict(seed=seed,replicates=n,prior_variance=v,PIT_decile_counts=bins.tolist(),coverage=cover,
                source='Synthetic draws from the declared prior and sampling distribution',
                scope='Calibration of exact posterior CDF for the full two-label toy hierarchy; not Cassini SBC',
                duplicate_information=dict(copies=[1,2,5,20],correct_Fisher=[1.,1.,1.,1.],false_independent_Fisher=[1.,2.,5.,20.]))


def tail_controls():
    """Normal location experiment at theta=0, fixed thresholds [-c,c].
    Conditional truncation, selection indicator, and two-sided censoring are
    different observation schemes. The index here is Fisher retention only.
    """
    rows=[]
    for mass in [.8,.9,.95,.99]:
        c=norm.ppf((1+mass)/2);phi=norm.pdf(c);tail=norm.sf(c)
        inside=mass-2*c*phi
        numeric,err=quad(lambda y:y*y*norm.pdf(y),-c,c,epsabs=1e-13)
        censor=inside+2*phi**2/tail
        assert abs(numeric-inside)<1e-12 and 0<inside<censor<=1
        rows.append(dict(central_probability=mass,threshold=c,retained_inside_Fisher_per_attempt=inside,
                         truncated_Fisher_per_retained_record=inside/mass,
                         two_sided_censored_Fisher_per_attempt=censor,
                         removed_tail_Fisher=1-inside,independent_quad_error=abs(numeric-inside)))
    return dict(rows=rows,full_Fisher=1.,
                scope='N(theta,1) at theta=0; thresholds fixed in measured coordinates, not centered on unknown theta',
                truncation='Conditional retained values with selected count: acceptance-indicator score vanishes at theta=0 only',
                censoring='Retain left/right censor labels and all trial counts; equivalent to clipping with endpoint atoms',
                transport='Known invertible theta-independent maps preserve experiment; mass conservation alone is insufficient')


def main():
    D=np.c_[np.eye(2),-np.ones(2)]
    S0=D@(.1*np.eye(3)+.9*np.ones((3,3)))@D.T;S1=D@D.T
    out=[]
    for orientation,A,B in [('coherent_to_independent',S0,S1),('independent_to_coherent',S1,S0)]:
        for order in [.5,1.,2.]:
            out.append(dict(orientation=orientation,**renyi_gaussian(np.zeros(2),A,np.zeros(2),B,order)))
    checks=[]
    for a in [.25,.5,1.,2.,3.]:
        x=renyi_gaussian([0],[[1]],[.7],[[1]],a)
        assert abs(x['value']-a*.7**2/2)<1e-10
    checks.append('equal-covariance analytic Renyi for five orders')
    H=gaussian_distance(np.zeros(2),S0,np.zeros(2),S1)['H2']
    assert abs(np.exp(-out[0]['value']/2)-(1-H))<1e-12
    checks.append('order-half equals independently computed Hellinger affinity')
    for a in [.5,2.]:
        A=.8;B=1.1;mu=.4
        direct=quad(lambda x: np.exp(a*norm.logpdf(x,scale=np.sqrt(A))+(1-a)*norm.logpdf(x,loc=mu,scale=np.sqrt(B))),-np.inf,np.inf,epsabs=2e-11)[0]
        formula=renyi_gaussian([0],[[A]],[mu],[[B]],a)['value']
        assert abs(np.log(direct)/(a-1)-formula)<1e-10
    checks.append('unequal-covariance overlap independently integrated')
    M=np.array([[1.,2.],[-1,1.]])
    for a in [.5,2.]:
        base=renyi_gaussian(np.zeros(2),S0,np.zeros(2),S1,a)['value']
        alt=renyi_gaussian(np.zeros(2),M@S0@M.T,np.zeros(2),M@S1@M.T,a)['value']
        tiny=renyi_gaussian(np.zeros(2),1e-24*S0,np.zeros(2),1e-24*S1,a)['value']
        assert abs(base-alt)<1e-10 and abs(base-tiny)<1e-10
    checks.append('invertible contrasts and covariance units preserve divergence')
    assert renyi_gaussian([0,0],np.diag([1.,0]),[0,1],np.diag([1.,0]),.5)['status']=='INFINITE'
    singular=renyi_gaussian([0,0],np.diag([1.,0]),[.7,0],np.diag([1.,0]),2.)
    assert abs(singular['value']-.49)<1e-12
    assert out[-1]['status']=='INFINITE' and out[-2]['status']=='FINITE'
    checks.append('common singular support and infinite second-moment case')
    result=dict(scope='Exact statistical controls and source-specific arithmetic; no new physical observations',
                logarithm_base='natural',renyi=out,mixture=noise_mixture(),rms=rms_control(),mantissas=mantissa_control(),calibration=calibration_control(),tails=tail_controls(),checks=checks)
    path=R/'research/revision_v18/results';path.mkdir(exist_ok=True)
    (path/'information_audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2,allow_nan=False))
    figures(result)


def figures(data=None):
    if data is None:data=json.loads((R/'research/revision_v18/results/information_audit.json').read_text())
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors=['#18A7D1','#7653C7','#9A5498','#17365D']
    for lang in ['en','es']:
        en=lang=='en';fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
        n=np.unique(np.round(np.geomspace(8,20000,140))).astype(int);h=3.;a=.05
        axs[0].semilogx(n,[ncx2.sf(chi2.isf(a,int(k)),int(k),h*h) for k in n],color=colors[1],label='RMS / energy' if en else 'RMS / energía')
        z=norm.isf(a/2);p=norm.sf(z-h)+norm.cdf(-z-h)
        axs[0].axhline(p,color=colors[0],label='Signed matched response' if en else 'Respuesta con signo')
        axs[0].axhline(a,color=colors[2],linestyle=':',label='Null level 0.05' if en else 'Nivel nulo 0,05')
        axs[0].set(xlabel='Number of observations n' if en else 'Número de observaciones n',ylabel='Power' if en else 'Potencia',title=r'$\theta=3/\sqrt{n}$ · known noise' if en else r'$\theta=3/\sqrt{n}$ · ruido conocido',ylim=(0,1))
        axs[0].legend(frameon=False,fontsize=8)
        vs=[x['prior_variance'] for x in data['mixture']['scenarios']]
        for key,label,c,ls in [('retained_K','Unknown mixture label' if en else 'Etiqueta de mezcla desconocida',colors[1],'-'),('oracle_K','Observed label (oracle)' if en else 'Etiqueta observada (oráculo)',colors[0],'--')]:
            axs[1].plot(vs,[x[key] for x in data['mixture']['scenarios']],marker='o',color=c,linestyle=ls,label=label)
        axs[1].set(xscale='log',xlabel='Prior variance v' if en else 'Varianza previa v',ylabel=r'$\mathcal{K}_{\pi,L}$',title='Mixture information and its oracle' if en else 'Información de la mezcla y su oráculo',ylim=(0,1))
        axs[1].legend(frameon=False,fontsize=8)
        for ax in axs:ax.spines[['top','right']].set_visible(False)
        fig.savefig(R/f'figures/information_retention_{lang}.pdf');plt.close(fig)
        fig,axs=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
        f=np.linspace(.05,.99,200)
        for beta,c in zip([.5,1,2],colors):
            axs[0].plot(f,(-np.log1p(-f))**(1/beta),color=c,label=rf'$\beta={beta}$')
        axs[0].axvline(.9,color=colors[3],linestyle=':');axs[0].axhspan(2.1,2.3,color=colors[2],alpha=.18)
        axs[0].set(xlabel='Saturation fraction' if en else 'Fracción de saturación',ylabel=r'$t/\tau$',ylim=(0,6),title=r'$1-\exp[-(t/\tau)^\beta]$');axs[0].legend(frameon=False)
        m=np.linspace(1,10,400);axs[1].plot(m,1/(m*np.log(10)),color=colors[1])
        axs[1].fill_between(m,0,1/(m*np.log(10)),where=(m>=2.1)&(m<=2.3),color=colors[2],alpha=.4)
        axs[1].set(xlabel='Significand m' if en else 'Significando m',ylabel='Density' if en else 'Densidad',title='Benford: monotone, no peak at 2.2' if en else 'Benford: monótona, sin máximo en 2,2')
        for ax in axs:ax.spines[['top','right']].set_visible(False)
        fig.savefig(R/f'figures/threshold_audit_{lang}.pdf');plt.close(fig)


if __name__ == '__main__':
    if '--figures' in sys.argv:figures()
    else:main()
