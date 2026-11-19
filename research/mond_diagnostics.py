"""Reproducible MOND inference controls; scope is explicit in each output.

Uses public SPARC reduced data for a one-amplitude statistical counterexample,
not for reproducing the 900-parameter RAR analysis or a gravity model ranking.
"""
from pathlib import Path
import json, hashlib, datetime
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import norm, t
from mond_quadrupole import quadrupole

R=Path(__file__).resolve().parent

def nu(y,delta):
    return np.exp(-np.log(-np.expm1(-np.asarray(y)**(delta/2)))/delta)

def physical_external(delta,ge=2.142e-10,a0=1.02e-10,n=192):
    eN=brentq(lambda y:y*nu(y,delta)-ge/a0,1e-7,10.,xtol=1e-13)
    result=quadrupole(n,eN=eN,delta=delta,a0=a0)
    result['physical_external_m_s2']=ge
    result['external_equation_relative_error']=abs(eN*nu(eN,delta)/(ge/a0)-1)
    return result

def main():
    source=R/'raw/sparc_massmodels.mrt'
    rows=[s.split() for s in source.read_text().splitlines()[25:] if len(s.split())==10]
    galaxy='NGC2403'  # chosen before fitting: well resolved, bulgeless, quality 1
    a=np.array([list(map(float,r[1:])) for r in rows if r[0]==galaxy])
    D,radius,v,sv,vg,vd,vb,_,_=a.T
    assert len(a)>20 and np.all(sv>0) and np.all(vb==0)
    a0=1.2e-10;ml=.5
    vN2=vg*np.abs(vg)+ml*vd*vd
    gN=vN2*1e6/(radius*3.085677581491367e19)
    template=np.sqrt(vN2*nu(gN/a0,1.))
    # y=A*template+epsilon; C diagonal and every physical nuisance fixed.
    tau=.1;residual=v-template;S=float(np.sum((template/sv)**2))
    b=float(np.sum(template*residual/sv**2))
    logL0=float(np.sum(norm.logpdf(v,template,sv)))
    variance=1/(S+tau**-2);mu=1+b*variance;sd=np.sqrt(variance)
    logBF=-.5*np.log1p(tau*tau*S)+.5*b*b*variance
    logZ=logL0+logBF
    def logjoint(A):
        return float(np.sum(norm.logpdf(v,A*template,sv))+norm.logpdf(A,1,tau))
    peak=logjoint(mu)
    exact_density_peak=norm.logpdf(mu,mu,sd)
    logZ_quad=peak+np.log(sd*quad(lambda u:np.exp(logjoint(mu+sd*u)-peak),-12,12,epsabs=1e-12,epsrel=1e-12)[0])
    coordinates=[]
    for scale in [1e-12,1e-6,1.,1e6,1e12]:
        # A'=s*A: same prior measure, same likelihood, same physical model.
        logjoint_prime=lambda ap:logjoint(ap/scale)-np.log(scale)
        peak_prime=logjoint_prime(scale*mu)
        integral=quad(lambda u:np.exp(logjoint_prime(scale*(mu+sd*u))-peak_prime),-12,12,epsabs=1e-12,epsrel=1e-12)[0]
        logZ_prime=peak_prime+np.log(scale*sd*integral)
        coordinates.append(dict(scale=scale,log_evidence=logZ_prime,
           BIC_logjoint=np.log(len(a))-2*peak_prime,
           BIC_normalized_posterior=np.log(len(a))-2*(exact_density_peak-np.log(scale)),
           coordinate_shift=2*np.log(scale)))
    assert abs(logZ-logZ_quad)<1e-9
    assert max(abs(c['log_evidence']-logZ) for c in coordinates)<1e-9
    mle=1+b/S;chi0=float(np.sum((residual/sv)**2));chi_fit=float(np.sum(((v-mu*template)/sv)**2))
    counterexample=dict(scope='One-amplitude conditional SPARC-data statistical control; not a gravity model comparison',
        galaxy=galaxy,galaxies_in_input=len(set(r[0] for r in rows)),rows_in_input=len(rows),observations=len(a),
        fixed=dict(distance_Mpc=float(D[0]),inclination_deg=63.,disk_ML=ml,gas_scale=1.,a0_m_s2=a0,delta=1.),
        assumptions='Independent tabulated velocity errors; no inclination, distance, disk geometry or intrinsic-scatter marginalization',
        prior='A ~ Normal(1,0.1^2); Gaussian calibration amplitude, not a gravitational coupling',
        posterior_mean=mu,posterior_sd=sd,mle=mle,chi2_fixed=chi0,chi2_at_posterior_mean=chi_fit,
        log_evidence=logZ,log_evidence_quadrature=logZ_quad,log_Bayes_factor_free_vs_fixed=logBF,
        BIC_likelihood=np.log(len(a))-2*float(np.sum(norm.logpdf(v,mle*template,sv))),
        coordinates=coordinates,sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    np.savetxt(R/'results/sparc_ngc2403_control.csv',np.column_stack([radius,v,sv,template,mu*template]),delimiter=',',
               header='radius_kpc,velocity_km_s,error_km_s,fixed_template_km_s,posterior_mean_template_km_s',comments='')
    # Fixed physical ge is different from fixed Newtonian eN when delta changes.
    scan=[physical_external(d) for d in [1.,1.5,2.,3.,4.,6.,8.,12.]]
    refined=[physical_external(d,n=384) for d in [1.,2.,4.,8.,12.]]
    errors=[abs(next(x['Q2_s_2'] for x in scan if x['delta']==r['delta'])/r['Q2_s_2']-1) for r in refined]
    assert max(errors)<5e-5
    qhat=1.6; sigma=1.8;upper=qhat+norm.ppf(.975)*sigma
    root=brentq(lambda d:physical_external(d,n=128)['Q2_s_2']/1e-27-upper,1.,12.,xtol=1e-7)
    root_fine=brentq(lambda d:physical_external(d,n=256)['Q2_s_2']/1e-27-upper,root-.03,root+.03,xtol=1e-7)
    assert abs(root-root_fine)<2e-4
    qpoint=scan[0]['Q2_s_2']/1e-27
    tail=[]
    for df in [4,8,30]:
        scale=sigma*np.sqrt((df-2)/df)  # match variance, not merely scale
        tail.append(dict(df=df,scale=scale,log_likelihood_ratio_point_vs_zero=float(t.logpdf(qhat-qpoint,df,scale=scale)-t.logpdf(qhat,df,scale=scale))))
    # Algebraic identity only: not a complete perturbation/stability calculation.
    phi=np.linspace(-.1,.1,1001);c13=1-np.exp(-4*phi)
    luminal_error=float(np.max(np.abs(np.exp(-4*phi)/(1-c13)-1)))
    assert luminal_error<1e-14
    j=dict(executed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        density_coordinate_counterexample=counterexample,
        fixed_physical_external_field=dict(scan=scan,refined=refined,max_relative_refinement_error=max(errors),
            ge_m_s2=2.142e-10,a0_m_s2=1.02e-10,upper_two_sided_gaussian_95_Q2_units_1e27=upper,
            delta_crossing=root_fine,scope='Conditional theory crossing of a published Gaussian summary band; not a delta confidence bound'),
        likelihood_tail_stress=dict(point_Q2_units_1e27=qpoint,
            gaussian_log_ratio=float(norm.logpdf(qhat,qpoint,sigma)-norm.logpdf(qhat,0,sigma)),
            variance_matched_student=tail,scope='Hypothetical tail sensitivity, not evidence for Student noise in Cassini'),
        luminal_identity_max_abs_error=luminal_error)
    (R/'results/mond_diagnostics.json').write_text(json.dumps(j,indent=2)+'\n')
    print(json.dumps(j,indent=2))
if __name__=='__main__':main()
