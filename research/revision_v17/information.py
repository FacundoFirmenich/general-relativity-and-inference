"""Decision-specific retained information; exact Gaussian controls and real designs.
All prior scales and cross-transit covariance scales are declared scenarios.
"""
from pathlib import Path
import sys,json,csv
import numpy as np
from scipy.special import ive
R=Path(__file__).resolve().parents[2];sys.path.insert(0,str(R));sys.path.insert(0,str(Path(__file__).resolve().parent))
from cassini import Model
from research.physical_checks import coefficient,gauss

def retained_gaussian(prior_variance,J):
    if prior_variance<=0 or J<0:raise ValueError('positive variance, nonnegative information required')
    return prior_variance*J/(1+prior_variance*J)

def main():
    assert retained_gaussian(1,0)==0
    assert abs(retained_gaussian(1,1)-.5)<1e-15
    # Independent Schur-complement calculation against Woodbury in existing class.
    cr=[]
    for ell in [300.,1800.,7200.,86400.]:
        m=Model(ell)
        for tau in [0.,.3,1.,3.]:
            d=1+tau*tau*m.ev;G=m.Q.T@(m.Q/d[:,None]);c=m.Q.T@(m.t/d)
            j=m.t@(m.t/d)-c@np.linalg.solve(G+np.eye(m.p)/1e4,c)
            jf=m.t@(m.t/d)-c@np.linalg.solve(G,c)
            if tau:assert abs(j-m.stats([0,np.log(tau)])[1])<1e-8
            cr.append(dict(ell_seconds=ell,s=1.,tau=tau*1e-13,J_a_with_beta_prior=float(j),J_a_profiled=float(jf),K_quadratic_conditional=retained_gaussian(9,j)))
    raw=list(csv.DictReader(open(R/'research/results/gaia_covariance_transits.csv')));gaia=[]
    for name in sorted({r['name'] for r in raw}):
        subset=[r for r in raw if r['name']==name];sig=np.array([float(r['sigma_shared_mas']) for r in subset]);J=float(np.sum(1/sig**2))
        gaia.append(dict(object=name,transits=len(sig),conditional_common_shift_J_per_mas2=J,between_transit_common_sd_mas=[0,.01,.1],sigma_common_shift_mas=[float(np.sqrt(1/J+k*k)) for k in [0,.01,.1]],conditional_K_at_prior_sd_0_1mas=[retained_gaussian(.01,J/(1+k*k*J)) for k in [0,.01,.1]],scope='Hypothetical common-shift design with declared common calibration covariance, not measured cross-transit covariance'))
    # No source-derived observational bound is changed: only sensitivity conversion.
    eccentric=[]
    for e,label in [(.2056,'Mercury-like benchmark'),(.6,'synthetic eccentric benchmark'),(.82697,'Icarus-like benchmark')]:
        for x in [.1,1.,3.,10.]:
            exact=coefficient(e,x);approx=np.pi*(x*(1-e*e))**2*np.exp(-x*(1-e*e));assert abs(exact-gauss(e,x))<2e-10
            eccentric.append(dict(benchmark=label,e=e,a_over_lambda=x,approximation_relative_error=approx/exact-1,alpha_limit_exact_over_approx=approx/exact,scope='Same external absolute precession limit and fixed mass/orbit; no published bound recomputed'))
    obj=dict(definition='K_pi,L(E)=1-R_pi,L(E)/R_pi,L(no observation), when 0<R0<infinity and no-observation action remains feasible',cassini=cr,gaia=gaia,eccentricity=eccentric,graph_xor=dict(I_target_single_bit_nats=0.,I_target_both_bits_nats=float(np.log(2)),scope='Independent fair bits, target XOR; graph adjacency alone does not assign information'),non_claims=['No universal Fisher threshold','No empirical fractal law','No physical force posterior from an optical template','No model-independent index'])
    out=R/'research/revision_v17/results';out.mkdir(exist_ok=True);(out/'information.json').write_text(json.dumps(obj,indent=2)+'\n');print(json.dumps(obj,indent=2))
if __name__=='__main__':main()
