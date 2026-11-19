"""Reproduce original notebook findings and test corrected independent adapters."""
from pathlib import Path
import sys, json, importlib.util, warnings, hashlib
import numpy as np
from scipy.integrate import quad
R = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R))
from ui.physics import *
from research.physical_checks import gauss

def main():
    p = R/'platform/originals/GraviTools_original.py'
    spec=importlib.util.spec_from_file_location('original_gravitools',p)
    old=importlib.util.module_from_spec(spec);sys.modules[spec.name]=old;spec.loader.exec_module(old)
    alpha=1e-4; lam=.1*AU; a=.387*AU; e=.2056
    op=old.YukawaParams(alpha,lam); new=YukawaParams(alpha,lam)
    before=float(old.perihelion_precession_yukawa(a,e,old.M_SUN,op))
    after=precession_per_orbit(a,e,new)
    quadrature=alpha*gauss(e,a/lam)
    assert abs(after-quadrature)<1e-14
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        zero=old.perihelion_precession_yukawa(a,0.,old.M_SUN,op)
    rho_old=float(old.rho_eff_yukawa_r(AU,old.M_SUN,op))
    rho_new=float(poisson_density(AU,old.M_SUN,new))
    assert np.isclose(rho_old,-rho_new,rtol=1e-13)
    # Independent radial Laplacian in dimensionless coordinates.
    x=1.4; h=1e-4
    f=lambda x: -np.exp(-x)/x
    lap=(f(x+h)-2*f(x)+f(x-h))/h**2 + (f(x+h)-f(x-h))/(h*x)
    exact=-np.exp(-x)/x
    assert np.isclose(lap,exact,rtol=3e-7)
    exterior=quad(lambda x:-alpha*x*np.exp(-x),0,np.inf)[0]
    assert abs(exterior+alpha)<1e-12
    tests=[]
    for e0 in [0.,.0549,.2056,.6,.82697,.95]:
        for x0 in [.1,1.,3.,10.]:
            y=precession_per_orbit(AU,e0,YukawaParams(alpha,AU/x0))
            if e0>0:assert abs(y-alpha*gauss(e0,x0))<1e-13
            tests.append({'e':e0,'a_over_lambda':x0,'rad_per_orbit':y})
    light=[]
    for b in [R_SUN if 'R_SUN' in globals() else old.R_SUN, 10*old.R_SUN]:
        for al in [0.,-1e-4,1e-4]:
            pp=YukawaParams(al,lam); got=float(light_deflection(b,G*old.M_SUN,pp))
            q,err=light_deflection_quadrature(b,G*old.M_SUN,pp)
            assert abs(got/q-1)<1e-12
            light.append(dict(b_m=b,alpha=al,analytic=got,quadrature=q,relative_difference=got/q-1))
    # Finite-difference force control and alpha=0 Newtonian limit.
    radius=.8*AU; h=radius*1e-4
    force_fd=-(potential(radius+h,MU_SUN,new)-potential(radius-h,MU_SUN,new))/(2*h)
    force=radial_acceleration(radius,MU_SUN,new)
    assert abs(force_fd/force-1)<2e-8
    assert radial_acceleration(radius,MU_SUN,YukawaParams(0,lam))==-MU_SUN/radius**2
    findings=[
      dict(id='G01',component='precession',finding='Original expression has units s^-1 but is reported rad/orbit; true-anomaly average also lacks the time weight',original=before,corrected_rad_per_orbit=after,independent_gauss=quadrature,original_e0_finite=bool(np.isfinite(zero))),
      dict(id='G02',component='Poisson density',finding='Exterior sign reversed; Poisson effective density is not an energy-condition verdict',original_kg_m3=rho_old,corrected_kg_m3=rho_new,laplacian_relative_error=lap/exact-1,exterior_mass_over_M=exterior),
      dict(id='G03',component='optical ansatz',finding='Fixed integration endpoints induce a truncation error; replacement checked against infinite-endpoint quadrature',checks=light),
      dict(id='G04',component='multi-probe Fisher',finding='Original uses the same Mercury placeholder for every named probe; each probe needs its own observable and joint covariance',status='NOT_IMPORTED_AS_PHYSICAL_RESULT'),
      dict(id='G05',component='constraint plot',finding='Flat Cassini alpha=2.3e-5 and other curves are schematic; no likelihood-derived force constraint',status='NOT_IMPORTED_AS_OBSERVATIONAL_LIMIT'),
      dict(id='G06',component='density plot',finding='Plotted denominator is a reference Yukawa density but label claims rho_GR=0',status='REPLACED_BY_SIGNED_DENSITY'),
      dict(id='G07',component='ppE, stability, energy conditions',finding='Static potential alone does not specify radiation, perturbation action or stress-energy pressures; no universal tachyon verdict from alpha sign',status='REQUIRES_ADDITIONAL_MODEL'),
      dict(id='G08',component='Bayes factors',finding='compute_bayes_factor returns log B whereas jeffreys_scale expects B; avoid incompatible scales',status='EXPLICIT_LOG_EVIDENCE_IN_NEW_UI')]
    out=dict(source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),source_read='Complete notebook cell and saved output',scope='Audit and numerical correction of supplied code; no new observations',findings=findings,precession_checks=tests,force_derivative_relative_error=float(force_fd/force-1),checks_passed=True)
    (R/'platform/audit_runs/verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'checks_passed':True,'findings':len(findings),'original_precession_number':before,'corrected_rad_per_orbit':after,'corrected_rho':rho_new},indent=2))

if __name__=='__main__':main()
