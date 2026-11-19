"""Physical benchmarks: finite-eccentricity Yukawa precession and mass normalization.

All examples use dimensionless mu=a=1. These are dynamical benchmarks,
not fitted ephemeris constraints or a reproduction of any published exclusion.
"""
from pathlib import Path
import json,datetime
import numpy as np
from scipy.special import ive
from scipy.integrate import quad,solve_ivp
R=Path(__file__).resolve().parent

def coefficient(e,x):
    if e==0:return np.pi*x*x*np.exp(-x)
    return 2*np.pi*np.sqrt(1-e*e)*x/e*np.exp(-x*(1-e))*ive(1,e*x)

def gauss(e,x):
    def integrand(f):
        u=x*(1-e*e)/(1+e*np.cos(f))
        return (1+u)*np.exp(-u)*np.cos(f)/e
    return quad(integrand,0,2*np.pi,epsabs=2e-12,epsrel=2e-12)[0]

def orbit(e,x,alpha):
    rp=1-e;vp=np.sqrt((1+e)/(1-e))
    def rhs(t,z):
        r=np.linalg.norm(z[:2]);q=r*x
        acc=-(1+alpha*(1+q)*np.exp(-q))*z[:2]/r**3
        return np.r_[z[2:],acc]
    def peri(t,z):return np.dot(z[:2],z[2:])
    peri.direction=1
    sol=solve_ivp(rhs,[0,10*np.pi],[rp,0,0,vp],method='DOP853',rtol=2e-13,atol=2e-15,
                  max_step=.04,events=peri)
    assert sol.success
    events=sol.y_events[0];angles=np.unwrap(np.arctan2(events[:,1],events[:,0]))
    return float(np.mean(np.diff(angles))),len(angles)

def main():
    rows=[]
    for e in [.0549,.2056,.6,.82697]:
        for x in [.1,1,3,10]:
            exact=coefficient(e,x);integral=gauss(e,x)
            near=np.pi*(x*(1-e*e))**2*np.exp(-x*(1-e*e))
            assert abs(integral-exact)<2e-10
            rows.append(dict(e=e,a_over_lambda=x,first_order_exact=exact,gauss_quadrature=integral,
                             semilatus_approximation=near,approximation_relative_error=near/exact-1))
    ode=[]
    for e,x in [(.2056,1),(.82697,1),(.82697,3)]:
        exact=coefficient(e,x)
        for h in [1e-4,5e-5]:
            pos,npos=orbit(e,x,h);neg,nneg=orbit(e,x,-h);slope=(pos-neg)/(2*h)
            ode.append(dict(e=e,a_over_lambda=x,step=h,coefficient=exact,orbit_derivative=slope,
                            relative_error=slope/exact-1,events_positive=npos,events_negative=nneg))
            assert abs(slope/exact-1)<1e-5
    # Exact dynamical normalization of a circular Newtonian-order orbit.
    alpha=1e-4;x=1.;b=alpha*(1+x)*np.exp(-x);c=alpha*np.exp(-x)
    mass=dict(alpha=alpha,r_cal_over_lambda=x,
              dynamical_to_bare_mass=1+b,potential_to_bare_mass=1+c,
              gamma_light_from_circular_calibration=2/(1+b)-1,
              gamma_from_potential_normalization=2/(1+c)-1,
              omitted_factor_in_leading_gamma_deviation=1+x)
    result=dict(executed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                scope='Perturbative two-body dynamics and conditional circular-orbit calibration; no observational refit',
                precession_grid=rows,orbit_checks=ode,mass_normalization=mass)
    (R/'results/physical_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
