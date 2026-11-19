"""Weak-field benchmarks. This module is not a planetary ephemeris integrator."""
import numpy as np
from scipy.constants import G, c

def _check(r, length):
    if np.any(np.asarray(r) <= 0) or length <= 0:
        raise ValueError('r and Yukawa length must be strictly positive')

def potential(r, mu, alpha=0., length=1.):
    """Specific potential in SI: mu=G M (m^3/s^2), r,length in metres."""
    _check(r, length)
    return -mu / np.asarray(r) * (1 + alpha * np.exp(-np.asarray(r)/length))

def radial_acceleration(r, mu, alpha=0., length=1.):
    _check(r, length)
    x=np.asarray(r)/length
    return -mu/np.asarray(r)**2 * (1+alpha*(1+x)*np.exp(-x))

def acceleration_1pn(position, velocity, mu):
    """Test-mass, central nonrotating body, harmonic coordinates, through 1PN.
    No N-body, spin, multipoles, radiation reaction or clock model is implied.
    """
    r=np.asarray(position, dtype=float); v=np.asarray(velocity,dtype=float)
    d=np.linalg.norm(r)
    if d <= 0: raise ValueError('Nonzero separation required')
    return -mu*r/d**3 + mu/(c*c*d**3)*((4*mu/d-v@v)*r+4*(r@v)*v)

def shapiro_one_way(r1,r2,separation,mu,gamma=1.):
    """Leading PPN static monopole light time, seconds. Coordinate distances SI."""
    if min(r1,r2,separation)<=0 or separation>=r1+r2:
        raise ValueError('Invalid light-path geometry')
    return (1+gamma)*mu/c**3*np.log1p(2*separation/(r1+r2-separation))

def validate():
    r=np.geomspace(.1,20,100); alpha=.2; lam=2.; mu=1.
    assert np.all(radial_acceleration(r,mu,0.,lam)==-mu/r**2)
    errs=[]
    # Centered differences, then Richardson extrapolation (not a copied force formula).
    for h in [.04,.02,.01,.005]:
        d1=-(potential(r+h,mu,alpha,lam)-potential(r-h,mu,alpha,lam))/(2*h)
        d2=-(potential(r+h/2,mu,alpha,lam)-potential(r-h/2,mu,alpha,lam))/h
        exact=radial_acceleration(r,mu,alpha,lam)
        errs.append(float(np.max(abs(((4*d2-d1)/3-exact)/exact))))
    assert errs[-1]<2e-6 and errs[-2]/errs[-1]>15
    # Original v3.1 formula: Newtonian recovery fails at alpha=0, r=lambda.
    old_ratio=2/np.e
    # The contact Yukawa multiplier, including an L2 operator bound.
    k=np.geomspace(.01,100,300); ls=np.geomspace(.001,.1,30)
    bound_errors=[]
    for l in ls:
        remainder=1/(1+(l*k)**2)-(1-(l*k)**2)
        bound_errors.append(float(np.max(remainder-(l*k)**4)))
    assert max(bound_errors)<1e-12
    from scipy.integrate import solve_ivp
    e=.20563; eps=1e-6
    def rhs(t,y):
        rr=y[:2]; vv=y[2:]; unit_length=1e9
        aa=acceleration_1pn(rr*unit_length,vv*c*np.sqrt(eps),c*c*eps*unit_length)*unit_length/(c*c*eps)
        return np.r_[vv,aa]
    def peri(t,y):return np.dot(y[:2],y[2:])
    peri.direction=1
    sol=solve_ivp(rhs,[0,8*2*np.pi],[1-e,0,0,np.sqrt((1+e)/(1-e))],method='DOP853',rtol=2e-12,atol=2e-13,events=peri,max_step=.1)
    assert sol.success
    angles=np.unwrap(np.arctan2(sol.y_events[0][:,1],sol.y_events[0][:,0]))
    observed=float(np.mean(np.diff(angles)));predicted=6*np.pi*eps/(1-e*e)
    assert abs(observed/predicted-1)<2e-4
    from scipy.integrate import quad
    b=1e9;left=-2e11;right=3e11;mu_s=1.3271244e20
    direct=2*mu_s/c**3*quad(lambda x:1/np.sqrt(b*b+x*x),left,right,epsabs=1e-11,epsrel=1e-12)[0]
    closed=shapiro_one_way(np.hypot(b,left),np.hypot(b,right),right-left,mu_s)
    assert abs(direct/closed-1)<1e-10
    return dict(richardson_relative_errors=errs,original_force_alpha_zero_ratio_at_r_equal_lambda=old_ratio,
        original_relative_error=old_ratio-1,contact_multiplier_bound_max_excess=max(bound_errors),
        periapsis_1pn=dict(eccentricity=e,dimensionless_mu_over_ac2=eps,integrated_radians_per_orbit=observed,leading_prediction=predicted,relative_error=observed/predicted-1),
        shapiro_quadrature_relative_error=direct/closed-1,
        scope='Analytic weak-field benchmarks; no claim of full GR ephemeris reproduction')
