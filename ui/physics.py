"""Audited adapters for the phenomenological additive Yukawa potential.

These functions do not specify a covariant action, radiative sector, screening
mechanism or stress-energy tensor. Such claims require additional physics.
"""
from dataclasses import dataclass
import numpy as np
from scipy.special import k1e
from scipy.integrate import quad
from research.physical_checks import coefficient

from ui.constants import G, C, AU, MU_SUN_NOMINAL
MU_SUN = MU_SUN_NOMINAL

@dataclass(frozen=True)
class YukawaParams:
    alpha: float
    lambda_m: float

    def __post_init__(self):
        if not np.isfinite(self.alpha) or not np.isfinite(self.lambda_m) or self.lambda_m <= 0:
            raise ValueError('Finite alpha and strictly positive finite lambda required')

def _positive(value, name):
    value = np.asarray(value, dtype=float)
    if not np.isfinite(value).all() or (value <= 0).any():
        raise ValueError(name + ' must be positive and finite')
    return value

def potential(r, mu, params):
    r = _positive(r, 'radius'); _positive(mu, 'GM')
    return -mu / r * (1 + params.alpha * np.exp(-r/params.lambda_m))

def radial_acceleration(r, mu, params):
    r = _positive(r, 'radius'); _positive(mu, 'GM')
    x = r / params.lambda_m
    return -mu / r**2 * (1 + params.alpha * (1+x) * np.exp(-x))

def precession_per_orbit(a, e, params):
    """First order in alpha; finite-eccentricity Kepler-orbit coefficient.

    e=0 is the epicyclic limit: a geometrical periapsis is undefined there.
    Validity requires a perturbatively small force correction on the orbit.
    """
    _positive(a, 'semimajor axis')
    if not np.isfinite(e) or not 0 <= e < 1:
        raise ValueError('0 <= e < 1 required')
    return float(params.alpha * coefficient(float(e), float(a/params.lambda_m)))

def poisson_density(r, mass, params):
    """Newtonian effective mass density for r>0, not a GR energy-condition test.

    Full distribution: M(1+alpha)delta^3(r) plus the exterior term below.
    The exterior integral is -M alpha; total distributional mass remains M.
    """
    r = _positive(r, 'radius'); _positive(mass, 'mass')
    return -mass*params.alpha*np.exp(-r/params.lambda_m)/(4*np.pi*r*params.lambda_m**2)

def light_deflection(b, mu, params, gamma_light=1.):
    """Straight-line, infinite-endpoint weak-field ansatz Psi=gamma_light Phi.

    This optional optical ansatz is independent of the force template; it is
    not a universal conversion from alpha to PPN gamma or to Cassini bounds.
    """
    b = _positive(b, 'impact parameter'); _positive(mu, 'GM')
    if not np.isfinite(gamma_light):
        raise ValueError('finite gamma_light required')
    x = b/params.lambda_m
    response = x * k1e(x) * np.exp(-x)
    return 2*(1+gamma_light)*mu/(C*C*b)*(1+params.alpha*response)

def light_deflection_quadrature(b, mu, params, gamma_light=1.):
    """Independent angular substitution z=b tan(t), no arbitrary z truncation."""
    _positive(b, 'impact parameter'); _positive(mu, 'GM')
    x = b/params.lambda_m
    f = lambda t: np.cos(t)*(1+params.alpha*(1+x/np.cos(t))*np.exp(-x/np.cos(t)))
    value, error = quad(f, -np.pi/2, np.pi/2, epsabs=2e-13, epsrel=2e-13)
    return (1+gamma_light)*mu/(C*C*b)*value, error
