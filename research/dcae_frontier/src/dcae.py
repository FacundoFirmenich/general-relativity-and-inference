"""DCAE v0.2: geometrically defined scale-asymmetry change.

This module performs algebra only; it does NOT assert that an effective
observational extent has been identified from an observational selection function.
Any physically meaningful application requires an explicit, reviewable contract.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import math
from typing import Any


@dataclass(frozen=True)
class ScaleContract:
    domain: str
    target_observable: str
    coordinate_frame: str
    length_unit: str
    physical_extent_definition: str
    observational_extent_definition: str
    source_identifier: str
    source_status: str
    observational_support_kind: str
    independent_validation: bool = False

    def validate(self) -> None:
        required = (self.domain, self.target_observable, self.coordinate_frame,
                    self.length_unit, self.physical_extent_definition,
                    self.observational_extent_definition, self.source_identifier,
                    self.source_status, self.observational_support_kind)
        if not all(x and isinstance(x, str) for x in required):
            raise ValueError('Incomplete contract: abstain; no implicit fields')


@dataclass(frozen=True)
class ScaleAsymmetry:
    delta: float
    increased_0_1: float
    decreased_0_1: float
    signed_minus1_1: float
    direction: str


def scale_asymmetry(physical_ref: float, observational_ref: float,
                    physical_now: float, observational_now: float) -> ScaleAsymmetry:
    """Exactly operationalizes log[(L_F/L_O)_now/(L_F/L_O)_ref].

    All four quantities must be strictly positive and expressed in the same
    physical coordinate and unit system. Comparable definitions are an input
    CONTRACT, not a condition that software can establish from four numbers.
    """
    vals = (physical_ref, observational_ref, physical_now, observational_now)
    if any(not isinstance(x, (int, float)) or not math.isfinite(x) or x <= 0 for x in vals):
        raise ValueError('Expected four finite, strictly positive lengths')
    # Separate logs are stable against product overflow/underflow.
    delta = (math.log(physical_now)-math.log(observational_now)
             -math.log(physical_ref)+math.log(observational_ref))
    plus = -math.expm1(-max(delta,0.0))
    minus = -math.expm1(min(delta,0.0))
    signed = math.tanh(delta / 2.0)
    direction = 'INCREASING' if delta > 1e-13 else ('DECREASING' if delta < -1e-13 else 'STABLE')
    return ScaleAsymmetry(delta, plus, minus, signed, direction)


def completeness_proxy(f_ref: float, f_now: float, d: int = 2) -> ScaleAsymmetry:
    """AN EXPLICIT SURROGATE: L_O = L_F*f**(1/d).

    f must be a VERIFIED fraction of the exact support population/measure
    of the particular estimator, not galaxy counts implicitly construed as
    angular completeness. Not a cosmological measurement of scale asymmetry.
    """
    if (not math.isfinite(f_ref) or not math.isfinite(f_now) or
        not 0.0 < f_ref <= 1.0 or not 0.0 < f_now <= 1.0 or
        not isinstance(d,int) or d <= 0):
        raise ValueError('Invalid fraction or dimension; no zero imputation')
    delta = math.log(f_ref/f_now)/d
    # Equivalent to scale_asymmetry for fixed support definitions.
    return scale_asymmetry(1.0, f_ref**(1.0/d), 1.0, f_now**(1.0/d))


def admissibility(contract: ScaleContract, *, has_physical_extent: bool,
                  has_observational_extent: bool, scale_matched: bool,
                  baseline_preregistered: bool) -> dict[str, Any]:
    """Non-compensatory observational claim gate; no score substitution.

    This is a proposed *research* contract; it does not grant instrument,
    operational or publication authority.
    """
    contract.validate()
    gates = {
        'G0_domain_and_observable': contract.domain in {
            'geophysics','heliophysics','astrophysics','galactic_physics',
            'cosmology','cosmic_web'},
        'G1_physical_extent_identified': bool(has_physical_extent),
        'G2_observational_extent_identified': bool(has_observational_extent),
        'G3_same_coordinate_and_definition': bool(scale_matched),
        'G4_baseline_fixed_ex_ante': bool(baseline_preregistered),
        'G5_primary_data_or_verifiable_derivation': contract.source_status in
            {'PRIMARY_VERIFIED','DERIVED_FROM_PRIMARY_VERIFIED'},
    }
    state = 'ADMIT_MEASUREMENT_LOCAL' if all(gates.values()) else 'NOT_ESTIMABLE_DCAE'
    if not contract.independent_validation and state == 'ADMIT_MEASUREMENT_LOCAL':
        state = 'RESEARCH_LOCAL_ONLY'
    return {'state':state, 'gates':gates, 'contract':asdict(contract),
            'no_authority_transfer':True}
