"""DCAE v0.2: finite positive scales, dimensionless log ratios, explicit status."""
import math

def from_logs(log_lf, log_lo, log_lf0, log_lo0):
    values = [log_lf, log_lo, log_lf0, log_lo0]
    if any(x is None for x in values):
        return {'status': 'NOT_ESTIMABLE', 'reason': 'missing log scale'}
    if not all(math.isfinite(x) for x in values):
        return {'status': 'INVALID_INPUT', 'reason': 'nonfinite log scale'}
    delta = math.fsum([log_lf, -log_lo, -log_lf0, log_lo0])
    positive = max(0., delta)
    value = -math.expm1(-positive)
    return {'status': 'OK', 'delta': delta, 'index': value,
            'log_remaining': -positive, 'log_ratio': log_lf-log_lo,
            'numerical_saturation': value == 1.}

def from_scales(lf, lo, lf0, lo0):
    values = [lf, lo, lf0, lo0]
    if any(x is None for x in values):
        return {'status': 'NOT_ESTIMABLE', 'reason': 'missing scale'}
    if not all(math.isfinite(x) and x > 0 for x in values):
        return {'status': 'INVALID_INPUT', 'reason': 'scales must be finite and positive'}
    return from_logs(*(math.log(x) for x in values))
