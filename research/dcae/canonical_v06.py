"""Componentwise DCAE v0.6 / DCAE v0.6 por componentes.

Finite scales are supplied by an explicit protocol; this module does not estimate
physical access or establish independence. / Las escalas dependen de un protocolo;
este módulo no estima acceso físico ni establece independencia.
"""
import math


def component(lf, lo, lf0, lo0, *, component_id, protocol_id, reference_id):
    """Return one component, preserving signed and boundary states.
    Devuelve una componente conservando signo y estados de frontera.
    None means unknown; a measured zero access is a different state.
    None significa desconocido; un acceso nulo establecido es otro estado.
    """
    metadata = dict(component_id=component_id, protocol_id=protocol_id,
                    reference_id=reference_id, canonical_version='0.6')
    if any(not isinstance(v, str) or not v.strip() for v in
           (component_id, protocol_id, reference_id)):
        return dict(metadata, status='INVALID_METADATA')
    scales = dict(lf=lf, lo=lo, lf0=lf0, lo0=lo0)
    result = dict(metadata, scales=scales)
    if any(v is None for v in scales.values()):
        return dict(result, status='NOT_ESTIMABLE')
    if any(isinstance(v, bool) or not isinstance(v, (int, float))
           or not math.isfinite(v) for v in scales.values()):
        return dict(metadata, status='INVALID_INPUT')
    if lf <= 0 or lf0 <= 0 or lo < 0 or lo0 < 0:
        return dict(result, status='INVALID_INPUT')
    if lo0 == 0:
        return dict(result, status='UNDEFINED_REFERENCE')
    if lo == 0:
        return dict(result, status='KNOWN_ZERO_ACCESS', access_ratio=0.,
                    access_ratio_state='EXACT_ZERO', delta=None,
                    delta_state='POSITIVE_INFINITY', i_up=1., i_down=0.,
                    signed_t=1., numerical_saturation=False)
    delta = math.fsum([math.log(lf), -math.log(lo),
                       -math.log(lf0), math.log(lo0)])
    up = -math.expm1(-max(0., delta))
    down = -math.expm1(-max(0., -delta))
    try:
        ratio = math.exp(-delta)
        ratio_state = 'FINITE_POSITIVE' if ratio > 0 else 'NUMERICAL_UNDERFLOW'
    except OverflowError:
        ratio, ratio_state = None, 'NUMERICAL_OVERFLOW'
    if ratio_state == 'NUMERICAL_UNDERFLOW':
        ratio = None  # Not physically zero / No es cero físico.
    return dict(result, status='OK', access_ratio=ratio,
                access_ratio_state=ratio_state, log_access_ratio=-delta,
                delta=delta, delta_state='FINITE', i_up=up, i_down=down,
                signed_t=math.tanh(delta/2),
                numerical_saturation=up == 1. or down == 1.)


def triangulation_registration(procedure_ids):
    """Count declared procedures, never certify identification or independence.
    Cuenta procedimientos declarados; no certifica identificación o independencia.
    """
    ids = list(procedure_ids)
    if any(not isinstance(v, str) or not v.strip() for v in ids):
        raise ValueError('Identificadores vacíos / Empty identifiers')
    n = len(set(ids))
    return dict(n_distinct_labels=n, duplicate_labels=len(ids)-n,
                triadic_entry_requirement_met=n >= 3,
                identifiability='NOT_ASSESSED', independence='NOT_ASSESSED',
                required_es='Mismo objetivo, ley conjunta y dependencias documentadas.',
                required_en='Same target, joint law and documented dependencies.')
