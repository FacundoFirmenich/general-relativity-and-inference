"""Mode-specific statistical-volume specialization of DCAE (FKP-type).

The FKP integral is established cosmological methodology. This module does
not infer P(k), does not estimate nbar or selection functions, and does not
assert that V_eff is a directly observed physical volume. It requires their
externally evidenced, dimensioned values for each non-overlapping cell.

For cumulative domains S: phi(S,k) = V_eff(S,k)/V_phys(S).
DCAE_FKP = [1 - (phi(S,k)/phi(S0,k))^(1/3)]_+.
This is an original *reparameterization* of published FKP information weights,
not a new estimator or an independent discovery.
"""
from __future__ import annotations
import csv
import hashlib
import math
import pathlib
from collections import defaultdict
from dataclasses import dataclass,asdict
from typing import Any

from dcae import scale_asymmetry

_REQUIRED={'cell_id','domain_bin','cell_volume_cMpc3','nbar_per_cMpc3','P_k_cMpc3','k_cMpc_inv'}

@dataclass(frozen=True)
class FkpCell:
    cell_id: str
    domain_bin: int
    volume: float
    number_density: float
    power: float
    k: float

    def __post_init__(self):
        if not self.cell_id or self.domain_bin < 1 or any(
          not math.isfinite(z) or z<=0 for z in
          (self.volume,self.number_density,self.power,self.k)):
            raise ValueError('Missing/invalid physical FKP inputs: fail closed')

    @property
    def information_fraction(self):
        npow=self.number_density*self.power
        # Avoid subtracting almost equal numbers when nP is small.
        # The reciprocal branch also has the correct limit at overflow.
        ratio = npow/(1.0+npow) if npow <= 1.0 else 1.0/(1.0+1.0/npow)
        return ratio*ratio

    @property
    def v_eff(self): return self.volume*self.information_fraction


def load_cells(path:pathlib.Path)->tuple[list[FkpCell],str]:
    raw=path.read_bytes();digest=hashlib.sha256(raw).hexdigest()
    with path.open(newline='',encoding='utf-8') as f:
        reader=csv.DictReader(f)
        if not _REQUIRED.issubset(reader.fieldnames or []):
            raise ValueError('Missing columns '+str(_REQUIRED-set(reader.fieldnames or [])))
        rows=[];seen=set()
        for entry in reader:
            name=entry['cell_id']
            if name in seen: raise ValueError('Duplicate cell_id; cumulative double counting prohibited')
            seen.add(name)
            rows.append(FkpCell(name,int(entry['domain_bin']),
                float(entry['cell_volume_cMpc3']),float(entry['nbar_per_cMpc3']),
                float(entry['P_k_cMpc3']),float(entry['k_cMpc_inv'])))
    if not rows: raise ValueError('No observed cells, NOT_ESTIMABLE')
    return rows,digest


def cumulative_fkp(cells:list[FkpCell])->list[dict[str,Any]]:
    if not cells: raise ValueError('No cells')
    # All reported k values must refer to the SAME physical mode across scale.
    kv={c.k for c in cells}
    if len(kv)!=1: raise ValueError('Mixed k modes: no scalar authority')
    bins=defaultdict(list)
    for c in cells:bins[c.domain_bin].append(c)
    idx=sorted(bins)
    if idx!=list(range(1,max(idx)+1)):
        raise ValueError('Nonconsecutive domain bins; cannot assume missing shells')
    v_phys=v_eff=0.0
    results=[]
    for b in idx:
        v_phys+=math.fsum(x.volume for x in bins[b]);v_eff+=math.fsum(x.v_eff for x in bins[b])
        f=v_eff/v_phys
        if f<=0 or f>1+1e-12:raise ValueError('Information volume inconsistent')
        results.append({'domain_bin':b,'k_cMpc_inv':next(iter(kv)),
            'N_cells':sum(len(bins[k]) for k in idx if k<=b),
            'V_phys_cMpc3':v_phys,'V_eff_cMpc3':v_eff,'phi_FKP':f,
            'dcae_FKP':None, 'interpretation':'reparameterized_FKP_information_efficiency'})
    first=results[0]
    for x in results:
        z=scale_asymmetry(first['V_phys_cMpc3']**(1/3),first['V_eff_cMpc3']**(1/3),
                          x['V_phys_cMpc3']**(1/3),x['V_eff_cMpc3']**(1/3))
        x['dcae_FKP']=asdict(z)
    return results
