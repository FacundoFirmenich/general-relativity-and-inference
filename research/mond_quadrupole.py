"""Recompute QUMOND external-field quadrupoles (Park et al. 2026, Eq. 9).

Specified interpolation family and external field; no Cassini refit.
The published e_N values are rounded. The supplied Gaussian Q2 summary is
not a likelihood reconstructed from spacecraft tracking records.
"""
from pathlib import Path
import json,datetime
import numpy as np
from scipy.special import roots_legendre
R=Path(__file__).resolve().parent
def quadrupole(n,eN=1.643,delta=1.,a0=1.02e-10,vmax=64.):
    z,w=roots_legendre(n)
    total=0.
    edges=[x for x in [0.,np.sqrt(eN),2.,4.,8.,16.,32.,64.] if x<vmax]+[vmax]
    for lo,hi in zip(edges[:-1],edges[1:]):
        v=((hi-lo)*z/2+(hi+lo)/2)[:,None];wv=(hi-lo)*w/2
        for l,h in [(-1.,-.99),(-.99,-.9),(-.9,0.),(0.,1.)]:
            xi=((h-l)*z/2+(h+l)/2)[None,:];wx=(h-l)*w/2
            # Positive expression avoids cancellation at the acceleration saddle.
            x=np.sqrt((v*v-eN)**2+2*eN*v*v*(1+xi))
            logbase=np.log(-np.expm1(-x**(delta/2)))
            nu_minus_one=np.expm1(-logbase/delta)
            integrand=nu_minus_one*(eN*(3*xi-5*xi**3)+v*v*(1-3*xi*xi))
            total+=float(wv@integrand@wx)
    q=1.5*total;mu=1.3271244e20
    Q=-1.5*a0**1.5/np.sqrt(mu)*q
    return dict(order=n,eN=eN,delta=delta,vmax=vmax,a0_m_s2=a0,mu_sun_m3_s2=mu,q=q,Q2_s_2=Q)
def main():
    baseline=[quadrupole(n) for n in [48,96,192,384]]
    corrected=[quadrupole(n,eN=1.659) for n in [96,192,384]]
    assert abs(baseline[-1]['Q2_s_2']/baseline[-2]['Q2_s_2']-1)<2e-5
    assert abs(corrected[-1]['Q2_s_2']/corrected[-2]['Q2_s_2']-1)<2e-5
    scan=[quadrupole(128,delta=d) for d in [1.,1.5,2.,3.,4.,6.,8.,12.]]
    refined=[quadrupole(256,delta=r['delta']) for r in scan]
    scan_error=max(abs(a['Q2_s_2']/b['Q2_s_2']-1) for a,b in zip(scan,refined))
    tail=[quadrupole(192,vmax=v) for v in [32.,64.,128.]]
    tail_error=max(abs(a['Q2_s_2']/tail[-1]['Q2_s_2']-1) for a in tail)
    assert scan_error<2e-5 and tail_error<2e-7
    j=dict(executed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
           scope='Specified QUMOND theoretical quadrupole, not raw-data fitting or full MOND model comparison',
           reference='Park et al., Phys Rev D 114 024066 (2026), equations 4b and 9, section IV.B.1',
           convention='mu is IAU nominal solar GM; eN and a0 match rounded published scenario values',
           v_integral_upper_limit=64.,baseline_convergence=baseline,surface_density_scenario=corrected,
           interpolation_shape_scan=refined,shape_scan_refinement_relative_error=scan_error,
           radial_truncation=tail,radial_truncation_relative_error=tail_error)
    (R/'results/mond_quadrupole.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j,indent=2))
if __name__=='__main__':main()
