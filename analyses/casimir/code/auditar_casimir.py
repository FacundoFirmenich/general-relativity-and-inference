"""Independent T=0 Drude Lifshitz quadrature and audit of public numerical table.

No experimental force or Yukawa bound is inferred from the theoretical table.
Constants intentionally match the archived generation code for a fair comparison.
"""
from pathlib import Path
import json,hashlib,csv
import numpy as np
import h5py
from scipy.special import roots_legendre
from scipy.optimize import curve_fit,brentq
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'sources/casimir/selected'; OUT=ROOT/'results'
HBAR=1.055e-34; C=2.9979e8; EV=1.602e-19; KB=1.381e-23
WP=13*EV/HBAR; GAMMA=.05*EV/HBAR; EPSR=1.03

def pressure(a,ny=384,nu=256,cutoff=None,ideal=False):
    """Integrate y=2aq, u=xi/(cq). Optional finite imaginary-frequency cutoff."""
    xu,wu=roots_legendre(nu); xu=(xu+1)/2;wu=wu/2
    cuts=[0.,60.]
    if cutoff is not None:
        ycut=2*a*cutoff/C
        if 0<ycut<60:cuts=[0.,ycut,60.]
    total=0.
    for lo,hi in zip(cuts[:-1],cuts[1:]):
        x,w=roots_legendre(ny);y=(lo+hi)/2+(hi-lo)*x/2;w=w*(hi-lo)/2
        umax=np.ones(len(y)) if cutoff is None else np.minimum(1,2*a*cutoff/(C*y))
        u=umax[:,None]*xu
        if ideal:
            z=2/np.expm1(y)[:,None]*np.ones_like(u)
        else:
            xi=C*y[:,None]*u/(2*a)
            eps=EPSR+WP**2/(xi*(xi+GAMMA))
            sq=np.sqrt(1+(eps-1)*u*u)
            te=((1-sq)/(1+sq))**2;tm=((eps-sq)/(eps+sq))**2
            decay=np.exp(-y)[:,None]
            z=te*decay/(1-te*decay)+tm*decay/(1-tm*decay)
        total+=np.dot(w*y**3,(z@wu)*umax)
    return -HBAR*C/(32*np.pi**2*a**4)*total

def finite_T_archived_cutoffs(a,nodes=128):
    """Independent Matsubara sum at 0.25 K with the archived exact truncations.

    Uses q, not transverse momentum, for numerical integration; the n=0
    regularization n=.005 and its half weight match the archived code.
    """
    x,w=roots_legendre(nodes);total=0.;T=.25
    for first in range(0,200001,2048):
        n=np.arange(first,min(first+2048,200001),dtype=float)
        if first==0:n[0]=.005
        xi=2*np.pi*KB*T*n/HBAR
        eps=EPSR+WP**2/(xi*(xi+GAMMA))
        low=2*a*xi/C
        high=np.minimum(np.sqrt(low*low+(2*a*1e9)**2),low+50)
        y=(low+high)[:,None]/2+(high-low)[:,None]*x/2
        u=low[:,None]/y
        sq=np.sqrt(1+(eps[:,None]-1)*u*u)
        te=((1-sq)/(1+sq))**2
        tm=((eps[:,None]-sq)/(eps[:,None]+sq))**2
        decay=np.exp(-y)
        terms=y*y*(te*decay/(1-te*decay)+tm*decay/(1-tm*decay))
        val=(terms@w)*(high-low)/2
        if first==0:val[0]*=.5
        total+=val.sum()
    return -KB*T/(8*np.pi*a**3)*total

def main():
    f=SRC/'Figure 5/Casimir_pressure_dist_v4.h5'
    with h5py.File(f) as h:
        a=h['data/block0_values'][...]
        labels=[v.decode() for v in h['data/block0_items'][...]]
    rows=[]; conv=[]; idealchecks=[]
    cutoff=2*np.pi*KB*.25*200000/HBAR
    for row in a:
        nm,praw,pbcs,pplasma=row
        dist=nm*1e-9
        p=pressure(dist)
        coarse=pressure(dist,ny=192,nu=128)
        if nm in (10,20,40,80,110):
            fine=pressure(dist,ny=768,nu=512)
            conv.append(dict(gap_nm=nm,p_384_256=p,p_768_512=fine,relative_difference=abs(p/fine-1)))
        idealcalc=pressure(dist,ideal=True)
        idealexact=-np.pi**2*HBAR*C/(240*dist**4)
        idealchecks.append(abs(idealcalc/idealexact-1))
        pc=pressure(dist,cutoff=cutoff)
        rows.append(dict(gap_nm=nm,table_Drude_Pa=praw,table_BCS_Pa=pbcs,
            independent_T0_Drude_Pa=p,independent_T0_frequency_cutoff_Pa=pc,
            pi_times_table_Drude_Pa=np.pi*praw,independent_over_table=p/praw,
            ratio_independent_over_pi_table=p/(np.pi*praw),
            ratio_frequency_cutoff_over_pi_table=pc/(np.pi*praw),
            relative_frequency_tail=abs((p-pc)/p),relative_coarse_difference=abs(p/coarse-1)))
    def power(x,amp,n):return amp/x**n
    pars,cov=curve_fit(power,a[:,0],np.pi*a[:,1],[-7e7,3.2],maxfev=20000)
    pred=power(a[:,0],*pars);res=np.pi*a[:,1]-pred
    r2=1-np.dot(res,res)/np.sum((np.pi*a[:,1]-np.mean(np.pi*a[:,1]))**2)
    rawpars,rawcov=curve_fit(power,a[:,0],a[:,1],[-7e7/np.pi,3.2],maxfev=20000)
    rawres=a[:,1]-power(a[:,0],*rawpars)
    rawr2=1-np.dot(rawres,rawres)/np.sum((a[:,1]-np.mean(a[:,1]))**2)
    # Propagate only the force normalization into the recovered static resonator
    # equations, keeping archived nuisance parameters fixed. This is not a refit.
    init=np.loadtxt(SRC/'Figure 1/Optimized_mech_freqs_v2.txt')
    gap0=18.;wm=np.interp(gap0,init[0],init[1]);mass=3.96e-14
    area=np.pi*(.5*11.3e-6)**2;exponent=3.19273178;static=[]
    for factor in (1.,1/np.pi):
        coeff=6.95778690e7*factor
        balance=lambda gap:mass*wm**2*(gap-gap0)*1e-9+area*coeff/gap**exponent
        gap=brentq(balance,exponent*gap0/(exponent+1),gap0,xtol=1e-12)
        stiffness=mass*wm**2-exponent*area*coeff/gap**(exponent+1)*1e9
        static.append(dict(coefficient_factor=factor,equilibrium_gap_nm=gap,
            frequency_Hz=np.sqrt(stiffness/mass)/(2*np.pi),force_balance_residual_N=balance(gap)))
    arrays=[]
    for af in SRC.rglob('*.npy'):
        z=np.load(af)
        arrays.append(dict(file=str(af.relative_to(ROOT)),shape=list(z.shape),finite=int(np.isfinite(z).sum()),
            total=int(z.size),sha256=hashlib.sha256(af.read_bytes()).hexdigest()))
    # At the minimum quasiparticle energy, quantify the unit mismatch in the
    # archived BCS tanh. This is a pointwise factor check, not a force-error bound.
    gap=1.764*KB*1.2*np.sqrt(1-.25/1.2)*(.9963+.7735*.25/1.2)/EV
    arg_correct=gap*EV/(2*KB*.25);arg_archived=gap/(2*KB*.25)
    finite=[]
    for nm in (10.,40.,110.):
        pa=finite_T_archived_cutoffs(nm*1e-9,96)
        pb=finite_T_archived_cutoffs(nm*1e-9,192)
        stored=a[np.argmin(abs(a[:,0]-nm)),1]
        finite.append(dict(gap_nm=nm,pressure_Pa=pb,table_pressure_Pa=stored,
            relative_quadrature_difference=abs(pa/pb-1),relative_difference_to_table=pb/stored-1))
    output=dict(scope='Independent Drude T=0 and finite-T theory audit, not a new experimental force fit',
        original_table_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),column_labels=labels,
        constants=dict(hbar=HBAR,c=C,eV=EV,kB=KB,Omega_p_eV=13,gamma_eV=.05,epsilon_r=EPSR),
        perfect_conductor_max_relative_error=max(idealchecks),convergence=conv,
        finite_T_025K_same_cutoffs=finite,
        power_fit_pi_times_Drude_table=dict(amplitude=pars[0],exponent=pars[1],R_squared=r2,
            meaning='Fit to numerical theory samples, not independent measurement uncertainty'),
        power_fit_raw_Drude_table=dict(amplitude=rawpars[0],exponent=rawpars[1],R_squared=rawr2),
        resonator_fixed_parameter_propagation=dict(unperturbed_gap_nm=gap0,
            unperturbed_frequency_Hz=wm/(2*np.pi),effective_mass_kg=mass,area_m2=area,
            results=static,scope='Static equilibrium and small-amplitude frequency at fixed archived parameters; neither a fit to measured sidebands nor a corrected experimental gap'),
        archived_matsubara_cutoff_energy_eV=cutoff*HBAR/EV,
        BCS_unit_audit=dict(energy_gap_eV=gap,correct_tanh_argument_at_gap=arg_correct,
            archived_tanh_argument_at_gap=arg_archived,correct_tanh_at_gap=np.tanh(arg_correct),
            caution='Dimensional inconsistency in archived expression; force impact not computed; current audit covers Drude only'),
        experimental_array_inventory=arrays,rows=rows)
    (OUT/'casimir_audit.json').write_text(json.dumps(output,indent=2)+'\n')
    with (OUT/'casimir_pressure_comparison.csv').open('w') as fp:
        wr=csv.DictWriter(fp,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
    assert max(idealchecks)<1e-10
    assert max(c['relative_difference'] for c in conv)<1e-5
    print(json.dumps({k:v for k,v in output.items() if k!='rows'},indent=2))
    print('Selected pressure comparisons:')
    print(json.dumps([r for r in rows if r['gap_nm'] in (10,20,40,80,110)],indent=2))

if __name__=='__main__':main()
