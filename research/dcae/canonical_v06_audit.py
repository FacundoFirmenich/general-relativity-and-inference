"""Finite controls of the canonical integration / Controles finitos canónicos.
No new empirical observations / Sin nuevas observaciones empíricas.
"""
from pathlib import Path
from collections import defaultdict
import datetime, hashlib, json, math, platform, time
from canonical_v06 import component, triangulation_registration
from dcae import from_scales
R = Path(__file__).resolve().parent

def H(rows, columns):
    mass = defaultdict(float)
    for values, p in rows:
        mass[tuple(values[k] for k in columns)] += p
    return -sum(p*math.log(p) for p in mass.values() if p > 0)

def cmi(rows, a, b, given=()):
    return H(rows,a+given)+H(rows,b+given)-H(rows,given)-H(rows,a+b+given)

def main():
    start=time.perf_counter()
    files=[Path(__file__),R/'canonical_v06.py',R/'dcae.py']
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    checks=[]
    def check(es,en,ok,details=None):
        assert ok,en
        checks.append(dict(es=es,en=en,passed=bool(ok),details=details))
    def f(*v,**kw):
        return component(*v,component_id='j',protocol_id='declared_Pi',reference_id='s0',**kw)
    close=lambda a,b:math.isclose(a,b,rel_tol=1e-12,abs_tol=2e-13)
    examples=[f(1.,r,1.,1.) for r in [.01,.05,.1,.5,1.,2.,10.]]
    check('Equivalencia con el núcleo histórico','Equivalence to the historical core',
          all(close(x['i_up'],from_scales(1.,r,1.,1.)['index']) for r,x in zip([.01,.05,.1,.5,1.,2.,10.],examples)))
    check('Inversa de la razón anterior','Reciprocal of the earlier ratio',
          all(close(x['access_ratio'],math.exp(-x['delta'])) for x in examples))
    check('Dos direcciones y transformación firmada','Both directions and signed transformation',
          all(close(x['i_up'],max(0.,1-x['access_ratio'])) and close(x['i_down'],max(0.,1-1/x['access_ratio'])) and close(x['signed_t'],(1-x['access_ratio'])/(1+x['access_ratio'])) for x in examples))
    check('Referencia idéntica','Identical baseline',f(7.,3.,7.,3.)['delta']==0 and f(7.,3.,7.,3.)['i_up']==0)
    a=f(11.,2.,5.,3.);b=f(11000.,2000.,5000.,3000.)
    check('Conversión consistente de unidades','Consistent unit conversion',close(a['delta'],b['delta']))
    zero=f(2.,0.,1.,1.);missing=f(2.,None,1.,1.)
    check('Cero establecido distinto de ausencia','Established zero differs from missingness',zero['status']=='KNOWN_ZERO_ACCESS' and zero['i_up']==1 and zero['delta_state']=='POSITIVE_INFINITY' and missing['status']=='NOT_ESTIMABLE')
    check('Referencia nula indefinida','Undefined zero baseline',f(2.,1.,1.,0.)['status']=='UNDEFINED_REFERENCE')
    check('Dominio y metadatos','Domain and metadata',f(-1.,1.,1.,1.)['status']=='INVALID_INPUT' and component(1,1,1,1,component_id='',protocol_id='P',reference_id='s0')['status']=='INVALID_METADATA')
    extreme=f(1e300,1e-300,1e-300,1e300)
    check('Saturación numérica no es acceso nulo','Numerical saturation is not zero access',extreme['status']=='OK' and extreme['delta_state']=='FINITE' and extreme['numerical_saturation'] and extreme['access_ratio_state']=='NUMERICAL_UNDERFLOW',extreme)
    thresholds=[dict(p=p,access_ratio=1-p,delta=-math.log1p(-p)) for p in [.9,.95]]
    check('Umbrales de normalización','Normalization thresholds',all(close(f(1,1-r['p'],1,1)['delta'],r['delta']) for r in thresholds),thresholds)
    # x=ln LF, ln LO=b*x+k*lambda(x), lambda=q*x; frozen reference x=0.
    beta,k,q,x=.5,.2,2.,.7;eps=1e-5
    gap=lambda x:f(math.exp(x),math.exp(beta*x+k*q*x),1,1)['delta']
    deriv=(gap(x+eps)-gap(x-eps))/(2*eps)
    check('Derivada total con cambio de protocolo','Total derivative under protocol change',abs(deriv-(1-beta-k*q))<1e-10,dict(computed=deriv,expected=1-beta-k*q,omitting_protocol_change=1-beta))
    check('Tres etiquetas no prueban identificabilidad','Three labels do not prove identifiability',triangulation_registration(['A','B','C'])['triadic_entry_requirement_met'] and triangulation_registration(['A','B','C'])['identifiability']=='NOT_ASSESSED' and not triangulation_registration(['A','A','B'])['triadic_entry_requirement_met'])
    # Explicit deterministic null direction, no Gaussian approximation.
    rows=[(1,0),(2,0),(3,0)];v=(0,1);g=(0,1)
    check('Tríada de rango insuficiente','Rank-deficient triad',all(sum(a*b for a,b in zip(r,v))==0 for r in rows) and sum(a*b for a,b in zip(g,v))!=0,dict(H=rows,null_vector=v,target=g))
    # X1,U independent fair bits: loss, stability, improvement.
    # Same expanded target; different observation channels.
    balances=[]
    for scenario in ["loss","stable","gain"]:
        law=[((x,u,x if scenario!="gain" else 0,0 if scenario=="loss" else (u if scenario=="stable" else 2*x+u)),.25) for x in [0,1] for u in [0,1]]
        h1=H(law,(0,));h=H(law,(0,1))-h1;j1=cmi(law,(0,),(2,));a=cmi(law,(1,),(2,),(0,));b=cmi(law,(0,1),(3,),(2,))
        phi1=j1/h1;phi2=cmi(law,(0,1),(2,3))/H(law,(0,1))
        balances.append(dict(scenario=scenario,phi1=phi1,phi2=phi2,identity_residual=phi2-phi1-(a+b-phi1*h)/(h1+h)))
    check('Balance de Shannon sin signo universal','Shannon balance without universal sign',all(abs(v['identity_residual'])<1e-14 for v in balances) and close(balances[0]['phi2'],.5) and close(balances[1]['phi2'],1) and close(balances[2]['phi2']-balances[2]['phi1'],1),balances)
    bounds=[.9*math.log(10),1.1*math.log(10)]
    mass=math.log10(bounds[1]/bounds[0])
    check('Banda convencional distinta de la histórica','Conventional band differs from historical band',bounds[0]<2.1 and bounds[1]>2.3 and close(mass,math.log10(11/9)),dict(band=bounds,conditional_Benford_mass=mass,scope='significands under an assumed Benford law; not a probability for log gaps'))
    # Interchanging reference gives signed gap minus original; boundary cases excluded.
    swapped=f(5.,3.,11.,2.)
    check('Inversión de referencia conserva el signo opuesto','Reversing the baseline reverses the sign',close(swapped['delta'],-f(11.,2.,5.,3.)['delta']))
    check('Salida JSON sin infinitos numéricos','JSON output without numeric infinities',bool(json.dumps([zero,missing,extreme],allow_nan=False)))
    assert all(hashlib.sha256((R/n).read_bytes()).hexdigest()==v for n,v in hashes.items())
    result=dict(all_passed=True,checks=checks,canonical_version='0.6',inputs_sha256=hashes,
                scope_es='Controles matemáticos finitos y estados de la interfaz; sin nuevos datos ni certificación de independencia.',
                scope_en='Finite mathematical controls and interface states; no new data or independence certification.',
                executed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),python=platform.python_version(),runtime_seconds=time.perf_counter()-start)
    (R/'results').mkdir(exist_ok=True)
    (R/'results/canonical_v06_checks.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(checks=len(checks),all_passed=True,scope='finite canonical v0.6 integration controls')))
if __name__=='__main__':main()
