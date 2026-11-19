"""DCAE normalization controls / Controles de normalización DCAE.
Finite implementation checks; analytic arguments are in the bilingual documentation.
"""
from pathlib import Path
import hashlib,json,math,platform,time,datetime
from dcae import from_scales,from_logs
R=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.perf_counter();inputs={p.name:sha(p) for p in [Path(__file__),R/'dcae.py']};checks=[]
 def check(es,en,value,details):
  assert value,en
  checks.append({'es':es,'en':en,'passed':bool(value),'details':details})
 ratios=[.01,.1,.5,1.,1.01,2.,5.,10.,20.,100.,10000.]
 err=max(abs(from_scales(r,1.,1.,1.)['index']-max(0.,1.-1./r)) for r in ratios)
 check('Cancelación del logaritmo','Logarithm cancellation',err<2e-15,{'max_abs_error':err,'ratios':ratios})
 rows=[]
 for c in [.5,.9,.95,.99]:
  r=1/(1-c);v=from_scales(r,1,1,1)
  rows.append({'chosen_index':c,'relative_asymmetry':r,'signed_log_ratio':v['delta'],'index_computed':v['index'],'derivative_dI_dR':r**-2})
 check('Umbrales fijados por definición','Thresholds fixed by definition',all(abs(x['index_computed']-x['chosen_index'])<1e-15 for x in rows),rows)
 err=0.
 for b in [2.,math.e,10.]:
  for r in ratios:
   delta_b=math.log(r)/math.log(b)
   consistent=-math.expm1(-math.log(b)*max(0.,delta_b))
   err=max(err,abs(consistent-from_scales(r,1.,1.,1.)['index']))
 check('Cambio de base consistente','Consistent base conversion',err<2e-15,{'bases':[2.,math.e,10.],'max_abs_error':err,'inconsistent_base10_at_R10':-math.expm1(-1)})
 err=max(abs((1/(1-from_scales(r,1,1,1)['index']))/r-1) for r in ratios if r>1)
 check('Inversa en la rama positiva','Positive-branch inverse',err<2e-12,{'max_relative_error':err,'scope':'finite ratios <= 10000; exact invertibility does not prevent floating-point loss'})
 low=[from_scales(r,1,1,1) for r in [.1,.5,1.]]
 check('El recorte conserva cero pero pierde la contracción','Clipping retains zero but loses contraction',all(x['index']==0 for x in low) and len({x['delta'] for x in low})==3,{'signed_gaps':[x['delta'] for x in low],'indices':[x['index'] for x in low]})
 errs=[]
 for r in [2.,10.,20.,100.]:
  h=r*1e-4;f=lambda x:from_scales(x,1,1,1)['index']
  numerical=(f(r+h)-f(r-h))/(2*h);errs.append(abs(numerical/(r**-2)-1))
 check('Sensibilidad de la coordenada','Coordinate sensitivity',max(errs)<2e-7,{'max_relative_error':max(errs),'finite_difference_relative_step':1e-4})
 sat=from_logs(1000.,0.,0.,0.)
 check('Saturación numérica explícita','Explicit numerical saturation',sat['numerical_saturation'] and sat['delta']==1000. and sat['log_remaining']==-1000.,sat)
 result={'scope_es':'Identidades y controles finitos de la normalización; ninguna medición física ni probabilidad de barrera.','scope_en':'Normalization identities and finite controls; no physical measurement or barrier probability.','all_passed':True,'checks':checks,'thresholds':rows,'inputs_sha256':inputs,'python':platform.python_version(),'runtime_seconds':time.perf_counter()-start,'executed_at':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 assert all(sha(R/k)==v for k,v in inputs.items())
 (R/'results').mkdir(exist_ok=True);(R/'results/normalization_checks.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
 print(json.dumps({'checks':len(checks),'all_passed':True,'scope':'normalization, not physical evidence'}))
if __name__=='__main__':main()
