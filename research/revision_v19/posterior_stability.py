"""Finite-grid posterior stability under bounded cost oscillation.
No statement about total variation of the laws of random posterior vectors.
"""
from pathlib import Path
import json
import numpy as np
from scipy.special import logsumexp
R=Path(__file__).resolve().parents[2]
def audit(q_reference,q_candidate,weights=None,N=2.):
 a=np.asarray(q_reference,dtype=float);b=np.asarray(q_candidate,dtype=float)
 if a.ndim!=1 or len(a)<2 or b.shape!=a.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):raise ValueError('Finite vectors on the same grid required')
 w=np.ones_like(a) if weights is None else np.asarray(weights,dtype=float)
 if w.shape!=a.shape or np.any(w<=0) or not np.all(np.isfinite(w)) or not np.isfinite(N) or N<=0:raise ValueError('Positive finite weights and normalization required')
 # Center before scaling; common additive cost is irrelevant to normalized measures.
 sa=(N/2)*(a-a[0]);sb=(N/2)*(b-b[0]);delta=sb-sa
 if not np.all(np.isfinite(delta)):raise ValueError('Effective contrast overflow')
 lp=np.log(w)-sa;lp-=logsumexp(lp);lq=np.log(w)-sb;lq-=logsumexp(lq)
 p=np.exp(lp);q=np.exp(lq);omega=float(np.ptp(delta))
 result={'cost_normalization_N':N,'omega':omega,'TV_actual':float(np.sum(np.abs(p-q))/2),'TV_bound':float(np.tanh(omega/4)),'KL_reference_to_candidate':float(np.dot(p,lp-lq)),'KL_candidate_to_reference':float(np.dot(q,lq-lp)),'KL_bound_each_direction':omega**2/8,'CDF_sup_actual':float(np.max(np.abs(np.cumsum(p)-np.cumsum(q))))}
 assert result['TV_actual']<=result['TV_bound']+2e-12
 assert max(result['KL_reference_to_candidate'],result['KL_candidate_to_reference'])<=result['KL_bound_each_direction']+2e-12
 assert result['CDF_sup_actual']<=result['TV_actual']+2e-12
 return result

def main():
 rng=np.random.default_rng(2026111819);random=[]
 for n in [2,3,17,101]:
  for scale in [.001,.1,1.,10.]:
   a=rng.normal(size=n);b=a+scale*rng.normal(size=n);w=np.exp(rng.normal(size=n));random.append(audit(a,b,w,N=3.7))
 sharp=[]
 for omega in [.01,.1,.4,1.,2.,4.]:
  r=1/(1+np.exp(omega/2));p=np.array([r,1-r]);a=-np.log(p)
  x=audit(a,a+[0,omega]);x['requested_omega']=omega;sharp.append(x)
  assert abs(x['TV_actual']-np.tanh(omega/4))<1e-13
 common=audit([.1,.4,1.1],[123.1,123.4,124.1],[.2,1.,3.],N=20)
 assert common['TV_actual']<1e-12
 # Independent counterexample: distinct deterministic vectors have mutually singular
 # Dirac laws, even when their within-grid distributions are arbitrarily close.
 near=audit([0.,0.],[0.,1e-6]);assert 0<near['TV_actual']<1e-6
 out=R/'research/revision_v19/results';out.mkdir(exist_ok=True)
 result={'scope':'Analytic finite-grid bounds and synthetic checks; no recovered INPOP cost matrix, no continuum interpolation theorem','seed':2026111819,'random_checks':random,'sharp_two_point_witnesses':sharp,'common_cost_shift':common,'different_deterministic_posterior_vectors':{'within_grid':near,'TV_between_Dirac_laws':1},'checks':{'16_nonuniform_weight_checks':len(random)==16,'sharpness_6_cases':len(sharp)==6,'common_shift_cancels':True,'random_law_vs_conditional_measure_distinguished':True}}
 (out/'posterior_stability.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'checks':result['checks'],'sharp_TV_at_omega_0_4':sharp[2]['TV_actual'],'KL_bound_at_omega_0_4':sharp[2]['KL_bound_each_direction']},indent=2))
if __name__=='__main__':main()
