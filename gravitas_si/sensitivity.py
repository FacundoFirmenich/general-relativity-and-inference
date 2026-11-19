"""Additional nuisance-prior and integration-domain checks, retained separately."""
from .inference import ConditionalModel,tensor_reference,summary,ROOT
import json,numpy as np
def main():
    m=ConditionalModel(); base=json.loads((ROOT/'results/hierarchical.json').read_text())
    expanded=tensor_reference(m,301,ubounds=(-2.5,3),vbounds=(-15,4))
    diff=expanded['logz']-base['posterior']['logz'];assert abs(diff)<.001
    rows=[]
    for B in [10.,1000.]:
        r=tensor_reference(m,301,B=B);r0=tensor_reference(m,301,sa=0,B=B)
        rows.append(dict(nuisance_prior_sd=B,log_B10=r['logz']-r0['logz'],**summary(r)))
    result=dict(expanded_bounds=[[-2.5,3],[-15,4]],expanded_domain_logz_change=diff,
        nuisance_prior_sensitivity=rows,scope='Numerical checks and declared sensitivity; no physical prior calibration')
    (ROOT/'results/hierarchical_sensitivity.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
