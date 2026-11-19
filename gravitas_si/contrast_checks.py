"""Numerical checks of contrast-space identities, not empirical validation."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.linalg import null_space
from scipy.special import softmax
from scipy.stats import chi2
R=Path(__file__).resolve().parents[1]

def main():
    rng=np.random.default_rng(2026111803)
    n=7;N=31.;w=np.exp(rng.normal(size=n));A=rng.normal(size=(n,n));K=A@A.T/12
    q=rng.multivariate_normal(np.zeros(n),K,size=1000)
    shift=rng.normal(size=(1000,1))
    p=softmax(np.log(w)-N*q/2,axis=1);p_shift=softmax(np.log(w)-N*(q+shift)/2,axis=1)
    D=np.c_[np.eye(n-1),-np.ones(n-1)]
    posterr=float(np.max(np.abs(p-p_shift)))
    cerr=float(np.max(np.abs(D@(K+np.ones((n,n)))@D.T-D@K@D.T)))
    diagonal_change=float(np.linalg.norm(D@(np.diag(np.diag(K))-K)@D.T))
    n=19;k=4
    X=rng.normal(size=(n,k));t=rng.normal(size=n)
    A=rng.normal(size=(n,n));C=A@A.T/n+np.eye(n)
    B=rng.normal(size=(n,5));K=B@B.T/5
    L=np.linalg.cholesky(C);Xw=np.linalg.solve(L,X);tw=np.linalg.solve(L,t)
    Kw=np.linalg.solve(L,np.linalg.solve(L,K).T).T
    U=null_space(Xw.T);u=U.T@tw;KU=U.T@Kw@U
    projected=float(u@np.linalg.solve(np.eye(len(u))+KU,u))
    solved=np.linalg.solve(C+K,np.c_[t,X]);ct=solved[:,0];cx=solved[:,1:]
    direct=float(t@ct-(t@cx)@np.linalg.solve(X.T@cx,X.T@ct))
    J0=float(u@u);lower=float(J0/(1+np.linalg.eigvalsh(KU)[-1]))
    # Printed V12 values only: this does not run an orbit solver.
    widths_i=np.array([2.873e-7,3.39e-7,1.623e-6,1.413e-5,1.284e-3])
    widths_c=np.array([1.3248e-5,1.0128e-7,.4321e-6,3.698e-6,3.3474e-4])
    ratios=widths_i/widths_c
    score=.5*((421.0/373.9)**2-(-340.9/373.9)**2)
    source=R/'inputs/author_v12/original_es.tex'
    if not source.exists():source=source.with_suffix('.txt')
    out=dict(seed=2026111803,posterior_common_shift_max_error=posterr,contrast_common_shift_max_error=cerr,
      independent_marginals_contrast_covariance_change=diagonal_change,
      information_direct=direct,information_quotient=projected,information_absolute_difference=abs(direct-projected),
      information_without_added_noise=J0,information_lower_bound=lower,
      pointwise_threshold=float(chi2.ppf(.95,1)),five_range_bonferroni_threshold=float(chi2.ppf(.99,1)),
      eighty_two_marginal_bonferroni_threshold=float(chi2.ppf(1-.05/82,1)),
      v12_rounded_endpoint_width_ratios=ratios.tolist(),v12_cacus_logscore_from_rounded_values=score,
      v12_crossover_relative_sigma_difference=abs(3.11656-3.11649)/3.11649,
      v12_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
      scope='Algebraic numerical checks and arithmetic on printed source values; no orbital or event-level replay')
    assert posterr<2e-14 and cerr<2e-14 and diagonal_change>0
    assert abs(direct-projected)<1e-11 and lower<=projected<=J0
    assert abs(score-.218)<.001
    (R/'results/contrast_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
