"""Source-motivated controls. All numerical examples are synthetic, not INPOP chains."""
from pathlib import Path
import json, datetime
import numpy as np
from scipy.signal import lfilter

R=Path(__file__).resolve().parent

def main():
    rng=np.random.default_rng(2026111803)
    rho=.9; n=15000; chains=5; repetitions=500
    lag=np.arange(1,n,dtype=float)
    inflation=float(1+2*np.sum((1-lag/n)*rho**lag))
    means=[]; example=None
    for _ in range(repetitions):
        innovation=rng.normal(size=(chains,n))*np.sqrt(1-rho*rho)
        innovation[:,0]=rng.normal(size=chains) # exact stationary initial law
        sample=lfilter([1],[1,-rho],innovation,axis=1)
        means.append(sample.mean())
        if example is None: example=sample
    flat=example.ravel(); shuffled=rng.permutation(flat)
    autocorr=lambda x: float(np.corrcoef(x[:-1],x[1:])[0,1])
    mcmc=dict(rho=rho,chains=chains,retained_per_chain=n,repetitions=repetitions,
              total_states=chains*n,exact_variance_inflation=inflation,
              exact_ess_for_mean=chains*n/inflation,
              exact_mcse=np.sqrt(inflation/(chains*n)),
              incorrect_iid_mcse=1/np.sqrt(chains*n),
              empirical_mcse=float(np.std(means,ddof=1)),
              pooled_lag1=autocorr(flat),shuffled_lag1=autocorr(shuffled),
              unchanged_mean_error=float(abs(flat.mean()-shuffled.mean())))
    assert mcmc['unchanged_mean_error']<1e-15
    # A shared systematic mode and N-inflated independent errors preserve
    # information about a common mean, but not within-transit contrasts.
    nccd=9; random_var=1.; systematic_var=.25
    C=random_var*np.eye(nccd)+systematic_var*np.ones((nccd,nccd))
    D=(random_var+nccd*systematic_var)*np.eye(nccd)
    common=np.ones(nccd); contrast=np.arange(nccd,dtype=float);contrast-=contrast.mean()
    contrast/=np.linalg.norm(contrast)
    info=lambda v,C:float(v@np.linalg.solve(C,v))
    covariance=dict(n_ccd=nccd,random_variance=random_var,systematic_variance=systematic_var,
                    common_information_block=info(common,C),common_information_inflated=info(common,D),
                    contrast_information_block=info(contrast,C),contrast_information_inflated=info(contrast,D))
    assert abs(covariance['common_information_block']-covariance['common_information_inflated'])<1e-12
    # Group-wise postfit scale expectation: trace term plus squared omitted signal.
    N=40; t=np.linspace(-1,1,N); X=np.c_[np.ones(N),t]
    C=.6**np.abs(np.subtract.outer(np.arange(N),np.arange(N)))
    Ci=np.linalg.inv(C); H=X@np.linalg.solve(X.T@Ci@X,X.T@Ci); M=np.eye(N)-H
    W=np.diag(np.r_[np.ones(12),np.zeros(N-12)])
    s=t*t; alpha=1.7
    expected_noise=float(np.trace(W@M@C@M.T))
    signal_term=float(alpha**2*s@M.T@W@M@s)
    samples=rng.multivariate_normal(alpha*s,C,size=20000)@M.T
    q=np.einsum('bi,ij,bj->b',samples,W,samples)
    weighting=dict(group_size=12,fitted_parameters=2,naive_df=10,
                   trace_noise=expected_noise,omitted_signal_term=signal_term,
                   expected_group_quadratic=expected_noise+signal_term,
                   simulated_group_quadratic=float(q.mean()),simulation_mcse=float(q.std(ddof=1)/np.sqrt(len(q))))
    assert abs(q.mean()-expected_noise-signal_term)<5*weighting['simulation_mcse']
    result=dict(executed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),seed=2026111803,
                scope='Synthetic controls with known probability laws; no external-chain reanalysis',
                mcmc_pooling=mcmc,transit_covariance=covariance,group_reweighting=weighting)
    (R/'results/inference_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
