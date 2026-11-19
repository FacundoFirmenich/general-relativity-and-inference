"""Registry of inference objects and a direct h versus hbar consistency check."""
from pathlib import Path
import json,csv
from scipy.constants import h,hbar,c,electron_volt
R=Path(__file__).resolve().parents[1];OUT=R/'results'
ref_c='https://doi.org/10.1038/nature01997';ref_m='https://arxiv.org/abs/2306.07069'
rows=[]
def add(symbol,value,units,kind,interpretation,source):
    rows.append(dict(symbol=symbol,value=value,units=units,kind=kind,interpretation=interpretation,source=source))
add('gamma_minus_1',2.1e-5,'dimensionless','published_parameter_estimate','Estimated PPN propagation parameter, not a direct measurement of Newton G',ref_c)
add('uncertainty_gamma',2.3e-5,'dimensionless','published_uncertainty','The published plus-minus; not GP interpolation error or a bound on every omitted model',ref_c)
add('GPR_mass_q997',.98e-24,'eV/c^2','posterior_quantile','Upper quantile under the selected likelihood/profile/prior and mean GP emulator',ref_m+' Table I')
add('GPUR_mass_q997',1.01e-24,'eV/c^2','posterior_quantile','Upper quantile of the uncertainty-propagation construction; not a measured nonzero mass',ref_m+' Table I')
add('GPR_mean_mass',.26e-24,'eV/c^2','posterior_mean','Positive posterior mean of a nonnegative parameter is not a detection',ref_m+' Table I')
add('GPUR_mean_mass',.34e-24,'eV/c^2','posterior_mean','Conditional posterior summary, not a direct instrument reading',ref_m+' Table I')
add('GPR_sample_max',1.45e-24,'eV/c^2','sample_maximum','Largest sampled value is not the endpoint of posterior support',ref_m+' Table I')
add('GPUR_sample_max',2.53e-24,'eV/c^2','sample_maximum','Depends on sampling; finite continuous likelihood preserves prior support',ref_m+' Table I')
add('prior_upper_mass',3.62e-23,'eV/c^2','prior_design','Analysis choice bounding the prior support; not a new measured mass',ref_m+' Section III')
add('GPUR_probability',.997,'probability','posterior_probability_level','Credible quantile level; frequentist coverage requires its own calibration',ref_m)
add('other_ephemeris_parameters',402,'count','model_dimension','Appendix A states 402 other parameters, while Section II.B states 401 for INPOP21a construction; definition/count discrepancy unresolved. Parameters are profiled, not jointly marginalized by MCMC',ref_m+' Appendix A and Section II.B')
add('GP_uncertainty_chains',300,'count','computational_design','Number of GP-perturbed chains, not physical observations',ref_m+' Section III.B.2')
add('quoted_lambda_lower',122.48e13,'km','derived_conversion','Quoted conversion to Compton length must match the declared hbar convention',ref_m+' abstract and equation (3)')
mass=1.01e-24;quoted=122.48e13
reduced=hbar*c/(mass*electron_volt)/1000;ordinary=h*c/(mass*electron_volt)/1000
mass_interval=[1.005e-24,1.015e-24]
reduced_interval=[hbar*c/(v*electron_volt)/1000 for v in mass_interval[::-1]]
ordinary_interval=[h*c/(v*electron_volt)/1000 for v in mass_interval[::-1]]
assert not reduced_interval[0]<=quoted<=reduced_interval[1]
assert ordinary_interval[0]<=quoted<=ordinary_interval[1]
check=dict(stated_mass_eV_c2=mass,reported_lambda_km=quoted,
    declared_convention='lambda = hbar/(m*c)',lambda_from_hbar_km=reduced,lambda_from_h_km=ordinary,
    reported_over_reduced=quoted/reduced,ordinary_over_reduced=ordinary/reduced,
    mass_rounding_interval_eV_c2=mass_interval,reduced_length_rounding_interval_km=reduced_interval,
    ordinary_length_rounding_interval_km=ordinary_interval,
    conclusion='Quoted length incompatible with the printed reduced-Compton convention even allowing rounding; compatible with ordinary Compton h convention',
    non_claim='Does not establish which constant was used in the underlying force code or a corrected numerical mass bound')
(OUT/'compton_consistency.json').write_text(json.dumps(check,indent=2)+'\n')
with (OUT/'number_registry.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(OUT/'number_registry.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(check,indent=2))
