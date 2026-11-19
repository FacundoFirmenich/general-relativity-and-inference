"""SI values and their distinct epistemic status; sources are explicit."""
C = 299792458.0  # m/s; exact SI definition
G = 6.67430e-11  # m^3 kg^-1 s^-2; CODATA 2022 estimate
G_STANDARD_UNCERTAINTY = 0.00015e-11
AU = 149597870700.0  # m; IAU 2012 Resolution B2, exact
MU_SUN_NOMINAL = 1.3271244e20  # m^3/s^2; exact IAU 2015 B3 conversion, not measured Sun GM
REGISTRY=[
 {'symbol':'c','value_si':C,'standard_uncertainty_si':0.0,'unit':'m s^-1','status_es':'Exacta por definición SI','status_en':'Exact by SI definition','source':'https://www.bipm.org/en/measurement-units/'},
 {'symbol':'G','value_si':G,'standard_uncertainty_si':G_STANDARD_UNCERTAINTY,'unit':'m^3 kg^-1 s^-2','status_es':'Estimación CODATA 2022','status_en':'CODATA 2022 estimate','source':'https://doi.org/10.1103/RevModPhys.97.025002'},
 {'symbol':'au','value_si':AU,'standard_uncertainty_si':0.0,'unit':'m','status_es':'Unidad astronómica exacta IAU 2012 B2','status_en':'Exact astronomical unit, IAU 2012 B2','source':'https://www.iau.org/static/resolutions/IAU2012_English.pdf'},
 {'symbol':'(GM)_sun_nominal','value_si':MU_SUN_NOMINAL,'standard_uncertainty_si':0.0,'unit':'m^3 s^-2','status_es':'Conversión nominal exacta IAU 2015 B3; no estimación solar','status_en':'Exact nominal IAU 2015 B3 conversion; not a solar estimate','source':'https://www.iau.org/static/resolutions/IAU2015_English.pdf'}]
