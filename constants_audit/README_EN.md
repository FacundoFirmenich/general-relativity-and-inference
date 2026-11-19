# SI constants audit

Sensitivity of T=0 Drude pressure with fixed material parameters; not an experimental fit or a Yukawa bound. hbar is derived from h/(2 pi); it has no final finite decimal digit.

Archived code is preserved: its rounded constants reproduce its table but are not current SI values. This execution compares both sets without modifying the original material.

c = 299792458 m/s; h = 6.62607015e-34 J s; e = 1.602176634e-19 C; k_B = 1.380649e-23 J/K.

| Gap (nm) | Relative pressure shift (%) | Relative quadrature error |
|---:|---:|---:|
| 10 | 0.00299835 | 1.974e-12 |
| 20 | -0.00128672 | 1.267e-12 |
| 40 | -0.00611386 | 7.432e-13 |
| 80 | -0.01312869 | -6.947e-13 |
| 110 | -0.01704784 | -3.194e-12 |

Reproduce: `python constants_audit/audit_si.py`

https://www.bipm.org/en/publications/si-brochure
