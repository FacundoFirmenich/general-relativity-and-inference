# DCAE v0.6 · software contract

This update integrates the canonical definition supplied by the author. The full formalization and manuscript remain outside the public repository.

## Conventions

For positive comparable scales and a fixed baseline, $\mathcal R_j=(L_{O,j}/L_F)/(L_{O,j,0}/L_{F,0})$, $\Delta_j=-\ln\mathcal R_j$ and $I_j^\uparrow=[1-\mathcal R_j]_+$. The earlier $R$ is **the reciprocal** of $\mathcal R_j$. We also retain $I_j^\downarrow=[1-1/\mathcal R_j]_+$ and $T_j=\tanh(\Delta_j/2)$.

Each component retains its protocol, baseline, coordinates and measurand. Indices are not aggregated to offset losses with improvements. The protocol declares question, instrumental architecture, calibration, acquisition/selection, reconstruction/interpretation and validation. Comparability is an obligation of physical analysis; the software does not certify it.

If $L_F$ increases with nonzero derivative, $d\Delta_j/d\ln L_F=1-\beta_j$. When the protocol changes, $\beta_j$ is a **total derivative** including that change. Its sign is not presumed.

## Run

```bash
python run.py dcae
```

`research/dcae/canonical_v06.py` provides `component(lf, lo, lf0, lo0, component_id=..., protocol_id=..., reference_id=...)`. It does not alter the historical `dcae.py` core. Inputs must use comparable units and coordinates. `None` represents an unknown quantity. Supply `lo=0` only when zero access is established by the protocol.

| State | Meaning |
|---|---|
| `OK` | Positive scales; signed gap and indices retained |
| `KNOWN_ZERO_ACCESS` | Known zero current access, positive baseline; a limit, not missing data |
| `NOT_ESTIMABLE` | A scale is missing |
| `UNDEFINED_REFERENCE` | Zero baseline access; undefined ratio |
| `INVALID_INPUT` / `INVALID_METADATA` | Invalid domain or procedure identification |

Numerical overflow, underflow and saturation are flagged; they never turn positive access into physical zero. Infinite states use explicit labels to preserve valid JSON.

`triangulation_registration` counts declared labels. $n\geq3$ meets the nominal triadic requirement; it **does not** establish independence or identifiability. These depend on the target and joint law. Procedures with shared calibration are not declared independent merely because they have different names.

## Normalization and conventional band

$I_j^\uparrow=0.90$ means $\mathcal R_j=0.1$ and $\Delta_j=\ln10$. This is an identity, not a probability or physical barrier. The v0.6 band $[0.9\ln10,1.1\ln10]$ is conventional and differs from the historical interval $[2.1,2.3]$. Its mass $\log_{10}(11/9)$ has a probabilistic interpretation only for significands under a justified Benford law. It does not automatically apply to logarithmic gaps or retrospectively selected constants. See [normalization](DCAE_NORMALIZATION_EN.md).

## Executed scope

`canonical_v06_audit.py` produces `results/canonical_v06_checks.json`: 17 finite checks of equivalence, signs, boundaries, baselines, protocol change, triangulation and the Shannon balance. Historical Gaussian/Fisher and FKP controls remain as conditional examples, not v0.6 foundations. This update does not rerun new experimental figures described in the source document, nor introduce new observations, gravity bounds or independence certification.
