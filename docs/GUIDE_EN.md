# Scientific and reproduction guide

## One entry point, several checks

`python run.py all` regenerates calculations, CSV tables, figures and the graph. `figures` redraws computed results; `calculations` recomputes them. `topology`, `fermi`, `information`, `stability`, `revision` and `frontier` run bounded blocks. Manifests in `records/` record commands, versions, seeds, before/after hashes and outcomes. Exit status zero means computational success, not physical acceptance.

## Cases and provenance

- Cassini: 1696 reconstructed epochs and 20 passes. Conditional inference; no independent Doppler telemetry reduction. Covariance mixtures retain a local predictive failure.
- Gaia: public-record covariance at conditioned geometry; no new global orbit solution.
- QUMOND/SPARC: quadrature and parameterization checks; no exclusion of every MOND theory.
- DESI: window-operator audit; no cosmological full-shape fit.
- Casimir: reproduction and normalization audit; historical constants are retained as provenance.
- Fermi-LAT: original simulations, HEASARC catalogue, archived photons and FT2 timing follow-up are distinguished. This is not Fermilab. It supplies neither a new LIV bound nor a global confirmation of GR.
- DCAE/COSMOS: scale identities and published-table reanalysis, not an independent cosmological measurement.
- FractalGo: excursion topology on finite graphs and temporal residuals. It estimates no fractal dimension and introduces no fractal gravity.

## Languages and originals

The interface and graph have Spanish/English selectors; notebooks and guides exist in both languages. JSON/CSV identifiers remain stable, while visible tables are translated. Original third-party files and historical outputs are not rewritten or represented as new translations. Numerical results are shared by both presentations.

## Installation and tested limits

Python 3.13 is the executed reference environment. Launchers create an isolated environment and download dependencies on first use. Windows and macOS launchers are provided, without claiming native validation on those systems. Colab provides browser access for iOS/iPadOS; there is no native iOS executable. Paper and supplementary PDFs are distributed outside the repository.

## Reading a result

Each screen links evidence and limitations. Distinguish mathematical identities, synthetic calculations, conditional reanalysis and external validation. Agreement between algorithms using the same data does not create independence. New relative observational triangulation remains for a documented future update.


[GRAVITAS identification and citation](DOI_EN.md).
