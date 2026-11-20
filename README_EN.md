# General Relativity and Inference Methods

> **Complete software distribution:** numerical inputs, results, interface and self-contained notebooks are included. Manuscripts are distributed separately.


[Download complete software](https://github.com/FacundoFirmenich/general-relativity-and-inference/archive/refs/tags/v1.11.zip) · [Cuaderno ES](https://raw.githubusercontent.com/FacundoFirmenich/general-relativity-and-inference/v1.11/GRAVITAS_ES.ipynb) · [Notebook EN](https://raw.githubusercontent.com/FacundoFirmenich/general-relativity-and-inference/v1.11/GRAVITAS_EN.ipynb)

**GRAVITAS 1.11 · Reproducible methods and visual software**

Facundo Firmenich · [ORCID](https://orcid.org/0009-0002-6578-3811) · f.firmenich@cedesur.org

[Castellano](README_ES.md) · [English](README_EN.md)

A visual platform for statistical and physical auditing: Cassini, GP emulators, Gaia, QUMOND, DCAE and a bounded topological application of FractalGo. Includes modules for processing inputs, computing results and regenerating figures. The paper and supplementary information are distributed separately.

## Start


- Windows: double-click `OPEN_GRAVITAS.bat`.
- Linux: `bash OPEN_GRAVITAS.bash`.
- macOS: open `OPEN_GRAVITAS.command`.
- Self-contained notebook: `GRAVITAS_EN.ipynb`; also usable in Colab from mobile.
- Graph without installation: [interactive graph EN](ui/epistemic/epistemic_graph_en.html).

Python 3.13 recommended. First installation needs Internet; the local server binds to 127.0.0.1. Windows/macOS launchers have not been executed on their native systems.

## Reproduce

```bash
python run.py all
python run.py figures
python run.py tables
python run.py topology
```

The `reproduce.bat`, `.bash` and `.command` launchers prepare the environment and execute `all`. The graph includes six routes, tour, play, pause, search, movable cards and downloadable sources. FractalGo provides prior, threshold, grid and draw controls and component export.

## Guides and scope

- [Scientific guide](docs/GUIDE_EN.md)
- [Physical assumptions](docs/PHYSICS_EN.md)
- [Information and compression](docs/INFORMATION_EN.md)
- [FractalGo](research/fractalgo/CONTRACT_EN.md)
- [Contributing](CONTRIBUTING_EN.md)

GRAVITAS connects mathematical verification, instrumental reconstruction and predictive diagnosis. Cassini uses reconstructed epochs and localizes pass-level discrepancies; FractalGo examines connectivity on declared graphs. Original data and technical keys are preserved; navigation, guides and user interpretation are provided in Spanish and English.


[GRAVITAS identification and citation](docs/DOI_EN.md) · DOI: [10.5281/zenodo.24381484](https://doi.org/10.5281/zenodo.24381484).

Post-v1.11 update: [canonical DCAE v0.6](docs/DCAE_V06_EN.md), with executable componentwise controls. The archived version tag and DOI have not been changed.

Independent FractalGo software: [version 0.1.0](https://github.com/FacundoFirmenich/fractal-go/tree/fe4fecdc44340fbadcbf1fdf3f226bad0a7661fa); DOI: 10.5281/zenodo.24383318.

FractalGo v0.1.0: [10.5281/zenodo.24383318](https://doi.org/10.5281/zenodo.24383318) · All versions: [10.5281/zenodo.24383317](https://doi.org/10.5281/zenodo.24383317).
