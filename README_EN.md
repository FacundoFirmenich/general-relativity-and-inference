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
- Graph without installation: `ui/epistemic/epistemic_graph_en.html`.

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

Numerical convergence does not prove new physics. Cassini uses reconstructed epochs; predictive discrepancies are retained. FractalGo here is a graph diagnostic, not a fractal dimension. Original data and technical keys are preserved; navigation, guides and user interpretation are provided in Spanish and English.
