# Relatividad general y métodos de inferencia

> **Distribución completa de software:** entradas numéricas, resultados, interfaz y cuadernos autónomos incluidos. Los manuscritos se distribuyen por separado.


[Descargar software completo](https://github.com/FacundoFirmenich/general-relativity-and-inference/archive/refs/tags/v1.11.zip) · [Cuaderno ES](https://raw.githubusercontent.com/FacundoFirmenich/general-relativity-and-inference/v1.11/GRAVITAS_ES.ipynb) · [Notebook EN](https://raw.githubusercontent.com/FacundoFirmenich/general-relativity-and-inference/v1.11/GRAVITAS_EN.ipynb)

**GRAVITAS 1.11 · Métodos reproducibles y software visual**

Facundo Firmenich · [ORCID](https://orcid.org/0009-0002-6578-3811) · f.firmenich@cedesur.org

[Castellano](README_ES.md) · [English](README_EN.md)

Plataforma visual de auditoría estadística y física: Cassini, emuladores GP, Gaia, QUMOND, DCAE y una aplicación topológica contenida de FractalGo. Incluye los módulos para procesar entradas, calcular resultados y regenerar figuras. El paper y la información suplementaria se distribuyen por separado.

## Empezar


- Windows: doble clic en `ABRIR_GRAVITAS.bat`.
- Linux: `bash ABRIR_GRAVITAS.bash`.
- macOS: abrir `ABRIR_GRAVITAS.command`.
- Cuaderno autónomo: `GRAVITAS_ES.ipynb`; también sirve en Colab desde móvil.
- Grafo sin instalación: `ui/epistemic/epistemic_graph_es.html`.

Python 3.13 recomendado. Primera instalación con Internet; servidor local en 127.0.0.1. Los lanzadores Windows/macOS no se han ejecutado en sus sistemas nativos.

## Reproducir

```bash
python run.py all
python run.py figures
python run.py tables
python run.py topology
```

Los autoejecutables `reproducir.bat`, `.bash` y `.command` preparan el entorno y ejecutan `all`. El grafo incluye seis recorridos, tour, play, pausa, búsqueda, cajas móviles y fuentes descargables. FractalGo permite cambiar prior, umbral, malla y realización y exportar componentes.

## Guías y alcance

- [Guía científica](docs/GUIDE_ES.md)
- [Hipótesis físicas](docs/PHYSICS_ES.md)
- [Información y compresión](docs/INFORMATION_ES.md)
- [FractalGo](research/fractalgo/CONTRACT_ES.md)
- [Contribuir](CONTRIBUTING_ES.md)

La convergencia numérica no prueba nueva física. Cassini utiliza épocas reconstruidas; sus discrepancias predictivas se conservan. FractalGo aquí es un diagnóstico de grafos, no una dimensión fractal. Los datos originales y claves técnicas se preservan; toda la navegación, guías e interpretación de usuario se ofrecen en castellano e inglés.

[Antecedente sismológico y alcance de la cita](docs/RELATED_WORK_ES.md) · DOI: [10.5281/zenodo.17875887](https://doi.org/10.5281/zenodo.17875887).
