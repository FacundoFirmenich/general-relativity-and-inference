# Relatividad general y métodos de inferencia

> **Edición de fuentes:** código, interfaz y documentación. La distribución pública de entradas numéricas y cuadernos autónomos está pendiente; el paquete completo se entrega por separado. Consultá `PUBLIC_DATA_PENDING.json`.

**GRAVITAS 1.11 · Métodos reproducibles y software visual**

Facundo Firmenich · [ORCID](https://orcid.org/0009-0002-6578-3811) · f.firmenich@cedesur.org

[Castellano](README_ES.md) · [English](README_EN.md)

Plataforma visual de auditoría estadística y física: Cassini, emuladores GP, Gaia, QUMOND, DCAE y una aplicación topológica contenida de FractalGo. Incluye los módulos para procesar entradas, calcular resultados y regenerar figuras. El paper y la información suplementaria se distribuyen por separado.

## Empezar

Los siguientes accesos requieren instalar el paquete numérico completo; esta edición de fuentes conserva su implementación.

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
