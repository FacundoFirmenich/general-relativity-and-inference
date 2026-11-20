# Guía científica y de reproducción

## Una entrada, varias comprobaciones

`python run.py all` regenera cálculos, tablas CSV, figuras y el grafo. `figures` vuelve a dibujar desde resultados computados; `calculations` recalcula. `topology`, `fermi`, `information`, `stability`, `revision` y `frontier` ejecutan bloques delimitados. Los manifiestos en `records/` registran comandos, versiones, semillas, hashes antes y después y resultados. Un código de salida cero es éxito computacional, no aceptación física.

## Casos y procedencia

- Cassini: 1696 épocas reconstruidas y 20 pasadas. Inferencia condicionada; no reducción independiente de telemetría Doppler. Las mezclas de covarianza conservan un fallo predictivo local.
- Gaia: covarianzas de registros públicos a geometría condicionada; no solución orbital global nueva.
- QUMOND/SPARC: cuadratura y controles de parametrización; no exclusión de todas las teorías MOND.
- DESI: auditoría del operador ventana; no ajuste full-shape cosmológico.
- Casimir: reproducción y auditoría de normalización; las constantes históricas se preservan como procedencia.
- Fermi-LAT: separación entre simulaciones originales, catálogo HEASARC, fotones archivados y seguimiento temporal con FT2. No es Fermilab. No proporciona una nueva cota LIV ni confirma globalmente RG.
- DCAE/COSMOS: identidad de escalas y reanálisis de tablas publicadas, no medición cosmológica independiente.
- FractalGo: topología de excursiones sobre grafos finitos y residuos temporales. No estima dimensión fractal ni introduce gravedad fractal.

## Idiomas y originales

La interfaz y el grafo tienen selector castellano/inglés; los cuadernos y guías existen en ambos idiomas. Los identificadores de los formatos JSON/CSV se conservan estables, y las tablas visibles se traducen. Los archivos originales de terceros y las salidas históricas no se reescriben ni se presentan como traducciones nuevas. Los resultados numéricos son comunes a ambas presentaciones.

## Instalación y límites comprobados

Python 3.13 es el entorno de referencia ejecutado. Los lanzadores crean un entorno aislado y descargan dependencias la primera vez. Windows y macOS tienen lanzadores preparados, pero no se afirma validación nativa en esos sistemas. Colab es la vía de navegador para iOS/iPadOS; no hay ejecutable iOS nativo. Los PDF del paper y supplementary se entregan fuera del repositorio.

## Leer un resultado

Cada pantalla enlaza datos y límites. Distinguí identidad matemática, cálculo sintético, reanálisis condicionado y validación externa. La concordancia entre algoritmos sobre los mismos datos no crea independencia. La integración canónica DCAE v0.6 explicita el protocolo y la triangulación; no incorpora como ejecutadas las nuevas afirmaciones empíricas del trabajo hermano.


[Identificación y cita de GRAVITAS](DOI_ES.md).

[DCAE: normalización y umbrales](DCAE_NORMALIZATION_ES.md).

[DCAE v0.6: definición por componentes y controles](DCAE_V06_ES.md).

[FractalGo: software independiente](https://github.com/FacundoFirmenich/fractal-go/tree/fe4fecdc44340fbadcbf1fdf3f226bad0a7661fa) · DOI: 10.5281/zenodo.24383318.
