# FractalGo · contrato computacional

Código y documentación bilingües; manuscrito y suplemento distribuidos por separado.

Comando: `python run.py topology`.

Afirmación finita: las bandas gaussianas marginales y los resúmenes invariantes ante permutaciones no determinan conectividad sobre un grafo fijo. Control exacto de 80 nodos; GP condicional suave frente a marginales independientes, 51/101/201 nodos y 1024 sorteos; excursiones residuales Cassini condicionales, 1696 épocas reconstruidas en 20 pasadas y 512 sorteos por modelo de ruido.

Grafo GP: camino ordenado fijo, masas previas positivas, cociente posterior/previa y umbrales 0,5/1/2. Priors uniforme y exp(-x), como escenarios explícitos. Sorteos comunes restringidos a mallas anidadas. Pesos trapezoidales nodales; no equivalen a la integración exacta de segmentos del análisis de cuantiles anterior. Sin afirmación de límite continuo.

Grafo Cassini: adyacencia cronológica dentro de pasada; huecos mayores que dos cadencias medianas cortan aristas. Mismo grafo para observado y replicado. Umbrales de signo 0/1/2 desviaciones estándar del ruido normalizado. Plantilla óptica, pesos y base nuisance heredados. Mezcla de cuatro núcleos y diagnósticos exploratorios. Fracciones predictivas no son valores p uniformemente calibrados ni detecciones de fuerza. Sin nueva estación, ajuste físico Yukawa ni SBC jerárquica completa.

Conservar IDs, filas, tiempos, aristas y hashes; sin capturas/eliminación, desplazamiento de colisiones, dimensiones inventadas ni factores de likelihood independientes derivados de los mismos datos. No se estima dimensión fractal ni se valida el constructo completo.


## English

# FractalGo · computational contract

Bilingual software and documentation; manuscript and supplement distributed separately.

Command: `python run.py topology`.

Finite assertion: fixed-graph excursion connectivity is not determined by marginal Gaussian bands or permutation-invariant scalar summaries. Test the exact 80-node witness; a smooth conditional GP versus independent marginals on 51/101/201 nodes with 1024 draws; and conditional Cassini residual excursions on 1696 reconstructed epochs in 20 passes, with 512 predictive draws for each noise model.

GP graph: fixed ordered path, strictly positive cell prior masses, posterior-to-prior ratios, thresholds 0.5/1/2. Uniform and exp(-x) priors are explicit sensitivity scenarios. Same GP draws are restricted to nested grids. Node masses use trapezoidal weights; this differs from exact segment integration used for the earlier quantile analysis. No continuous-limit claim.

Cassini graph: chronological adjacency within pass; gaps exceeding twice that pass's median cadence break edges. Same graph for observed and replicated residuals. Signed thresholds 0/1/2 whitened-noise standard deviations. The inherited optical template, weights and nuisance basis remain conditioned upon. The four-kernel mixture and all new diagnostics are exploratory, not prospectively selected. Posterior predictive tail fractions are not uniformly calibrated p-values or force detections. No new station data, no physical Yukawa fit, no complete hierarchical SBC.

Retain every node ID, original row, timestamp, edge and input hash. No capture/deletion, collision displacement, invented dimensions or independent likelihood factors from derived summaries. The prototype does not estimate a fractal dimension or validate the entire FractalGo construct.
