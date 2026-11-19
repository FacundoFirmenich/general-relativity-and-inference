# FractalGo · contrato computacional

Código y documentación bilingües; manuscrito y suplemento distribuidos por separado.

Comando: `python run.py topology`.

Afirmación finita: las bandas gaussianas marginales y los resúmenes invariantes ante permutaciones no determinan conectividad sobre un grafo fijo. Control exacto de 80 nodos; GP condicional suave frente a marginales independientes, 51/101/201 nodos y 1024 sorteos; excursiones residuales Cassini condicionales, 1696 épocas reconstruidas en 20 pasadas y 512 sorteos por modelo de ruido.

Grafo GP: camino ordenado fijo, masas previas positivas, cociente posterior/previa y umbrales 0,5/1/2. Priors uniforme y exp(-x), como escenarios explícitos. Sorteos comunes restringidos a mallas anidadas. Pesos trapezoidales nodales; no equivalen a la integración exacta de segmentos del análisis de cuantiles anterior. Sin afirmación de límite continuo.

Grafo Cassini: adyacencia cronológica dentro de pasada; huecos mayores que dos cadencias medianas cortan aristas. Mismo grafo para observado y replicado. Umbrales de signo 0/1/2 desviaciones estándar del ruido normalizado. Plantilla óptica, pesos y base nuisance heredados. Mezcla de cuatro núcleos y diagnósticos exploratorios. Fracciones predictivas no son valores p uniformemente calibrados ni detecciones de fuerza. Sin nueva estación, ajuste físico Yukawa ni SBC jerárquica completa.

Conservar IDs, filas, tiempos, aristas y hashes; sin capturas/eliminación, desplazamiento de colisiones, dimensiones inventadas ni factores de likelihood independientes derivados de los mismos datos. No se estima dimensión fractal ni se valida el constructo completo.
