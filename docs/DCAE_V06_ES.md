# DCAE v0.6 · contrato del software

Esta actualización integra la definición canónica aportada por el autor. La formalización completa y el manuscrito permanecen fuera del repositorio público.

## Convenciones

Para escalas comparables positivas y referencia fija, $\mathcal R_j=(L_{O,j}/L_F)/(L_{O,j,0}/L_{F,0})$, $\Delta_j=-\ln\mathcal R_j$ e $I_j^\uparrow=[1-\mathcal R_j]_+$. El antiguo $R$ es **el inverso** de $\mathcal R_j$. Se conservan además $I_j^\downarrow=[1-1/\mathcal R_j]_+$ y $T_j=\tanh(\Delta_j/2)$.

Cada componente conserva su protocolo, referencia, coordenadas y mensurando. Los índices no se agregan para compensar pérdidas con mejoras. El protocolo declara pregunta, arquitectura instrumental, calibración, adquisición/selección, reconstrucción/interpretación y validación. La comparabilidad es una obligación del análisis físico; el software no la certifica.

Si $L_F$ crece con derivada no nula, $d\Delta_j/d\ln L_F=1-\beta_j$. Cuando cambia el protocolo, $\beta_j$ es una **derivada total**, que incluye ese cambio. El signo no se presume.

## Ejecutar

```bash
python run.py dcae
```

`research/dcae/canonical_v06.py` ofrece `component(lf, lo, lf0, lo0, component_id=..., protocol_id=..., reference_id=...)`. No altera el núcleo histórico `dcae.py`. La entrada debe usar unidades y coordenadas comparables. `None` representa una cantidad desconocida. `lo=0` sólo se introduce cuando el acceso nulo está establecido por el protocolo.

| Estado | Lectura |
|---|---|
| `OK` | Escalas positivas; se conservan desfase firmado e indicadores |
| `KNOWN_ZERO_ACCESS` | Acceso actual conocido nulo, referencia positiva; límite, no ausencia de datos |
| `NOT_ESTIMABLE` | Falta una escala |
| `UNDEFINED_REFERENCE` | Acceso de referencia nulo; cociente no definido |
| `INVALID_INPUT` / `INVALID_METADATA` | Dominio o identificación del procedimiento incorrectos |

Los desbordamientos, subdesbordamientos y saturaciones numéricas se señalan; nunca convierten un acceso positivo en cero físico. Los estados infinitos se representan por etiquetas explícitas para conservar JSON válido.

`triangulation_registration` cuenta etiquetas declaradas. $n\geq3$ habilita el requisito nominal triádico, **no** acredita independencia ni identificabilidad. Estas dependen del objetivo y de la ley conjunta. Procedimientos con calibración compartida no se declaran independientes por tener distintos nombres.

## Normalización y banda convencional

$I_j^\uparrow=0.90$ equivale a $\mathcal R_j=0.1$ y $\Delta_j=\ln10$. Es una identidad, no una probabilidad ni una barrera física. La banda v0.6 $[0.9\ln10,1.1\ln10]$ es convencional y distinta del intervalo histórico $[2.1,2.3]$. Su masa $\log_{10}(11/9)$ sólo tiene interpretación probabilística para significandos bajo una ley Benford justificada. No se aplica automáticamente a desfases logarítmicos o a constantes escogidas retrospectivamente. Véase [normalización](DCAE_NORMALIZATION_ES.md).

## Alcance ejecutado

`canonical_v06_audit.py` produce `results/canonical_v06_checks.json`: 17 controles finitos de equivalencia, signos, límites, referencias, cambio de protocolo, triangulación y balance de Shannon. Se conservan los controles históricos gaussianos/Fisher y FKP como ejemplos condicionados, no fundamentos de v0.6. Esta actualización no reejecuta las nuevas cifras experimentales narradas en el documento fuente, ni introduce nuevas observaciones, cotas gravitatorias o certificación de independencia.
