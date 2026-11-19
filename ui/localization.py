"""Human-readable ES/EN labels; exported scientific keys remain stable."""
import pandas as pd
LABELS={
'model':('Modelo','Model'),'pass_id':('ID de pasada','Pass ID'),'epochs':('Épocas','Epochs'),'threshold':('Umbral','Threshold'),'broken_edges':('Aristas interrumpidas','Broken edges'),
'components_observed_median':('Componentes: mediana observada','Components: observed median'),'components_predictive_025':('Componentes: predictivo 2,5 %','Components: predictive 2.5%'),'components_predictive_975':('Componentes: predictivo 97,5 %','Components: predictive 97.5%'),
'longest_observed_median':('Mayor componente: mediana observada','Longest component: observed median'),'longest_predictive_025':('Mayor componente: predictivo 2,5 %','Longest component: predictive 2.5%'),'longest_predictive_975':('Mayor componente: predictivo 97,5 %','Longest component: predictive 97.5%'),
'coverage_observed_median':('Cobertura: mediana observada','Coverage: observed median'),'coverage_predictive_median':('Cobertura: mediana predictiva','Coverage: predictive median'),
'components_lower_tail_count':('Recuento: cola inferior de componentes','Count: component lower tail'),'longest_upper_tail_count':('Recuento: cola superior de longitud','Count: longest upper tail'),'draws':('Realizaciones','Draws'),'components_lower_tail':('Fracción: cola inferior de componentes','Fraction: component lower tail'),'longest_upper_tail':('Fracción: cola superior de longitud','Fraction: longest upper tail'),
'prior_variance':('Varianza previa','Prior variance'),'Fisher_mixture':('Fisher de la mezcla','Mixture Fisher'),'prior_Fisher':('Fisher previa','Prior Fisher'),'Bayes_risk':('Riesgo de Bayes','Bayes risk'),'van_Trees_lower_bound':('Cota inferior de van Trees','Van Trees lower bound'),'oracle_Bayes_risk':('Riesgo de Bayes oracular','Oracle Bayes risk'),'retained_K':('Retención K','Retained K'),'oracle_K':('K oracular','Oracle K'),'risk_quad_error':('Error de cuadratura del riesgo','Risk quadrature error'),'variance_identity_error':('Error de identidad de varianza','Variance identity error')}
def table(rows,language):
 d=pd.DataFrame(rows);idx=0 if language=='Español' else 1
 if 'model' in d:d['model']=d['model'].replace({'reference_7200':('Referencia 7200 s','7200 s reference')[idx],'mixture':('Mezcla','Mixture')[idx]})
 return d.rename(columns={k:v[idx] for k,v in LABELS.items()})
FINDINGS_ES={
'G01':('Precesión','La expresión original tiene unidades de s⁻¹, pero se presenta en rad/órbita; el promedio en anomalía verdadera también omite el peso temporal.'),
'G02':('Densidad de Poisson','El signo exterior está invertido; la densidad efectiva de Poisson no constituye un dictamen sobre condiciones de energía.'),
'G03':('Hipótesis óptica','Los extremos de integración fijos introducen truncamiento; se comprueba la sustitución mediante cuadratura con extremos infinitos.'),
'G04':('Fisher multisondeo','El original reutiliza una plantilla provisional de Mercurio en todos los sondeos; cada uno requiere su observable y covarianza conjunta.'),
'G05':('Gráfico de restricciones','La curva plana Cassini alfa = 2,3×10⁻⁵ y otras curvas son esquemáticas; no son restricciones derivadas de una verosimilitud.'),
'G06':('Gráfico de densidad','El denominador es una densidad Yukawa de referencia, aunque la etiqueta afirma rho_GR = 0.'),
'G07':('ppE, estabilidad y energía','El potencial estático no especifica radiación, acción perturbativa ni presiones; el signo de alfa no permite un dictamen universal sobre taquiones.'),
'G08':('Factores de Bayes','compute_bayes_factor devuelve log B, pero jeffreys_scale espera B; las escalas no son intercambiables.')}
