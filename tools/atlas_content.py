"""Plain-language narrative, coupled to the graph's existing evidence."""
from html import escape

def bi(es,en): return {'es':es,'en':en}

def enhance(g):
    g['experience_version']='2.0'
    g['hero']=bi('Una prueba de gravedad es tan sólida como la información que conserva.','A gravity test is only as strong as the information it preserves.')
    g['lead']=bi('Seguimos la señal desde la respuesta física hasta la conclusión. Seis recorridos muestran qué calculamos, cómo lo comprobamos y qué cambia al interpretar el resultado.','Follow the signal from physical response to conclusion. Six paths show what we calculated, how we checked it, and what changes when the result is interpreted.')
    g['chapters']=[
      dict(route='gp',focus='gp_result',question=bi('¿Basta con reproducir las barras de error?','Is reproducing the error bars enough?'),
       answer=bi('La información decisiva está también entre los puntos.','Decisive information also lies between the points.'),
       explanation=bi('Dos emuladores pueden tener las mismas bandas marginales y producir leyes posteriores distintas. En una malla finita, la equivalencia se decide exactamente con las medias y covarianzas de contraste: Db y DKDᵀ.','Two emulators can have identical marginal bands and produce different posterior laws. On a finite grid, equivalence is determined exactly by contrast means and covariances: Db and DKDᵀ.'),
       why=bi('Esto convierte una impresión visual en una comparación verificable. Las distancias entre leyes cuantifican el efecto; FractalGo permite explorar la conectividad sobre un grafo fijo.','This turns a visual impression into a testable comparison. Distances between laws quantify the effect; FractalGo explores connectivity on a fixed graph.'),
       evidence=bi('Demostración en malla finita y controles sintéticos de covarianza.','Finite-grid proof and synthetic covariance controls.'),
       metric=bi('Db · DKDᵀ','Db · DKDᵀ'),label=bi('Criterio exacto de contraste','Exact contrast criterion')),
      dict(route='physics',focus='physical_check',question=bi('¿Qué respuesta física estamos ajustando?','Which physical response are we fitting?'),
       answer=bi('La fuerza y la calibración deben salir del mismo modelo.','Force and calibration must come from the same model.'),
       explanation=bi('Derivar el potencial Yukawa añade el factor (1 + r/λ) a la aceleración. Conservar la excentricidad cambia también la precesión. Contrastamos la fórmula de Bessel con cuadratura e integración orbital.','Differentiating the Yukawa potential adds the factor (1 + r/λ) to the acceleration. Retaining eccentricity also changes precession. We compare the Bessel formula with quadrature and orbital integration.'),
       why=bi('Una plantilla precisa permite atribuir una diferencia al ajuste correcto. Las acciones métricas, híbridas y de Palatini se distinguen por los modos físicos que admiten.','A precise template lets a difference be assigned to the correct fit. Metric, hybrid and Palatini actions are distinguished by the physical modes they support.'),
       evidence=bi('Respuesta fenomenológica y verificación perturbativa de dos cuerpos.','Phenomenological response and two-body perturbative verification.'),
       metric=bi('Bessel · cuadratura · órbita','Bessel · quadrature · orbit'),label=bi('Tres rutas de verificación','Three verification routes')),
      dict(route='cassini',focus='cassini_check',question=bi('Si los algoritmos coinciden, ¿qué queda por comprobar?','When the algorithms agree, what remains to be checked?'),
       answer=bi('La predicción local pone a prueba al modelo.','Local prediction puts the model to the test.'),
       explanation=bi('La jerarquía integra amplitud óptica, parámetros auxiliares y ruido correlacionado en 1696 épocas reconstruidas. Distintos motores numéricos coinciden. Una mezcla exploratoria mejora la predicción agregada y el análisis por pasadas localiza una discrepancia persistente.','The hierarchy integrates optical amplitude, nuisance parameters and correlated noise across 1696 reconstructed epochs. Different numerical engines agree. An exploratory mixture improves aggregate prediction, while pass-level analysis localizes a persistent discrepancy.'),
       why=bi('La pasada del 24–25 de junio de 2002 señala dónde investigar el modelo temporal y la reconstrucción instrumental. La concordancia computacional y la adecuación predictiva responden preguntas diferentes.','The June 24–25, 2002 pass identifies where to investigate the temporal model and instrumental reconstruction. Computational agreement and predictive adequacy answer different questions.'),
       evidence=bi('Reconstrucción condicionada; diagnóstico exploratorio con 1024 réplicas. Cero excedencias no equivale a probabilidad cero.','Conditional reconstruction; exploratory diagnostic with 1024 replicates. Zero exceedances is not zero probability.'),
       metric=bi('0 / 1024','0 / 1024'),label=bi('Excedencias en la pasada 12','Exceedances in pass 12')),
      dict(route='information',focus='information_check',question=bi('¿Qué información necesita esta decisión?','Which information does this decision need?'),
       answer=bi('La pérdida se mide respecto de un objetivo.','Loss is measured against a target.'),
       explanation=bi('Un resumen puede conservar una dirección útil y perder otra. Usamos Fisher, Rényi y riesgo bayesiano para distinguir esas pérdidas. K expresa la reducción de riesgo frente a decidir sin observaciones, para un prior y una pérdida declarados.','A summary can preserve one useful direction and lose another. We use Fisher, Rényi and Bayes risk to distinguish these losses. K expresses risk reduction relative to deciding without observations, for a declared prior and loss.'),
       why=bi('Podemos comparar transformaciones por lo que permiten decidir, y conservar la información local que un indicador agregado dejaría oculta.','We can compare transformations by the decisions they support, retaining local information that an aggregate indicator would hide.'),
       evidence=bi('Identidades bajo supuestos explícitos y controles generativos de mezcla.','Identities under explicit assumptions and generative mixture controls.'),
       metric=bi('K = 1 − R(E)/R(∅)','K = 1 − R(E)/R(∅)'),label=bi('Reducción relativa de riesgo','Relative risk reduction')),
      dict(route='cases',focus='cases_operator',question=bi('¿Cómo llega esta geometría a los datos?','How does this geometry reach the data?'),
       answer=bi('La covarianza decide qué contrastes sobreviven.','Covariance determines which contrasts survive.'),
       explanation=bi('En Gaia reconstruimos la covarianza por tránsito desde registros públicos. Una inflación diagonal puede conservar el modo común y alterar otros contrastes. En QUMOND verificamos por cuadratura una respuesta teórica al campo externo especificado.','For Gaia we reconstruct transit covariance from public records. Diagonal inflation can preserve the common mode while changing other contrasts. For QUMOND we verify a theoretical response to the specified external field by quadrature.'),
       why=bi('Los casos conectan la auditoría estadística con objetos físicos concretos. Cada comparación conserva su observable, su diseño y sus condiciones de contorno.','These cases connect statistical auditing to concrete physical objects. Each comparison retains its observable, design and boundary conditions.'),
       evidence=bi('Registros ESA a diseño fijo y cuadratura teórica QUMOND; controles complementarios en el suplemento.','ESA records at fixed design and theoretical QUMOND quadrature; complementary controls in the supplement.'),
       metric=bi('Gaia · QUMOND','Gaia · QUMOND'),label=bi('Dos geometrías, dos contrastes','Two geometries, two tests')),
      dict(route='dcae',focus='dcae_check',question=bi('¿Ver más equivale a inferir mejor?','Does seeing more mean inferring better?'),
       answer=bi('Cobertura y capacidad inferencial son preguntas distintas.','Coverage and inferential ability are different questions.'),
       explanation=bi('DCAE describe el cambio relativo entre la extensión física y el acceso observacional definido por un procedimiento. En las tablas COSMOS, coberturas muy parecidas producen resultados cuyo orden cambia con el observable y la tolerancia.','DCAE describes relative change between physical extent and observational access defined by a procedure. In the COSMOS tables, very similar coverage produces results whose ordering changes with the observable and tolerance.'),
       why=bi('Conservar componente, referencia y signo permite distinguir crecimiento y recuperación del acceso. La identificabilidad se examina además con el operador y la ley conjunta.','Retaining component, baseline and sign distinguishes growing disparity from recovery of access. Identifiability is additionally examined through the operator and joint law.'),
       evidence=bi('Formalización v0.6 y reanálisis de tablas. La medición cosmológica física del DCAE sigue pendiente.','Version 0.6 formalization and table reanalysis. Physical cosmological DCAE measurement remains pending.'),
       metric=bi('I↑ⱼ = [1 − ℛⱼ]₊','I↑ⱼ = [1 − ℛⱼ]₊'),label=bi('Acceso relativo por procedimiento','Procedure-specific relative access'))
    ]
    for chapter in g['chapters']:
        assert any(n['id']==chapter['focus'] for n in g['nodes'])
    # Keep the original formal tour fields usable for existing consumers.
    g['tour']=[dict(route=c['route'],focus=c['focus'],title=c['question'],body={l:c['answer'][l]+' '+c['explanation'][l] for l in ('es','en')}) for c in g['chapters']]

def static_cards(g,lang):
    bits=[]
    for c in g['chapters']:
        bits.append('<article class="chapter-card" data-route="'+c['route']+'"><span class="topic">'+escape(g['routes'][c['route']][lang])+'</span><h2>'+escape(c['question'][lang])+'</h2><p class="answer">'+escape(c['answer'][lang])+'</p><p>'+escape(c['explanation'][lang])+'</p><div class="card-bottom"><span>'+escape(c['metric'][lang])+'</span><button class="open-route" data-open="'+c['route']+'" disabled>'+('Explorar' if lang=='es' else 'Explore')+'</button></div><details class="offline-detail"><summary>'+('Por qué importa y evidencia' if lang=='es' else 'Why it matters and evidence')+'</summary><p>'+escape(c['why'][lang])+'</p><p>'+escape(c['evidence'][lang])+'</p></details></article>')
    return '\n'.join(bits)
