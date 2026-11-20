"""Gravitas · visual companion to the Yukawa paper. Start via ABRIR_GRAVITAS."""
from pathlib import Path
import sys, json, io, zipfile, time
import numpy as np
import plotly.graph_objects as go
import streamlit as st
from scipy.special import softmax

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ui.physics import YukawaParams, radial_acceleration, precession_per_orbit, AU, MU_SUN
from research.revision_v17.audit import gaussian_distance
from research.revision_v18.information_audit import renyi_gaussian
from ui.jobs import execute

st.set_page_config(page_title='Gravitas · General Relativity', page_icon='◈', layout='wide')
st.markdown('''<style>
.stApp{background:#f8f9fc;color:#192037}
[data-testid="stSidebar"]{background:#f0edf6;border-right:1px solid #e1deea}
h1,h2,h3{font-family:Georgia,serif!important;letter-spacing:-.025em}
h1{font-size:2.4rem!important} h2{font-size:1.6rem!important}
.hero{padding:30px 34px;background:linear-gradient(118deg,#182940,#473766 65%,#715170);border-radius:16px;color:white;margin:8px 0 26px}
.hero h1{color:white}.hero p{color:#e2e5f2;max-width:810px;font-size:1.06rem}
.eyebrow{font:600 .73rem sans-serif;letter-spacing:.18em;text-transform:uppercase;color:#d7d0ee}
.small-note{color:#596176;font-size:.88rem;line-height:1.6}
[data-testid="stMetric"]{background:white;border:1px solid #e2e5ed;border-radius:12px;padding:14px}
.stButton>button{border-radius:9px;border-color:#756397}
</style>''',unsafe_allow_html=True)

lang=st.sidebar.selectbox('Idioma / Language',['Español','English'])
def T(es,en):return es if lang=='Español' else en
st.sidebar.markdown('## GRAVITAS\nFacundo Firmenich')
st.sidebar.caption(T('Relatividad general · laboratorio visual','General relativity · visual laboratory'))
pages=[T(a,b) for a,b in [('Inicio','Overview'),('Dinámica Yukawa','Yukawa dynamics'),('Cassini','Cassini'),('Emuladores GP','GP emulators'),('Información','Information'),('Reproducir','Reproduce'),('Documentos y procedencia','Documents & provenance'),('Grafo epistémico','Epistemic graph'),('Fermi-LAT · auditoría','Fermi-LAT · audit'),('FractalGo · topología','FractalGo · topology')]]
page=st.sidebar.radio(T('Explorar','Explore'),pages,label_visibility='collapsed')
st.sidebar.markdown('[ORCID ↗](https://orcid.org/0009-0002-6578-3811)')
st.sidebar.caption(T('Módulo del paper. Se conserva la plataforma original.','Paper module. The original platform is preserved.'))

COLORS=['#3c4777','#875d96','#298ca0','#a1769e','#667eaa']
def chart(fig,height=410):
    fig.update_layout(template='plotly_white',height=height,font=dict(family='Arial',color='#242e43'),
                      colorway=COLORS,margin=dict(l=35,r=25,t=35,b=40),paper_bgcolor='rgba(0,0,0,0)',
                      legend=dict(orientation='h',y=1.08))
    st.plotly_chart(fig,use_container_width=True,config={'displaylogo':False,'toImageButtonOptions':{'format':'svg'}})

def read(relative):return json.loads((ROOT/relative).read_text())
def download_json(data,name,label=None):
    st.download_button(label or T('Descargar resultado JSON','Download JSON result'),json.dumps(data,indent=2,ensure_ascii=False),file_name=name,mime='application/json')

if page==pages[0]:
    st.markdown('<div class="hero"><div class="eyebrow">GRAVITAS / '+T('LABORATORIO DE INVESTIGACIÓN','RESEARCH COMPANION')+'</div><h1>'+T('De la gravedad a la evidencia.','From gravity to evidence.')+'</h1><p>'+T('Explorá las ecuaciones, examiná los datos y reproducí los resultados del paper. Cada vista conserva la diferencia entre un cálculo físico, un diagnóstico estadístico y una observación.','Explore equations, inspect data, and reproduce the paper’s results. Every view distinguishes physical calculations, statistical diagnostics and observations.')+'</p></div>',unsafe_allow_html=True)
    a,b,c=st.columns(3)
    a.metric(T('Épocas reconstruidas','Reconstructed epochs'),'1,696')
    b.metric(T('Pasos temporales Cassini','Cassini observing passes'),'20')
    c.metric(T('Componentes de ruido','Noise components'),'4')
    st.caption(T('Reconstrucción condicional DSS-25; no es una nueva reducción de la telemetría original.','Conditional DSS-25 reconstruction; not a new reduction of the original telemetry.'))
    left,right=st.columns([1.2,1])
    with left:
        st.subheader(T('Una ruta de lectura','A reading route'))
        st.markdown(T('**1 · Dinámica.** Cambiá el alcance y la amplitud de la fuerza efectiva.\n\n**2 · Inferencia.** Inspeccioná la jerarquía Cassini y sus discrepancias por paso.\n\n**3 · Información.** Compará covarianzas, compresión y tratamiento de colas.\n\n**4 · Reproducción.** Ejecutá cálculos, tablas o figuras y descargá el registro.','**1 · Dynamics.** Change the effective force range and amplitude.\n\n**2 · Inference.** Inspect the Cassini hierarchy and pass-level discrepancies.\n\n**3 · Information.** Compare covariance, compression and tail handling.\n\n**4 · Reproduction.** Run calculations, tables or figures and download the record.'))
    with right:
        e=np.linspace(0,.9,28);x=np.logspace(-1,1.2,36)
        z=np.array([[precession_per_orbit(AU,ee,YukawaParams(1,AU/xx)) for xx in x] for ee in e])
        fig=go.Figure(go.Surface(x=np.log10(x),y=e,z=z,colorscale=[[0,'#e9edf8'],[.45,'#6c8ba7'],[1,'#745284']],showscale=False))
        fig.update_layout(scene=dict(xaxis_title='log₁₀(a/λ)',yaxis_title='e',zaxis_title='Δϖ/α',camera=dict(eye=dict(x=1.6,y=-1.6,z=1.1))))
        chart(fig,340)
        st.caption(T('Coeficiente perturbativo de precesión; no es una superficie observada.','Perturbative precession coefficient; not an observed surface.'))

elif page==pages[1]:
    st.title(T('Dinámica Yukawa','Yukawa dynamics'))
    st.write(T('Plantilla fenomenológica aditiva: el límite α = 0 recupera Newton.','Additive phenomenological template: α = 0 recovers Newton.'))
    st.latex(r'\Phi(r)=-\frac{\mu}{r}\left(1+\alpha e^{-r/\lambda}\right),\quad a_r=-\frac{\mu}{r^2}\left[1+\alpha(1+r/\lambda)e^{-r/\lambda}\right]')
    a,b,c=st.columns(3)
    alpha=a.number_input('α',min_value=-.01,max_value=.01,value=.0001,step=.0001,format='%.5f')
    loglam=b.slider('log₁₀(λ / AU)',-3.,3.,-1.,.05)
    ecc=c.slider(T('Excentricidad e','Eccentricity e'),0.,.95,.2056,.0001,format='%.4f')
    semimajor=st.number_input('a / AU',min_value=.001,max_value=100.,value=.387,step=.01)
    pars=YukawaParams(alpha,AU*10**loglam)
    radius=np.logspace(-3,3,400)*AU
    correction=radial_acceleration(radius,MU_SUN,pars)/(-MU_SUN/radius**2)-1
    fig=go.Figure(go.Scatter(x=radius/AU,y=correction,mode='lines',line=dict(color=COLORS[1],width=3),name='Δa / aN'))
    fig.update_xaxes(type='log',title='r / AU');fig.update_yaxes(title='Δa / aN')
    chart(fig)
    prec=precession_per_orbit(semimajor*AU,ecc,pars)
    st.metric(T('Precesión de primer orden · rad/órbita','First-order precession · rad/orbit'),f'{prec:.9g}')
    st.caption(T('Aproximación perturbativa de dos cuerpos, a excentricidad finita. Para e = 0 se muestra el límite epicíclico. Esta fuerza no determina por sí sola deflexión luminosa, ppE ni estabilidad de una teoría covariante.','Perturbative two-body approximation at finite eccentricity. At e = 0 this is the epicyclic limit. This force alone does not determine light deflection, ppE or covariant-theory stability.'))
    download_json(dict(alpha=alpha,lambda_m=pars.lambda_m,a_m=semimajor*AU,e=ecc,precession_rad_per_orbit=prec,scope='First order in alpha; illustrative inputs; not fitted bounds'),'yukawa_calculation.json')

elif page==pages[2]:
    st.title('Cassini · '+T('jerarquía y diagnóstico','hierarchy & diagnostics'))
    j=read('research/revision_v17/results/cassini.json')
    st.write(T('Resultados ejecutados sobre 1.696 épocas reconstruidas. La amplitud óptica A no se identifica con α, γ ni una masa de gravitón.','Executed results for 1,696 reconstructed epochs. Optical amplitude A is not identified with α, γ or a graviton mass.'))
    st.latex(r'y=X\beta+a\,t+\epsilon,\quad \epsilon\sim\mathcal N(0,C_{s,\tau,\ell}),\quad p(\ell)=1/4')
    a,b=st.columns([1,1.6]); probs=j['mixture_ppc']['kernel_posterior_probabilities']
    with a:
        st.subheader(T('Mezcla de escalas temporales','Temporal scale mixture'))
        chart(go.Figure(go.Bar(x=[str(int(x))+' s' for x in j['kernel_lengths_seconds']],y=probs,marker_color=COLORS[:4],name='p(ℓ | y)')),300)
        st.caption(T('Pesos posteriores condicionales. La mezcla se exploró después del diagnóstico inicial.','Conditional posterior weights. The mixture was explored after the initial diagnostic.'))
    with b:
        model=st.selectbox(T('Control predictivo','Predictive check'),[T('Mezcla','Mixture'),T('Referencia 7200 s','7200 s reference')])
        q=j['mixture_ppc'] if model in ['Mezcla','Mixture'] else j['reference_ppc']
        x=np.array([r['pass_id'] for r in j['pass_predictions']]);mid=(np.array(q['pass_replicated_975'])+q['pass_replicated_025'])/2
        fig=go.Figure(go.Scatter(x=x,y=mid,mode='markers',error_y=dict(type='data',array=np.array(q['pass_replicated_975'])-mid),name=T('Intervalo predictivo 95%','95% predictive interval'),marker_color=COLORS[2]))
        fig.add_scatter(x=x,y=q['pass_observed_median'],mode='markers',name=T('Mediana observada','Observed median'),marker=dict(color=COLORS[1],size=9,symbol='diamond'))
        fig.update_xaxes(title=T('ID de pasada','Pass ID'),dtick=1);fig.update_yaxes(title=T('Mediana del residuo normalizado','Median standardized residual'))
        chart(fig,340)
    st.subheader(T('Inspección por paso','Pass-level inspection'))
    pass_ids=[r['pass_id'] for r in j['pass_predictions']]
    selected=st.selectbox(T('Pasada Cassini','Cassini pass'),pass_ids,index=pass_ids.index(12)); row=next(r for r in j['pass_predictions'] if r['pass_id']==selected)
    st.write(f"{row['start_utc']} → {row['end_utc']} · n = {row['n']}")
    l,m,n=st.columns(3);l.metric(T('Fracción predictiva · mezcla','Predictive fraction · mixture'),f"{row['mixture_ppc_tail']:.6f}");m.metric(T('Cambio de log score','Log-score change'),f"{row['difference']:.4f}");n.metric(T('Simulaciones por modelo','Replicates per model'),'1,024')
    st.caption(T('Una fracción simulada cero significa cero excedencias en 1.024 réplicas, no probabilidad física cero. Borrar un paso conserva el diseño reconstruido y los pesos empíricos: es validación interna condicionada.','A zero simulated fraction means zero exceedances in 1,024 replicates, not zero physical probability. Pass deletion retains the reconstructed design and empirical weights: this is conditional internal validation.'))
    with st.expander(T('Supuestos y límites completos','Full assumptions and limitations')):
        st.write(T('Condicionado a residuos híbridos DSS-25 reconstruidos, pesos fijos, alcance óptico fijo y base auxiliar del diseño completo. La exclusión de pasadas no es validación externa independiente.',j['scope']))
        limits_es=['Sin comparación entre estaciones: el diseño corresponde a DSS-25.','Pesos y reconstrucción previos al diagnóstico; sin validación cruzada de datos crudos.','Alternativas temporales elegidas tras observar el fallo original; no confirmación prerregistrada.','El prior del diseño dependiente de datos permanece fijo entre particiones.','Sin posterior físico de alfa/lambda ni jerarquía de todos los parámetros de efemérides.']
        for es,en in zip(limits_es,j['limits']):st.write('• '+T(es,en))
    download_json(j,'cassini_hierarchy.json')

elif page==pages[3]:
    st.title(T('Auditoría de emuladores GP','GP emulator audit'))
    st.write(T('Mismas bandas marginales, distinta ley conjunta. Explorador gaussiano finito de tres nodos; no reproduce una matriz INPOP publicada.','Same marginal bands, different joint laws. A finite three-node Gaussian explorer; not a reproduction of a published INPOP matrix.'))
    rho=st.slider(T('Correlación compartida ρ','Shared correlation ρ'),0.,.99,.9,.01)
    D=np.array([[1,0,-1],[0,1,-1]],float);K=(1-rho)*np.eye(3)+rho*np.ones((3,3));S0=D@K@D.T;S1=D@D.T
    distances=gaussian_distance(np.zeros(2),S0,np.zeros(2),S1)
    a,b,c=st.columns(3);a.metric('KL · nats',f"{distances['KL']:.6f}");b.metric('H²',f"{distances['H2']:.6f}");c.metric(T('Razón de varianzas de contraste','Contrast variance ratio'),f'{1-rho:.3f}')
    t=np.linspace(0,2*np.pi,200);unit=np.array([np.cos(t),np.sin(t)])
    fig=go.Figure()
    for S,label,color in [(S0,T('Covarianza coherente','Coherent covariance'),COLORS[1]),(S1,T('Independencia postulada','Assumed independence'),COLORS[2])]:
        xy=np.linalg.cholesky(S)@unit
        fig.add_scatter(x=xy[0],y=xy[1],mode='lines',name=label,line=dict(color=color,width=3))
    fig.update_xaxes(title='q₁ − q₃');fig.update_yaxes(title='q₂ − q₃',scaleanchor='x',scaleratio=1);chart(fig)
    renyi=[renyi_gaussian(np.zeros(2),S0,np.zeros(2),S1,x) for x in [.5,1.,2.]]
    st.latex(r'p_i=\frac{w_i e^{-q_i}}{\sum_j w_j e^{-q_j}},\qquad \mathcal L(p)\longleftrightarrow(Db,DKD^T)')
    st.caption(T('Contornos de Mahalanobis unidad, no regiones de probabilidad 68%. La equivalencia presupone pesos positivos fijos; tolerancia numérica de rango 10⁻¹⁰.','Unit Mahalanobis contours, not 68% probability regions. Equivalence assumes fixed positive weights; numerical rank tolerance 10⁻¹⁰.'))
    st.subheader(T('Presupuesto de error del posterior','Posterior error budget'))
    omega=st.slider(T('Oscilación del error de coste efectivo ω','Effective-cost error oscillation ω'),0.,4.,.4,.01)
    x,y=st.columns(2);x.metric(T('Cota de variación total','Total-variation bound'),f'{np.tanh(omega/4):.6f}');y.metric('KL · '+T('cota en nats','bound in nats'),f'{omega**2/8:.6f}')
    st.latex(r'\delta_i=\frac N2(\widetilde q_i-q_i),\quad\omega=\max_i\delta_i-\min_i\delta_i')
    st.caption(T('Cota analítica entre dos distribuciones normalizadas sobre la malla. No es distancia entre las leyes de vectores posteriores aleatorios ni un error INPOP medido.','Analytic bound between two normalized distributions on the grid. Not a distance between laws of random posterior vectors or a measured INPOP error.'))
    download_json(read('research/revision_v19/results/posterior_stability.json'),'posterior_stability.json')
    download_json(dict(rho=rho,contrast_covariance_reference=S0.tolist(),contrast_covariance_alternative=S1.tolist(),distances=distances,renyi=renyi),'gp_contrast_audit.json')

elif page==pages[4]:
    with st.expander(T('DCAE v0.6 · acceso por procedimiento', 'DCAE v0.6 · procedure-specific access')):
        st.latex(r'\mathcal R_j=\frac{L_{O,j}/L_F}{L_{O,j,0}/L_{F,0}},\quad \Delta_j=-\ln\mathcal R_j,\quad I_j^\uparrow=[1-\mathcal R_j]_+')
        st.write(T('La razón anterior R es 1/ℛⱼ. Se conservan las componentes, el signo, el protocolo y las dependencias. Tres procedimientos inician una tríada; no prueban independencia ni identificabilidad.', 'The earlier ratio R equals 1/ℛⱼ. Components, sign, protocol and dependencies are retained. Three procedures initiate a triad; they do not prove independence or identifiability.'))
        st.caption(T('Controles matemáticos ejecutados; la actualización no incorpora nuevos datos cosmológicos ni nuevas cotas de gravedad.', 'Executed mathematical controls; this update adds no new cosmological data or gravity bounds.'))
        download_json(read('research/dcae/results/canonical_v06_checks.json'),'dcae_v06_checks.json')
    st.info(T('DCAE: I = 0.90 equivale por definición a R = 10 y Δ = ln(10). No es una probabilidad ni una barrera física. Conservá R, el desfase firmado y sus incertidumbres.', 'DCAE: I = 0.90 means R = 10 and Δ = ln(10) by definition. It is neither a probability nor a physical barrier. Retain R, the signed gap and their uncertainties.'))

    st.title(T('Información: qué se conserva','Information: what is preserved'))
    j=read('research/revision_v18/results/information_audit.json')
    tabs=st.tabs([T('Señal y RMS','Signal & RMS'),T('Colas y censura','Tails & censoring'),T('Jerarquía y Fisher','Hierarchy & Fisher'),T('Umbrales','Thresholds')])
    with tabs[0]:
        rows=j['rms']['rows'];fig=go.Figure()
        for key,label,color in [('matched_power',T('Contraste firmado bilateral','Two-sided signed test'),COLORS[1]),('RMS_power','RMS',COLORS[2])]:
            fig.add_scatter(x=[r['n'] for r in rows],y=[r[key] for r in rows],mode='lines+markers',name=label,line=dict(color=color))
        fig.update_xaxes(type='log',title='n');fig.update_yaxes(title=T('Potencia','Power'),range=[0,1]);chart(fig)
        st.caption(T('Control exacto: Yi ~ N(θ,1), θ = 3/√n, nivel 0,05. La RMS pierde el signo y tiene información de Fisher local nula en θ = 0. No es una afirmación sobre toda función de coste cuadrática.','Exact control: Yi ~ N(θ,1), θ = 3/√n, level 0.05. RMS loses the sign and has zero local Fisher information at θ = 0. This is not a claim about every quadratic cost function.'))
    with tabs[1]:
        rows=j['tails']['rows'];fig=go.Figure()
        for key,label,color in [('retained_inside_Fisher_per_attempt',T('Sólo interior + conteos','Interior only + counts'),COLORS[1]),('two_sided_censored_Fisher_per_attempt',T('Interior + etiqueta de cola','Interior + tail label'),COLORS[2])]:
            fig.add_scatter(x=[r['central_probability'] for r in rows],y=[r[key] for r in rows],mode='lines+markers',name=label,line=dict(color=color))
        fig.update_xaxes(title=T('Masa central conservada','Retained central mass'));fig.update_yaxes(title=T('Fisher por intento','Fisher per attempt'),range=[0,1.02]);chart(fig)
        st.write(T('Conservar el 95% central retiene 0,72090 de la información de localización por intento en este modelo. Conservar además qué cola se censuró retiene 0,99417.','Retaining the central 95% retains 0.72090 of location information per attempt in this model. Also retaining which tail was censored retains 0.99417.'))
        st.caption(T('N(θ,1) evaluada en θ = 0; cortes fijos en coordenadas medidas. La información por observación retenida tiene otro denominador. Transportar masa mediante un mapa invertible conocido conserva el experimento; conservar sólo la masa total no lo garantiza.','N(θ,1) evaluated at θ = 0; fixed cuts in measured coordinates. Information per retained observation has a different denominator. Transport through a known invertible map preserves the experiment; total-mass conservation alone does not guarantee this.'))
    with tabs[2]:
        st.latex(r'\mathcal K_{\pi,L}=1-\frac{R^*_{\pi,L}(Y)}{R^*_{\pi,L}(\varnothing)}')
        st.write(T('Retención de capacidad decisoria relativa al prior y a la pérdida declarados. No es un índice universal de verdad.','Decision-capacity retention relative to the declared prior and loss. Not a universal truth index.'))
        from ui.localization import table
        st.dataframe(table(j['mixture']['scenarios'],lang),use_container_width=True)
        st.caption(T('Mezcla sintética de ruido N(0,1) / N(0,9) con etiqueta latente. La información con etiqueta conocida es un comparador oracular, no la información disponible. SBC de esta jerarquía de juguete: 4.000 réplicas; no valida toda la jerarquía Cassini.','Synthetic N(0,1) / N(0,9) noise mixture with latent label. Known-label information is an oracle comparator, not available information. SBC of this toy hierarchy: 4,000 replicates; it does not validate the entire Cassini hierarchy.'))
    with tabs[3]:
        fraction=st.slider(T('Fracción de saturación','Saturation fraction'),.50,.99,.90,.01)
        st.metric('t / τ · 1 − exp(−t/τ)',f'{-np.log1p(-fraction):.8f}')
        st.caption(T('Identidad condicional del modelo exponencial. A 95% resulta ln(20), no ln(10). La coincidencia de mantisas no establece que Cassini usara un criterio de parada del 90%.','Conditional identity of the exponential model. At 95% the result is ln(20), not ln(10). A mantissa coincidence does not establish that Cassini used a 90% stopping rule.'))
    download_json(j,'information_audit.json')

elif page==pages[5]:
    st.title(T('Reproducir en un clic','One-click reproduction'))
    st.write(T('Los botones ejecutan el código del paquete con sus entradas incluidas. Cada ejecución conserva su log y un registro propio.','Buttons execute the packaged code with included inputs. Every execution preserves its log and its own record.'))
    options={T('FractalGo: topología','FractalGo: topology'):'topology',T('Fermi-LAT: catálogo y eventos','Fermi-LAT: catalogue and events'):'fermi',T('Estabilidad del posterior','Posterior stability'):'stability',T('Todo: cálculos, tablas y figuras','Everything: calculations, tables & figures'):'all',T('Sólo cálculos','Calculations only'):'calculations',T('Sólo figuras','Figures only'):'figures',T('Sólo tablas','Tables only'):'tables',T('Auditoría GraviTools','GraviTools audit'):'platform',T('Controles de información','Information controls'):'information',T('Jerarquía Cassini y revisión','Cassini hierarchy & revision'):'revision',T('Compilar paper y suplemento (requiere TeX)','Compile paper & supplement (requires TeX)'):'documents'}
    if not (ROOT/'publication/en/paper.md').exists():
        options={k:v for k,v in options.items() if v!='documents'}
    choice=st.selectbox(T('Qué ejecutar','What to run'),list(options))
    st.caption(T('Los cálculos completos pueden tardar varios minutos. PDF requiere además Pandoc y XeLaTeX.','Full calculations may take several minutes. PDFs additionally require Pandoc and XeLaTeX.'))
    if st.button(T('Ejecutar ahora','Run now'),type='primary'):
        box=st.empty();lines=[]
        def report(line):
            lines.append(line);box.code(''.join(lines[-18:]),language='text')
        with st.spinner(T('Calculando…','Computing…')):
            try:st.session_state['last_job']=execute(options[choice],report)
            except Exception as exc:st.error(str(exc))
    if 'last_job' in st.session_state:
        job=st.session_state['last_job'];st.write(T('Estado: ','Status: ')+({'completed':'Completado','failed':'Falló'}.get(job['status'],job['status']) if lang=='Español' else job['status']))
        download_json(job,'execution_record.json')
        log=ROOT/job['log'];st.download_button(T('Descargar log completo','Download full log'),log.read_bytes(),file_name=log.name,mime='text/plain')

elif page==pages[6]:
    st.title(T('Documentos y procedencia','Documents & provenance'))
    st.write(T('Paper, supplementary y fuentes originales se conservan separados de los resultados de exploración.','Paper, supplementary and original sources are kept separate from exploration results.'))
    edition=read('ui/edition.json')
    st.caption(T('Edición documental: ','Document edition: ')+edition['pdf_version'])
    if not (ROOT/'PAPER_ES.pdf').exists():st.info(T('El manuscrito y su suplemento se entregan por separado; no están incluidos en este repositorio de software.','The manuscript and supplement are delivered separately; they are not included in this software repository.'))
    a,b=st.columns(2)
    for i,(name,label) in enumerate([('PAPER_ES.pdf','Paper · Español'),('PAPER_EN.pdf','Paper · English'),('SUPPLEMENTARY_ES.pdf','Supplementary · Español'),('SUPPLEMENTARY_EN.pdf','Supplementary · English')]):
        p=ROOT/name
        if p.exists():[a,b][i%2].download_button(label,p.read_bytes(),file_name=name,mime='application/pdf')
    st.subheader(T('Constantes y Sistema Internacional','Constants and the International System'))
    from ui.constants import REGISTRY
    for row in REGISTRY:
        st.markdown(f"**{row['symbol']}** = {format(row['value_si'], '.5e') if row['symbol']=='G' else format(row['value_si'], '.12g')} {row['unit']} · "+row['status_es' if lang=='Español' else 'status_en']+f" · u = {row['standard_uncertainty_si']:.12g}")
    download_json(REGISTRY,'si_constants.json')
    st.subheader(T('GraviTools: recuperación y verificación','GraviTools: recovery & verification'))
    audit=read('platform/audit_runs/verification.json')
    for f in audit['findings']:
        from ui.localization import FINDINGS_ES
        title,body=FINDINGS_ES[f['id']]
        with st.expander(f['id']+' · '+T(title,f['component'])):
            st.write(T(body,f['finding']))
            st.caption(T('Registro original con claves técnicas estables:','Original record with stable technical keys:'));st.json(f)
    download_json(audit,'gravitools_verification.json')
    st.markdown('[GraviTools original · Colab](https://colab.research.google.com/drive/1itRdXq6j8HqIRTWDGI0mgAeJjOCDtPl1)')
    st.caption(T('Nexus y Gravitas se reconocen como desarrollos previos del autor. La documentación Nexus está recuperada; su motor Python no se ha importado ni certificado en este módulo.','Nexus and Gravitas are recognized as the author’s prior developments. Nexus documentation has been recovered; its Python engine has not been imported or certified in this module.'))


if page==pages[7]:
    st.title(T('Grafo epistémico','Epistemic graph'))
    st.caption(T('Tour guiado, reproducción con pausa, fuentes descargables y límites de cada afirmación.','Guided tour, pausable playback, downloadable sources and claim boundaries.'))
    from streamlit.components.v1 import html as embed_html
    graph_path=ROOT/'ui/epistemic'/('epistemic_graph_es.html' if lang=='Español' else 'epistemic_graph_en.html')
    st.download_button(T('Abrir fuera de Gravitas: descargar HTML autónomo','Open outside Gravitas: download standalone HTML'),graph_path.read_bytes(),file_name=graph_path.name,mime='text/html')
    embed_html(graph_path.read_text(),height=1050,scrolling=True)


if page==pages[8]:
    import pandas as pd
    st.title(T('Fermi-LAT: del registro al resultado','Fermi-LAT: from record to result'))
    st.info(T('Control complementario de procedencia e inferencia. No es una prueba Yukawa ni una nueva confirmación o refutación de RG.','Complementary provenance/inference control. Not a Yukawa test or new confirmation or refutation of GR.'))
    base=ROOT/'platform/fermi_review'
    failure=read('platform/fermi_review/results/hunter_exact_replay.json')
    st.subheader(T('El fallo reproducido','The reproduced failure'))
    st.write(T('Un delimitador desplazaba las columnas: REDSHIFT recibía distancia de luminosidad y TL0 recibía TL100. La salida −0.0000 ± 0.0000 no era una medición nula.','A delimiter shifted columns: REDSHIFT received luminosity distance and TL0 received TL100. The −0.0000 ± 0.0000 output was not a null measurement.'))
    download_json(failure,'hunter_exact_replay.json')
    j=read('platform/fermi_review/photon_followup/results/execution.json')
    end=st.selectbox(T('Ventana después del trigger (s)','Window after trigger (s)'),[100,10,1000])
    rows=[r for r in j['contrasts'] if r['window_end_s']==end]
    frame=pd.DataFrame(rows).drop(columns=['window_end_s','status'])
    frame.columns=[T('GRB','GRB'),T('N baja','N low'),T('N alta','N high'),T('Diferencia de medianas (s)','Median difference (s)'),T('Bootstrap 2.5% (s)','Bootstrap 2.5% (s)'),T('Bootstrap 97.5% (s)','Bootstrap 97.5% (s)')]
    st.dataframe(frame,hide_index=True,use_container_width=True)
    st.caption(T('4096 réplicas; intervalo condicional a los eventos seleccionados. Sin IRF/exposición efectiva, fondo ni emisión intrínseca modelados. Los valores vacíos indican insuficiencia, no cero.','4096 replicates; interval conditional on selected events. No modelled IRF/effective exposure, background or intrinsic emission. Empty values indicate insufficient events, not zero.'))
    preview=base/'photon_followup/results'/('photon_timing_ES.png' if lang=='Español' else 'photon_timing_EN.png')
    if preview.exists():st.image(str(preview),use_container_width=True)
    else:st.info(T('La figura se regenera con el paquete numérico completo: python run.py fermi.','Regenerate the figure with the complete numerical package: python run.py fermi.'))
    if st.button(T('Reejecutar selección y controles','Rerun selection and checks'),key='fermi_run'):
        with st.spinner(T('Ejecutando…','Running…')):result=execute('fermi')
        st.write(result);download_json(result,'fermi_execution.json')
    download_json(j,'photon_timing_execution.json')

if page in ['FractalGo · topología','FractalGo · topology']:
    st.link_button(T('Abrir repositorio FractalGo', 'Open FractalGo repository'), 'https://github.com/FacundoFirmenich/fractal-go')
    from ui.topology import render
    render(ROOT,T,chart,download_json)
