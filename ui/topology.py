"""Bilingual visible interface to saved, reproducible finite-graph calculations."""
from pathlib import Path
import json
import numpy as np
import plotly.graph_objects as go
import streamlit as st
from research.fractalgo.core import posterior,excursion,path_edges,cell_weights

def render(root,T,chart,download_json):
    st.title(T('FractalGo · topología de la inferencia','FractalGo · inference topology'))
    st.write(T('Aplicación contenida: conectividad sobre grafos declarados. No estima dimensión fractal ni añade datos independientes.','Bounded application: connectivity on declared graphs. It estimates no fractal dimension and adds no independent data.'))
    base=root/'research/fractalgo/results'
    tabs=st.tabs([T('Explorar GP','Explore GP'),T('Cassini por pasada','Cassini by pass'),T('Contrato y reproducción','Contract & reproduction')])
    with tabs[0]:
        mode=st.selectbox(T('Experimento','Experiment'),[T('Control exacto: mismos resúmenes','Exact control: same summaries'),T('Emulador condicional','Conditional emulator')])
        threshold=st.slider(T('Umbral posterior / previa','Posterior / prior threshold'),.5,2.,1.,.1)
        if mode in ['Control exacto: mismos resúmenes','Exact control: same summaries']:
            z=st.slider('Z',-3.,3.,1.,.1);n=80;x=np.arange(n);pi=np.ones(n)/n
            curves=[('A',z*np.r_[np.ones(40),-np.ones(40)]),('B',z*np.tile([1.,-1.],40))]
            st.caption(T('Ley q=Zv, Z normal estándar; misma marginal en cada nodo. El grafo camino permanece fijo.','Law q=Zv, Z standard normal; identical marginal at each node. The path graph stays fixed.'))
        else:
            if not (base/'gp_inputs_and_draws.npz').exists():
                st.info(T('Falta la carga numérica completa. Este repositorio conserva el código; instalá las entradas autorizadas para explorar las realizaciones guardadas.','The complete numerical payload is missing. This repository retains the code; install authorized inputs to explore saved draws.'));return
            data=np.load(base/'gp_inputs_and_draws.npz');count=st.selectbox(T('Nodos','Nodes'),[51,101,201],index=2)
            draw=st.slider(T('Realización guardada','Saved draw'),0,len(data['joint_delta'])-1,0)
            prior=st.selectbox(T('Prior declarado','Declared prior'),['uniform_x','exp_minus_x'])
            ids=np.linspace(0,200,count,dtype=int);x=data['x'][ids];pi=cell_weights(x)*(1 if prior=='uniform_x' else np.exp(-x));pi/=pi.sum()
            curves=[(T('GP conjunto','Joint GP'),data['mean'][ids]+data['joint_delta'][draw,ids]),(T('Marginales independientes','Independent marginals'),data['mean'][ids]+data['independent_delta'][draw,ids])]
            st.caption(T('Control condicionado al perfil heredado; no son costes originales recuperados de INPOP.','Control conditioned on the inherited profile; these are not recovered original INPOP costs.'))
        fig=go.Figure();export=[]
        for (label,q),color in zip(curves,['#7857B5','#DB9663']):
            r=excursion(q,pi,path_edges(len(x)),threshold);p=r['posterior'];ratio=np.exp(r['log_ratio'])
            fig.add_scatter(x=x,y=ratio,name=label,mode='lines+markers',marker=dict(size=4),line=dict(color=color))
            st.metric(label+' · '+T('componentes','components'),r['components'])
            export.append(dict(law=label,threshold=threshold,node_ids=list(range(len(x))),coordinates=x.tolist(),prior_mass=pi.tolist(),cost=q.tolist(),posterior_mass=p.tolist(),components=r['node_groups'],component_masses=r['component_masses']))
        fig.add_hline(y=threshold,line_dash='dot',line_color='#667085');fig.update_layout(xaxis_title=T('Nodo / coordenada','Node / coordinate'),yaxis_title=T('Posterior / previa','Posterior / prior'));chart(fig)
        st.caption(T('Mover un umbral es exploración, no una prueba confirmatoria. Todos los nodos permanecen en la exportación.','Moving a threshold is exploration, not a confirmatory test. Every node remains in the export.'))
        download_json(export,'topology_exploration.json')
        st.image(str(base/('topology_gp_es.png' if T('es','en')=='es' else 'topology_gp_en.png')),use_container_width=True)
    with tabs[1]:
        data=json.loads((base/'cassini.json').read_text())
        pass_ids=sorted({r['pass_id'] for r in data['rows']})
        pid=st.selectbox(T('ID de pasada','Pass ID'),pass_ids,index=pass_ids.index(12))
        threshold_c=st.selectbox(T('Umbral residual','Residual threshold'),[0.,1.,2.])
        rows=[r for r in data['rows'] if r['pass_id']==pid and r['threshold']==threshold_c]
        from ui.localization import table
        st.dataframe(table(rows,T('Español','English')),hide_index=True,use_container_width=True)
        st.caption(T('Fracciones predictivas exploratorias; no son valores p calibrados. Se conserva cada hueco y el modelo instrumental heredado. Cero excedencias no equivale a probabilidad cero.','Exploratory predictive fractions, not calibrated p-values. Gaps and the inherited instrumental model are retained. Zero exceedances does not imply zero probability.'))
        st.image(str(base/('topology_cassini_es.png' if T('es','en')=='es' else 'topology_cassini_en.png')),use_container_width=True)
        download_json(data,'cassini_topology.json')
    with tabs[2]:
        st.code('python run.py topology',language='bash')
        st.markdown((root/'research/fractalgo'/T('CONTRACT_ES.md','CONTRACT_EN.md')).read_text())
        st.write(T('Validación independiente Gaia y medición COSMOS: pendientes; no se anuncian como ejecutadas.','Independent Gaia validation and COSMOS measurement: pending; neither is claimed as executed.'))
