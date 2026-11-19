#!/usr/bin/env python3
"""Figures from published rows. Never simulate observations or interpolate trends."""
import csv, pathlib, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[1]
DATA=ROOT/'data'; OUT=ROOT/'figures'; OUT.mkdir(exist_ok=True)

def rows(f):
    with (DATA/f).open(newline='',encoding='utf-8') as handle: return list(csv.DictReader(handle))
walk=rows('cosmos_2026_walk_table2.csv'); loss=rows('cosmos_2026_loss_table3.csv')
fig,ax=plt.subplots(figsize=(7.2,4.4))
ax.plot([float(r['coverage']) for r in walk], [float(r['H_1.54']) for r in walk], marker='o', ms=4, linewidth=1.7,label='Recorrido espacial · H(ε=1,54)')
ax.plot([float(r['coverage']) for r in loss], [float(r['S_1.54']) for r in loss], marker='s',ms=4,linewidth=1.7,label='Pérdida independiente · S(ε=1,54)')
ax.plot([0,1],[0,1],ls=':',linewidth=1,label='Identidad cobertura = recuperación')
ax.set(xlabel='Fracción observada del campo (no escala cosmológica)',ylabel='Probabilidad de recuperación reportada',ylim=(-0.03,1.04),xlim=(0,1.02))
ax.grid(alpha=.2);ax.legend(fontsize=8,loc='upper left');fig.tight_layout()
fig.savefig(OUT/'fig01_inferibilidad_cosmos.png',dpi=190)
fig.savefig(OUT/'fig01_inferibilidad_cosmos.pdf')
plt.close(fig)

fig,ax=plt.subplots(figsize=(7.2,4.35))
c=rows('cosmos_2026_components_table4.csv');xs=np.arange(3);width=.36
ax.bar(xs-width/2,[100*float(c[0]['rel_rmse_'+k]) for k in ['delta_omega','f_hd','L_hd']],width,label='Recorrido espacial')
ax.bar(xs+width/2,[100*float(c[1]['rel_rmse_'+k]) for k in ['delta_omega','f_hd','L_hd']],width,label='Pérdida independiente')
ax.set_xticks(xs, ['Contraste de\ndensidad','Fracción de\nalta densidad','Extensión\nconectada'])
ax.set(ylabel='Error relativo (%), reportado en COSMOS',ylim=(0,37))
ax.grid(axis='y',alpha=.2);ax.legend(fontsize=8,loc='upper left');fig.tight_layout()
fig.savefig(OUT/'fig02_errores_por_observable.png',dpi=190)
fig.savefig(OUT/'fig02_errores_por_observable.pdf')
plt.close(fig)

fig,ax=plt.subplots(figsize=(7.2,4.4))
r=rows('scsn_degradation.csv');x=[float(z['sigma_h_km']) for z in r];y=[float(z['Pbnd_pct']) for z in r]
ax.plot(x,y,marker='o',linewidth=1.8)
ax.set_xscale('log');ax.set(xlabel='Incertidumbre de localización horizontal σ_h (km)',ylabel='Masa posterior en borde (%), CSV Sysmic',ylim=(-5,105))
ax.grid(alpha=.2);fig.tight_layout()
fig.savefig(OUT/'fig03_sysmic_resolucion_no_escala.png',dpi=190)
fig.savefig(OUT/'fig03_sysmic_resolucion_no_escala.pdf')
plt.close(fig)
print('3 figures in',OUT)
