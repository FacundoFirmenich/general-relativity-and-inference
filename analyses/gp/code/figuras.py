"""Publication figures drawn only from the archived numerical results."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from gravitas_si.figure_language import translate_figure
R=Path(__file__).resolve().parents[1]; O=R/'figures';O.mkdir(exist_ok=True)
lang=os.environ.get('GRAVITAS_FIGLANG','es')
if lang=='en':O=O/'en';O.mkdir(exist_ok=True)
import sys
sys.path.insert(0,str(R.parents[1]))
from gravitas_si.visuals import style
style()
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
    'savefig.bbox':'tight','figure.dpi':130,'pdf.fonttype':42})
blue='#2455D6';violet='#7653C7';navy='#17365D'
def save(fig,name):
    translate_figure(fig,lang)
    fig.savefig(O/(name+'.pdf'));fig.savefig(O/(name+'.png'),dpi=180);plt.close(fig)
z=np.load(R/'results/gp_realizations_curves.npz'); j=json.loads((R/'results/gp_realizations.json').read_text())
fig,axs=plt.subplots(1,3,figsize=(11.2,3.25),constrained_layout=True)
for ax,key,title,col in [(axs[0],'joint_draw_examples','Ley conjunta del GP',blue),
                        (axs[1],'independent_draw_examples','Marginales independientes',violet)]:
    for v in z[key][:3]: ax.plot(z['x'],v,color=col,lw=.7,alpha=.75)
    ax.plot(z['x'],2*z['pointwise_sd'],'--',color=navy,lw=.8)
    ax.plot(z['x'],-2*z['pointwise_sd'],'--',color=navy,lw=.8)
    ax.set(xlim=(0,1.5),ylim=(-7,7),xlabel='Coordenada diagnóstica x',title=title)
axs[0].set_ylabel('Error de la función de coste')
for law,label,col in [('joint_GP','Conjunta',blue),('independent_marginals','Independiente',violet)]:
    rows=[r for r in j['results'] if r['law']==law]
    axs[2].semilogx([r['grid_points'] for r in rows],[r['quantile_997_std'] for r in rows],'o-',label=label,color=col)
axs[2].set(xlabel='Nodos de evaluación (sin datos nuevos)',ylabel='DE del cuantil posterior 99,7 %',title='La malla no debe inventar precisión')
axs[2].legend(frameon=False);axs[2].grid(alpha=.2)
save(fig,'ley_conjunta')

z=np.load(R/'results/gp_emulator_curves.npz');j=json.loads((R/'results/gp_emulator.json').read_text())
fig,axs=plt.subplots(1,2,figsize=(9.6,3.35),constrained_layout=True)
x=z['x'];axs[0].plot(x,z['delta_Q_exact'],color=navy,lw=2,label='Perfil exacto')
for key,label,col in [('x','GP en x',violet),('x_squared','GP en x²',blue)]:
    axs[0].plot(x,z['mean__4__'+key],label=label,color=col)
    rows=[r for r in j['results'] if r['coordinate']==key]
    axs[1].semilogy([r['n_training'] for r in rows],[100*abs(r['relative_quantile_error']) for r in rows],'o-',color=col,label=label)
nodes=np.arange(4);axs[0].scatter(nodes,nodes**4-2*j['z_hat']*nodes**2,s=22,color=navy)
axs[0].set(xlim=(0,1.7),ylim=(-2,17),xlabel='x (sin unidades; no es masa)',ylabel='ΔQ',title='Cuatro nodos: estructura del borde')
axs[1].set(xlabel='Nodos de entrenamiento',ylabel='Error absoluto relativo del cuantil (%)',title='También se conserva el caso exitoso',xticks=[4,7,13])
for ax in axs:ax.legend(frameon=False);ax.grid(alpha=.2)
save(fig,'emulador_fisico')

j=json.loads((R/'results/gp_cassini_refits.json').read_text())
rows=[r for r in j['results'] if r['model']=='constant' and r['ell_seconds']==7200]
fig,ax=plt.subplots(figsize=(7.6,2.9),constrained_layout=True)
xx=np.arange(len(rows));ax.bar(xx-.17,[100*r['naive_injection_recovery'] for r in rows],width=.34,color=violet,label='Sustraer GP bajo señal nula')
ax.bar(xx+.17,[100*r['joint_injection_recovery'] for r in rows],width=.34,color=blue,label='Ajuste conjunto GLS')
ax.set(xticks=xx,xticklabels=['3 × 10⁻¹⁴','1 × 10⁻¹³','3 × 10⁻¹³'],ylim=(0,120),xlabel='DE de frecuencia fraccional del GP supuesto',ylabel='Señal inyectada recuperada (%)')
ax.legend(frameon=False,loc='upper right');ax.grid(axis='y',alpha=.2)
for i,r in enumerate(rows):ax.text(i-.17,100*r['naive_injection_recovery']+2,f"{100*r['naive_injection_recovery']:.1f}",ha='center',fontsize=9)
save(fig,'cassini_inyeccion')
print('Saved 3 figures as PDF and PNG')
