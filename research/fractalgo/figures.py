"""Regenerate both languages from executed results, never hardcoded ordinates."""
from pathlib import Path
import json
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];OUT=R/'research/fractalgo/results'
def main():
    gp=json.loads((OUT/'gp.json').read_text());cs=json.loads((OUT/'cassini.json').read_text())
    for lang in ['es','en']:
        en=lang=='en';fig,axes=plt.subplots(1,2,figsize=(10,3.6),layout='constrained')
        for prior,ax in zip(['uniform_x','exp_minus_x'],axes):
            for law,color,label in [('joint_GP','#7857B5','Joint GP' if en else 'GP conjunto'),('independent_marginals','#DB9663','Independent marginals' if en else 'Marginales independientes')]:
                rows=[r for r in gp['rows'] if r['prior']==prior and r['law']==law and r['threshold']==1.]
                ax.errorbar([r['nodes'] for r in rows],[r['components']['mean'] for r in rows],yerr=[2*r['components']['mcse'] for r in rows],fmt='o-',color=color,label=label,capsize=3)
            ax.set(title=('Uniform prior' if en else 'Prior uniforme') if prior=='uniform_x' else r'Prior $\propto e^{-x}$',xlabel='Nodes' if en else 'Nodos',ylabel='Mean components' if en else 'Componentes medios');ax.spines[['top','right']].set_visible(False)
        axes[0].legend(frameon=False,fontsize=8)
        fig.suptitle('Conditional GP control · threshold 1 · bars: 2 MCSE' if en else 'Control GP condicional · umbral 1 · barras: 2 MCSE',fontsize=11)
        fig.savefig(R/f'figures/topology_gp_{lang}.pdf');fig.savefig(OUT/f'topology_gp_{lang}.png',dpi=160);plt.close(fig)
        fig,axes=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
        for name,ax in zip(['reference_7200','mixture'],axes):
            rows=[r for r in cs['rows'] if r['model']==name and r['threshold']==0.]
            x=[r['pass_id'] for r in rows]
            ax.fill_between(x,[r['components_predictive_025'] for r in rows],[r['components_predictive_975'] for r in rows],color='#7857B5',alpha=.18,label='Predictive 95%' if en else 'Predictivo 95%')
            ax.plot(x,[r['components_observed_median'] for r in rows],'o',color='#DB9663',label='Observed residual median' if en else 'Mediana residual observada')
            ax.set(title=('7200 s' if name=='reference_7200' else ('Mixture' if en else 'Mezcla')),xlabel='Pass ID' if en else 'ID de pasada',ylabel='Signed excursion components' if en else 'Componentes de excursión con signo');ax.spines[['top','right']].set_visible(False)
        axes[0].legend(frameon=False,fontsize=8)
        fig.suptitle('Cassini · conditional exploratory check · not a force detection' if en else 'Cassini · control exploratorio condicional · no detección de fuerza',fontsize=11)
        fig.savefig(R/f'figures/topology_cassini_{lang}.pdf');fig.savefig(OUT/f'topology_cassini_{lang}.png',dpi=160);plt.close(fig)
    print('Four bilingual topology figures regenerated')
if __name__=='__main__':main()
