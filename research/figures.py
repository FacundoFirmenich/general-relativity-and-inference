"""Publication figures from executed public-data and synthetic controls."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from physical_checks import coefficient
R=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                    'axes.spines.right':False,'axes.labelcolor':'#25334D','text.color':'#25334D',
                    'xtick.color':'#25334D','ytick.color':'#25334D','pdf.fonttype':42})
colors=['#6E52A3','#338DC2','#9C5296','#39A9C8']
rows=list(csv.DictReader((R/'research/results/gaia_covariance_transits.csv').open()))
j=json.loads((R/'research/results/inference_checks.json').read_text())
for lang in ['en','es']:
    en=lang=='en'
    fig,ax=plt.subplots(1,2,figsize=(10,3.45),layout='constrained',gridspec_kw={'width_ratios':[1.6,1]})
    for number,name,color in zip([21,1917,1566,161989],['Lutetia','Cuyo','Icarus' if en else 'Ícaro','Cacus'],colors):
        r=np.sort([float(x['ratio']) for x in rows if int(x['number_mp'])==number]); p=(np.arange(len(r))+.5)/len(r)
        ax[0].plot(r,p,'o-',ms=3,lw=1.1,color=color,label=f'{name} (n={len(r)})')
    ax[0].set_xlabel(r'$\sigma_{\rm shared}/\sigma_{\rm independent}$')
    ax[0].set_ylabel('Fraction of accepted transits' if en else 'Fracción de tránsitos aceptados')
    ax[0].legend(frameon=False,fontsize=8.5,loc='lower right')
    ax[0].set_title('ESA Gaia FPR · public observations' if en else 'ESA Gaia FPR · observaciones públicas',loc='left',fontsize=10,weight='bold')
    vals=[1,j['transit_covariance']['contrast_information_inflated']]
    ax[1].bar([0,1],vals,color=[colors[0],colors[2]],width=.58)
    ax[1].set_xticks([0,1],['Full block' if en else 'Bloque completo','Inflated diagonal' if en else 'Diagonal inflada'],fontsize=8)
    ax[1].set_ylabel('Information in a unit contrast' if en else 'Información de un contraste unitario')
    ax[1].set_ylim(0,1.17)
    for i,v in enumerate(vals):ax[1].text(i,v+.035,f'{v:.3f}',ha='center')
    ax[1].set_title('Synthetic covariance control' if en else 'Control sintético de covarianza',loc='left',fontsize=10,weight='bold')
    for a in ax:a.grid(axis='y',alpha=.12);a.set_axisbelow(True)
    fig.savefig(R/f'figures/{lang}/gaia_reconstruction.pdf');fig.savefig(R/f'figures/{lang}/gaia_reconstruction.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.7,3.1),layout='constrained')
    d=j['mcmc_pooling']; v=[d['incorrect_iid_mcse'],d['exact_mcse'],d['empirical_mcse']]
    labels=['IID assumption','Exact AR(1)','500 simulations'] if en else ['Supuesto IID','AR(1) exacto','500 simulaciones']
    ax.barh(labels,v,color=[colors[1],colors[0],colors[2]],height=.55)
    for i,x in enumerate(v):ax.text(x+.0003,i,f'{x:.6f}',va='center',fontsize=9)
    ax.set_xlim(0,.0203);ax.set_xlabel('Monte Carlo standard deviation of the pooled mean' if en else 'Desviación Monte Carlo de la media agrupada')
    ax.set_title(r'5 chains × 15000 states · $\rho_k=0.9^k$' if en else r'5 cadenas × 15000 estados · $\rho_k=0.9^k$',loc='left',fontsize=11,weight='bold')
    ax.grid(axis='x',alpha=.12);ax.set_axisbelow(True)
    fig.savefig(R/f'figures/{lang}/pooling_control.pdf');plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(10,3.5),layout='constrained')
    x=np.geomspace(.05,12,180)
    for e,c in zip([.0549,.2056,.6,.82697],colors):
        exact=coefficient(e,x);approx=np.pi*(x*(1-e*e))**2*np.exp(-x*(1-e*e))
        ax[0].semilogx(x,100*(approx/exact-1),color=c,label=f'e = {e:g}',lw=1.7)
    ax[0].axhline(0,color='#8692A6',lw=.7);ax[0].legend(frameon=False,fontsize=8)
    ax[0].set_xlabel(r'$a/\lambda$');ax[0].set_ylabel('Approximation error (%)' if en else 'Error de aproximación (%)')
    ax[0].set_title('Finite eccentricity matters' if en else 'Importa la excentricidad finita',loc='left',fontsize=11,weight='bold')
    xx=np.linspace(0,5,200)
    ax[1].plot(xx,(1+xx)*np.exp(-xx),color=colors[0],lw=2,label='Circular dynamics' if en else 'Dinámica circular')
    ax[1].plot(xx,np.exp(-xx),color=colors[2],lw=2,ls='--',label='Potential value' if en else 'Valor del potencial')
    ax[1].fill_between(xx,np.exp(-xx),(1+xx)*np.exp(-xx),color=colors[0],alpha=.12)
    ax[1].set_xlabel(r'$r_c/\lambda$');ax[1].set_ylabel(r'$-(\gamma_{\rm light}-1)/(2\alpha)$, $\alpha\to0$')
    ax[1].set_title('Calibration changes the coefficient' if en else 'La calibración cambia el coeficiente',loc='left',fontsize=11,weight='bold')
    ax[1].legend(frameon=False,fontsize=9)
    for a in ax:a.grid(alpha=.13);a.set_axisbelow(True)
    fig.savefig(R/f'figures/{lang}/physical_benchmarks.pdf');fig.savefig(R/f'figures/{lang}/physical_benchmarks.png',dpi=180);plt.close(fig)
    mond=json.loads((R/'research/results/mond_quadrupole.json').read_text())
    data=mond['interpolation_shape_scan'];xx=[r['delta'] for r in data];yy=[r['Q2_s_2']/1e-27 for r in data]
    fig,ax=plt.subplots(figsize=(9.3,3.7),layout='constrained')
    ax.axhspan(1.6-3.6,1.6+3.6,color=colors[1],alpha=.10,label=r'Cassini $\pm2\sigma$ (published)' if en else r'Cassini $\pm2\sigma$ (publicado)')
    ax.axhspan(1.6-1.8,1.6+1.8,color=colors[1],alpha=.24,label=r'Cassini $\pm1\sigma$ (published)' if en else r'Cassini $\pm1\sigma$ (publicado)')
    ax.axhline(1.6,color=colors[1],lw=1)
    ax.plot(xx,yy,'o-',color=colors[0],lw=1.8,ms=5,label=r'QUMOND · $e_N=1.643$')
    ax.scatter([1],[mond['surface_density_scenario'][-1]['Q2_s_2']/1e-27],marker='D',s=48,facecolors='none',edgecolors=colors[2],label=r'$e_N=1.659$, $\delta=1$',zorder=5)
    ax.annotate('33.868',xy=(1,yy[0]),xytext=(1.45,32),color=colors[0],fontsize=9)
    ax.set_xlim(.7,12.35);ax.set_ylim(-3,37)
    ax.set_xlabel(r'Interpolation shape $\delta$' if en else r'Forma de interpolación $\delta$')
    ax.set_ylabel(r'$Q_2\;[10^{-27}\,\mathrm{s}^{-2}]$')
    ax.set_title('A conditional QUMOND prediction · no ephemeris refit' if en else 'Predicción QUMOND condicional · sin reajustar efemérides',loc='left',fontsize=11,weight='bold')
    ax.legend(frameon=False,fontsize=9,loc='upper right');ax.grid(alpha=.12);ax.set_axisbelow(True)
    fig.savefig(R/f'figures/{lang}/mond_quadrupole.pdf');fig.savefig(R/f'figures/{lang}/mond_quadrupole.png',dpi=180);plt.close(fig)
