"""Publication figures from executed MOND control records."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1]
def main():
    j=json.loads((R/'research/results/mond_diagnostics.json').read_text())
    old=json.loads((R/'research/results/mond_quadrupole.json').read_text())
    navy,cyan,violet,plum='#23385D','#258DA8','#7462B3','#AB5D99'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    for lang in ['en','es']:
        en=lang=='en';fig,ax=plt.subplots(1,2,figsize=(10,3.6),layout='constrained')
        c=j['density_coordinate_counterexample']['coordinates'];x=[np.log10(r['scale']) for r in c]
        ax[0].plot(x,[r['coordinate_shift'] for r in c],'o-',color=plum,lw=2,label=r'$\Delta\mathrm{BIC}(P)=2\ln s$')
        ax[0].plot(x,[r['log_evidence']-c[2]['log_evidence'] for r in c],'s--',color=navy,lw=2,label=r'$\Delta\ln Z=0$')
        ax[0].set_xlabel(r'Coordinate scale $\log_{10}s$' if en else r'Escala de coordenada $\log_{10}s$')
        ax[0].set_ylabel('Change under identical probability measure' if en else 'Cambio con idéntica medida probabilística')
        ax[0].set_title('Posterior height depends on units' if en else 'El pico posterior depende de las unidades',loc='left',weight='bold',fontsize=10)
        ax[0].legend(frameon=False,fontsize=9)
        a=j['fixed_physical_external_field']['scan'];b=old['interpolation_shape_scan']
        upper=j['fixed_physical_external_field']['upper_two_sided_gaussian_95_Q2_units_1e27']
        ax[1].axhspan(3.2-upper,upper,color=cyan,alpha=.16,label='Cassini 95% (Gaussian)' if en else 'Cassini 95% (gaussiano)')
        ax[1].plot([r['delta'] for r in a],[r['Q2_s_2']/1e-27 for r in a],'o-',color=violet,lw=2,label=r'Fixed $g_e$' if en else r'$g_e$ fijo')
        ax[1].plot([r['delta'] for r in b],[r['Q2_s_2']/1e-27 for r in b],'s--',color=navy,lw=1.6,label=r'Fixed $e_N$' if en else r'$e_N$ fijo')
        ax[1].set_xlabel(r'Interpolation shape $\delta$' if en else r'Forma de interpolación $\delta$')
        ax[1].set_ylabel(r'$Q_2\;[10^{-27}\,\mathrm{s}^{-2}]$')
        ax[1].set_title('Two different physical paths' if en else 'Dos recorridos físicos diferentes',loc='left',weight='bold',fontsize=10)
        ax[1].legend(frameon=False,fontsize=8)
        for aa in ax:aa.grid(alpha=.13);aa.set_axisbelow(True)
        for ext in ['pdf','png']:fig.savefig(R/f'figures/{lang}/mond_inference_controls.{ext}',dpi=180)
        plt.close(fig)
if __name__=='__main__':main()
