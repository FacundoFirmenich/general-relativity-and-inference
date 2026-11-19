"""Bilingual publication figures from the checked published tables only."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1]
d=json.loads((R/'results/REANALYSIS.json').read_text())['cosmos_matched_selection']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'axes.labelcolor':'#24324A','xtick.color':'#24324A','ytick.color':'#24324A'})
colors=['#7254A3','#168EAE']
for lang in ['es','en']:
 es=lang=='es';fig,ax=plt.subplots(1,2,figsize=(10.3,3.45),gridspec_kw={'width_ratios':[1,1.05]})
 labels=['Recorridos · 89,8 %','Pérdida · 90,0 %'] if es else ['Walks · 89.8%','Loss · 90.0%']
 eps=list(d['tolerance_contrasts']);x=np.array(list(map(float,eps)))
 for j,key in enumerate(['walk_success','independent_loss_success']):
  y=[100*d['tolerance_contrasts'][e][key] for e in eps]
  ax[0].plot(x,y,color=colors[j],marker=['o','D'][j],lw=1.7,ms=5,label=labels[j])
 ax[0].set(xlabel='Tolerancia nominal' if es else 'Nominal tolerance',ylabel='Recuperación estimada (%)' if es else 'Estimated recovery (%)',xticks=x,ylim=(30,105))
 ax[0].legend(frameon=False,fontsize=8,loc='upper left')
 ax[0].set_title('A   El orden depende del umbral' if es else 'A   Ordering depends on threshold',loc='left',fontweight='bold',fontsize=10)
 keys=['delta_omega','f_hd','L_hd'];yy=np.arange(3)
 for j,key in enumerate(['walk','independent_loss']):
  vals=[100*d['component_relative_rmse'][k][key] for k in keys]
  pos=yy+(j-.5)*.28
  ax[1].barh(pos,vals,height=.25,color=colors[j])
  for a,v in zip(pos,vals):ax[1].text(v+.5,a,f'{v:.1f}',va='center',fontsize=8,color='#24324A')
 ax[1].set(yticks=yy,yticklabels=['Densidad','Fracción HD','Conectividad'] if es else ['Density','HD fraction','Connectivity'],xlabel='RMSE relativo (%)' if es else 'Relative RMSE (%)',xlim=(0,36))
 ax[1].invert_yaxis();ax[1].set_title('B   Ningún mecanismo domina' if es else 'B   Neither mechanism dominates',loc='left',fontweight='bold',fontsize=10)
 for a in ax:a.grid(axis='x',alpha=.12);a.set_axisbelow(True)
 fig.tight_layout(w_pad=2)
 (R/'figures').mkdir(exist_ok=True)
 for ext in ['pdf','png']:fig.savefig(R/f'figures/cosmos_{lang}.{ext}',bbox_inches='tight',dpi=180,metadata={'Creator':'GRAVITAS-SI; published COSMOS Tables 2–4'})
 plt.close(fig)
print('Two bilingual COSMOS figures generated from published table values.')
