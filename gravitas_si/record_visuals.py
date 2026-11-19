"""Visualize printed V12 summaries; never imply a new orbital reduction."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1]
def main():
 colors=['#7653C7','#18A7D1','#9A5498','#17365D']
 plt.rcParams.update({'font.size':10,'text.color':'#17365D','axes.labelcolor':'#17365D','xtick.color':'#17365D','ytick.color':'#17365D','axes.edgecolor':'#B7B8CF','font.family':'DejaVu Sans','pdf.fonttype':42})
 for lang in ['en','es']:
  en=lang=='en';tr=lambda a,b:a if en else b
  fig,ax=plt.subplots(1,2,figsize=(11,3.7),gridspec_kw={'width_ratios':[.9,1.4]});fig.patch.set_facecolor('white')
  ranges=np.array([.1,1,3,10,100]);ratio=np.array([2.873e-7,3.39e-7,1.623e-6,1.413e-5,1.284e-3])/np.array([1.3248e-5,1.0128e-7,.4321e-6,3.698e-6,3.3474e-4])
  ax[0].axhline(1,color=colors[3],lw=.8,ls='--');ax[0].scatter(ranges,ratio,s=65,c=colors[0],zorder=3)
  ax[0].set(xscale='log',yscale='log',xlabel=tr('Yukawa range (au)','Alcance Yukawa (UA)'),ylabel=tr('Interval width: Icarus / Cacus','Anchura del intervalo: Ícaro / Cacus'),ylim=(.012,8))
  ax[0].set_title(tr('A  Five discrete ranges','A  Cinco alcances discretos'),loc='left',weight='bold',pad=12)
  xs=[.25409,.26758,.31976];ys=[2,1,0]
  for x,y,c in zip(xs,ys,colors):
   ax[1].plot([.249,x],[y,y],color=c,alpha=.45,lw=2);ax[1].scatter(x,y,color=c,s=70,zorder=3);ax[1].annotate(f'{x:.5f}',(x,y),xytext=(0,10),textcoords='offset points',ha='center',color=c,fontsize=9,weight='bold')
  ax[1].axvline(.25409,color=colors[3],ls=':',lw=1)
  ax[1].set_yticks(ys,[tr('Nominal','Nominal'),tr('Cacus Gaia variance ×1.490','Varianza Gaia Cacus ×1.490'),tr('Cacus: omit 2 Sep 2022','Cacus: retirar 2 sep 2022')])
  ax[1].set(xlim=(.247,.332),ylim=(-.5,2.55),xlabel=tr('Equal-precision crossing (au)','Cruce de igual precisión (UA)'))
  ax[1].set_title(tr('B  Sensitivity to one radar date','B  Sensibilidad a una fecha radar'),loc='left',weight='bold',pad=12)
  for a in ax:
   a.spines[['top','right']].set_visible(False);a.grid(axis='x',alpha=.12,color=colors[0]);a.tick_params(labelsize=9)
  fig.subplots_adjust(left=.085,right=.98,top=.83,bottom=.21,wspace=.98)
  fig.text(.085,.035,tr('Archived V12 summaries. Panel A uses printed endpoints; panel B shows design sensitivity, not confidence limits.','Resúmenes archivados V12. A usa extremos impresos; B muestra sensibilidad del diseño, no límites de confianza.'),fontsize=8,color=colors[3])
  folder=R/'figures'/lang;folder.mkdir(exist_ok=True)
  fig.savefig(folder/'orbital_sensitivity.pdf',bbox_inches='tight');fig.savefig(folder/'orbital_sensitivity.png',dpi=180,bbox_inches='tight');plt.close(fig)
if __name__=='__main__':main()
