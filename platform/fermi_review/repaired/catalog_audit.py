"""Audited Fermi GRB catalogue exploration, not a Lorentz-invariance fit."""
from pathlib import Path
import gzip,re,json,hashlib
import numpy as np
import pandas as pd
from scipy.integrate import quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,LogNorm
ROOT=Path(__file__).resolve().parents[1]
SOURCE='https://heasarc.gsfc.nasa.gov/FTP/heasarc/dbase/dump/heasarc_fermilgrb.tdat.gz'
def parse_catalog(raw):
 text=gzip.decompress(raw).decode()
 if 'table_name = heasarc_fermilgrb' not in text:raise ValueError('Wrong catalogue / Catálogo incorrecto')
 header,body=text.split('<DATA>',1)
 names=re.findall(r'^field\[([^]]+)\]',header,re.M)
 required={'gcn_name','trigger_met','redshift','tl0','irfs'}
 if not required.issubset(names) or '4FGL' in names:raise ValueError('Invalid GRB schema / Esquema GRB inválido')
 rows=[]
 for line in body.split('<END>')[0].split('</DATA>')[0].splitlines():
  if not line.strip() or line.startswith('#'):continue
  values=line.removesuffix('|').split('|')
  if len(values)!=len(names):raise ValueError(f'Column mismatch {len(values)} != {len(names)}')
  rows.append(values)
 expected=int(re.search(r'TOTAL ROWS:\s*(\d+)',header).group(1))
 if len(rows)!=expected:raise ValueError('Incomplete catalogue / Catálogo incompleto')
 df=pd.DataFrame(rows,columns=names)
 for name in ['trigger_met','redshift','tl0','ra','dec']:df[name]=pd.to_numeric(df[name],errors='coerce')
 return df,header

def K(z,n=1,omega_m=.315):
 if not np.isfinite(z) or z<0 or n not in (1,2):raise ValueError('Invalid z or order')
 return quad(lambda u:(1+u)**n/np.sqrt(omega_m*(1+u)**3+1-omega_m),0,z,epsabs=1e-12,epsrel=1e-12)[0]

def liv_delay(E_high,E_low,z,E_QG,H0_s,n=1,sign=1):
 # Observed energies in the same unit. H0 is an explicitly supplied cosmological parameter.
 if min(E_high,E_low)<0 or E_QG<=0 or H0_s<=0 or sign not in (-1,1):raise ValueError('Invalid physical parameter')
 return sign*(n+1)/(2*H0_s)*(E_high**n-E_low**n)/E_QG**n*K(z,n)

def main():
 raw=ROOT/'raw/heasarc_fermilgrb.tdat.gz'
 if not raw.exists():
  import requests
  r=requests.get(SOURCE,timeout=60);r.raise_for_status();raw.parent.mkdir(parents=True,exist_ok=True);raw.write_bytes(r.content)
 df,header=parse_catalog(raw.read_bytes());out=ROOT/'results';out.mkdir(exist_ok=True)
 df.to_csv(out/'NASA_HEASARC_FERMILGRB.csv',index=False)
 good=np.isfinite(df.redshift)&(df.redshift>0)&np.isfinite(df.tl0)&(df.tl0>0)
 selected=df.loc[good,['gcn_name','redshift','tl0','trigger_met','ra','dec','irfs']].copy()
 selected['K1_dimensionless']=[K(z) for z in selected.redshift]
 selected.to_csv(out/'catalog_onset_redshift_selected.csv',index=False)
 cmap=LinearSegmentedColormap.from_list('violet_apricot',['#6043A3','#AD67A5','#EDAA78'])
 for lang in ['ES','EN']:
  es=lang=='ES'
  with plt.rc_context({'figure.facecolor':'white','axes.facecolor':'white','text.color':'#172B4D','axes.labelcolor':'#172B4D','xtick.color':'#172B4D','ytick.color':'#172B4D','font.size':11}):
   fig,ax=plt.subplots(figsize=(10,6));sc=ax.scatter(selected.K1_dimensionless,selected.tl0,c=selected.tl0,cmap=cmap,norm=LogNorm(),s=55,edgecolors='white',linewidths=.6)
   ax.set_yscale('log');ax.set_xlabel(r'$K_1(z)=\int_0^z (1+u)/\sqrt{\Omega_m(1+u)^3+1-\Omega_m}\,du$')
   ax.set_ylabel('Inicio estimado de emisión LAT, TL0 (s)' if es else 'Estimated LAT emission onset, TL0 (s)')
   ax.set_title('Fermi-LAT · inicio de emisión y corrimiento al rojo' if es else 'Fermi-LAT · emission onset and redshift',pad=16)
   ax.grid(alpha=.15,color='#7788A5');fig.colorbar(sc,ax=ax,label='TL0 (s)')
   for _,row in selected.nlargest(4,'redshift').iterrows():ax.annotate(row.gcn_name,(row.K1_dimensionless,row.tl0),xytext=(-6,-18) if row.K1_dimensionless==selected.K1_dimensionless.max() else ((-6,8) if row.K1_dimensionless>selected.K1_dimensionless.max()*.9 else (5,8)),ha='right' if row.K1_dimensionless>selected.K1_dimensionless.max()*.9 else 'left',textcoords='offset points',fontsize=8,color='#895380')
   fig.text(.07,.025,'Catálogo observado; TL0 no es un retardo espectral. Sin ajuste LIV ni confirmación de RG.' if es else 'Observed catalogue; TL0 is not a spectral lag. No LIV fit or confirmation of GR.',fontsize=10,color='#5C5773')
   fig.tight_layout(rect=(0,.065,1,1));fig.savefig(out/f'Fermi_catalog_{lang}.png',dpi=170);fig.savefig(out/f'Fermi_catalog_{lang}.pdf');plt.close(fig)
 checks={'published_schema':True,'row_count':True,'finite_K':bool(np.isfinite(selected.K1_dimensionless).all()),'K_zero':K(0)==0,'equal_energy_zero_delay':liv_delay(1,1,1,1e19,2e-18)==0,'no_synthetic_lags_used':True}
 try:parse_catalog(gzip.compress(b'table_name = 4FGL\n<DATA>\n'))
 except ValueError:checks['reject_wrong_catalog']=True
 else:checks['reject_wrong_catalog']=False
 result={'source':SOURCE,'sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),'catalog_rows':len(df),'redshift_positive':int((df.redshift>0).sum()),'plotted_positive_onset_and_redshift':len(selected),'cosmology':{'omega_m':.315,'flat':True,'scope':'conditioning parameter, not an SI constant'},'checks':checks,'status':'DESCRIPTIVE_ONLY','notebooks_originals_changed':False,'non_claims':['No raw-photon reduction','No spectral-lag measurement','No LIV or GR confirmation','No synthetic values used as observations']}
 (out/'catalog_audit.json').write_text(json.dumps(result,indent=2)+'\n');assert all(checks.values());print(json.dumps(result,indent=2))
if __name__=='__main__':main()
