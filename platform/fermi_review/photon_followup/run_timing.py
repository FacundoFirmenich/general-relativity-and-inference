"""Auditable SOURCE-event screening and descriptive timing, not a LIV likelihood."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];HERE=Path(__file__).resolve().parent
OUT=HERE/'results';OUT.mkdir(exist_ok=True)
P=json.loads((HERE/'protocol.json').read_text())
meta={x['target']:x for x in json.loads((R/'results/weekly_recovery.json').read_text())}
scmeta={x['target']:x for x in json.loads((R/'raw/spacecraft/manifest.json').read_text())}
def bitmask(bits):
 b=np.asarray(bits)
 if b.ndim!=2 or b.shape[1]!=32:raise ValueError('Expected FITS 32X field')
 return (b.astype(np.uint64)*(np.uint64(1)<<np.arange(31,-1,-1,dtype=np.uint64))).sum(axis=1)
def inside(t,starts,stops):
 idx=np.searchsorted(starts,t,side='right')-1
 safe=np.clip(idx,0,len(starts)-1)
 return (idx>=0)&(t>=starts[safe])&(t<stops[safe]),safe
# Independent orientation witness: FITS 32X is most-significant-bit first.
w=np.zeros((2,32),bool);w[0,24]=True;w[1,31]=True
assert np.array_equal(bitmask(w),[128,1])
rows=[];cuts=[];live=[];frames={};sources=[];rng=np.random.default_rng(P['bootstrap_seed'])
for name in P['targets']:
 m=meta[name];eventfile=R/'results/event_subsets'/f'{name}_audit_subset.fits';scfile=R/scmeta[name]['path']
 assert hashlib.sha256(eventfile.read_bytes()).hexdigest()==m['subset_sha256']
 assert hashlib.sha256(scfile.read_bytes()).hexdigest()==scmeta[name]['sha256']
 with fits.open(eventfile,checksum=True) as eh,fits.open(scfile,checksum=True) as sh:
  eh.verify('exception');sh.verify('exception')
  e=eh['EVENTS'].data;s=sh['SC_DATA'].data;g=eh['GTI'].data
  assert eh['EVENTS'].header['TIMESYS']==sh['SC_DATA'].header['TIMESYS']=='TT'
  assert np.all(s['STOP']>s['START']) and np.all(s['START'][1:]>=s['STOP'][:-1])
  assert np.all(g['STOP']>g['START']) and np.all(g['START'][1:]>=g['STOP'][:-1])
  t=np.array(e['TIME']);dt=t-m['trigger_MET'];energy=np.array(e['ENERGY']);classmask=bitmask(e['EVENT_CLASS']);typemask=bitmask(e['EVENT_TYPE'])
  ingti,_=inside(t,g['START'],g['STOP']);insc,ix=inside(t,s['START'],s['STOP'])
  quality=insc&(s['DATA_QUAL'][ix]>0)&(s['LAT_CONFIG'][ix]==1)
  dist=SkyCoord(m['ra_deg']*u.deg,m['dec_deg']*u.deg).separation(SkyCoord(e['RA']*u.deg,e['DEC']*u.deg)).deg
  masks=[('time_0_1000', (dt>=0)&(dt<1000)),('energy_100_100000_MeV',(energy>=100)&(energy<=100000)),('SOURCE_class',(classmask&P['event_class_mask'])!=0),('FRONT_or_BACK',(typemask&P['event_type_mask'])!=0),('ROI_3deg',dist<=P['ROI_deg']),('zenith_90deg',e['ZENITH_ANGLE']<P['zenith_max_deg']),('photon_GTI',ingti),('spacecraft_quality',quality)]
  keep=np.ones(len(e),bool)
  for label,mask in masks:
   keep&=mask;cuts.append({'target':name,'cut':label,'remaining':int(keep.sum())})
  assert np.all(insc[keep]) and np.all(ingti[keep])
  ids=np.column_stack([e['RUN_ID'][keep],e['EVENT_ID'][keep]])
  assert len(np.unique(ids,axis=0))==len(ids)
  provenance=pd.read_csv(R/'results/event_subsets'/f'{name}_audit_subset.csv')
  assert np.allclose(provenance.time_MET_s,t,rtol=0,atol=1e-7)
  frame=pd.DataFrame({'source_weekly_row_zero_based':provenance.source_row_zero_based.to_numpy()[keep],'run_id':e['RUN_ID'][keep],'event_id':e['EVENT_ID'][keep],'time_MET_s':t[keep],'dt_s':dt[keep],'energy_MeV':energy[keep],'separation_deg':dist[keep],'zenith_deg':e['ZENITH_ANGLE'][keep],'event_class_mask':classmask[keep],'event_type_mask':typemask[keep],'spacecraft_row':ix[keep]})
  frame=pd.DataFrame({col:frame[col].to_numpy().astype(frame[col].dtype.newbyteorder('='),copy=False) for col in frame.columns})
  frame.to_csv(OUT/f'{name}_screened_events.csv',index=False);frames[name]=frame
  # Clip spacecraft rows to each analysis window and photon GTI. LIVETIME/duration
  # is assumed uniform within each SC row; this is not effective area/exposure.
  scq=(s['DATA_QUAL']>0)&(s['LAT_CONFIG']==1)
  for end in P['sensitivity_end_s']:
   wall=0.;lt=0.
   for a,b in zip(g['START'],g['STOP']):
    length=np.maximum(0,np.minimum(s['STOP'],min(b,m['trigger_MET']+end))-np.maximum(s['START'],max(a,m['trigger_MET'])))
    wall+=float(length[scq].sum());lt+=float((length[scq]*s['LIVETIME'][scq]/(s['STOP'][scq]-s['START'][scq])).sum())
   assert -1e-8<=lt<=wall+1e-6<=end+1e-5
   live.append({'target':name,'window_end_s':end,'quality_GTI_wall_s':wall,'livetime_approx_s':lt,'method':'FT2 fractional-row approximation; not effective exposure'})
   sub=frame[frame.dt_s<end];a=sub.loc[sub.energy_MeV<1000,'dt_s'].to_numpy();b=sub.loc[sub.energy_MeV>=1000,'dt_s'].to_numpy()
   result={'target':name,'window_end_s':end,'n_low':len(a),'n_high':len(b),'median_high_minus_low_s':None,'bootstrap_95_low_s':None,'bootstrap_95_high_s':None,'status':'INSUFFICIENT_EVENTS'}
   if min(len(a),len(b))>=P['minimum_events_per_band']:
    draws=np.median(rng.choice(b,(P['bootstrap_draws'],len(b)),replace=True),axis=1)-np.median(rng.choice(a,(P['bootstrap_draws'],len(a)),replace=True),axis=1)
    lo,hi=np.quantile(draws,[.025,.975]);result.update(median_high_minus_low_s=float(np.median(b)-np.median(a)),bootstrap_95_low_s=float(lo),bootstrap_95_high_s=float(hi),status='DESCRIPTIVE_CONDITIONAL_ONLY')
    np.save(OUT/f'{name}_{end}s_bootstrap.npy',draws)
   rows.append(result)
  sources.append({'target':name,'event_subset_sha256':m['subset_sha256'],'spacecraft_sha256':scmeta[name]['sha256'],'weekly_source':m['source_url'],'weekly_source_sha256':m['source_sha256']})
for filename,data in [('timing_contrasts',rows),('cutflow',cuts),('livetime',live)]:pd.DataFrame(data).to_csv(OUT/(filename+'.csv'),index=False)
for lang in ['ES','EN']:
 es=lang=='ES'
 with plt.rc_context({'font.size':10,'figure.facecolor':'white','axes.facecolor':'white','text.color':'#172B4D','axes.labelcolor':'#172B4D','xtick.color':'#172B4D','ytick.color':'#172B4D'}):
  fig,axes=plt.subplots(2,2,figsize=(11,7),sharex=True,sharey=True)
  for ax,name in zip(axes.flat,P['targets']):
   f=frames[name];f=f[f.dt_s<100]
   for high,color,label in [(False,'#7652AC','0.1–1 GeV'),(True,'#EAA571','1–100 GeV')]:
    x=f[(f.energy_MeV>=1000)==high];ax.scatter(x.dt_s,x.energy_MeV/1000,s=25,color=color,label=f'{label} (N={len(x)})',alpha=.85,edgecolors='white',linewidths=.3)
   ax.set_yscale('log');ax.set_xlim(0,100);ax.grid(alpha=.16);ax.set_title(name);ax.legend(fontsize=8)
  for ax in axes[-1,:]:ax.set_xlabel('Tiempo desde el trigger (s)' if es else 'Time since trigger (s)')
  for ax in axes[:,0]:ax.set_ylabel('Energía reconstruida (GeV)' if es else 'Reconstructed energy (GeV)')
  fig.suptitle('Fermi-LAT · eventos seleccionados por calidad e instrumento' if es else 'Fermi-LAT · quality- and instrument-screened events',fontsize=15,y=.98)
  fig.text(.06,.02,'Selección SOURCE; incluye posible fondo. No es una medición de c ni una cota LIV.' if es else 'SOURCE selection; possible background remains. Not a measurement of c or a LIV bound.',fontsize=10,color='#635379')
  fig.tight_layout(rect=(0,.05,1,.95));fig.savefig(OUT/f'photon_timing_{lang}.png',dpi=170);fig.savefig(OUT/f'photon_timing_{lang}.pdf');plt.close(fig)
result={'protocol_sha256':hashlib.sha256((HERE/'protocol.json').read_bytes()).hexdigest(),'sources':sources,'contrasts':rows,'checks':{'FITS_bit_orientation':True,'input_hashes_verified':True,'SC_GTI_intervals_ordered':True,'event_row_provenance':True,'selected_unique_event_ids':True,'livetime_not_exceed_walltime':True},'status':'DESCRIPTIVE_TIMING_EXECUTED','physical_inference':'NOT_IDENTIFIED','reason':'Unmodelled intrinsic emission, energy-dependent effective area/PSF/dispersion and background; selected events are not labelled as source photons.'}
(OUT/'execution.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({'contrasts':rows,'status':result['status']},indent=2))
