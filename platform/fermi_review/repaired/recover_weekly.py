"""Recover archived LAT events. Geometric/time selection only; NOT a likelihood-ready sample."""
from pathlib import Path
import requests,re,json,hashlib,math,concurrent.futures
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
R=Path(__file__).resolve().parents[1]
base='https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly/photon/'
cache=R/'raw/weekly';cache.mkdir(exist_ok=True)
out=R/'results/event_subsets';out.mkdir(exist_ok=True)
catalog=pd.read_csv(R/'results/NASA_HEASARC_FERMILGRB.csv')
targets=json.loads((R/'raw/archived_targets.json').read_text())
index=(R/'raw/weekly_photon_index.html').read_text()
# Week 9 ends at this MET, as verified from its official FITS header.
week9stop=239760269.
jobs=[]
for name,old in targets.items():
 match=catalog[catalog['name'].str.startswith(name[:9])]
 if len(match)!=1:continue
 row=match.iloc[0];met=float(row.trigger_met)
 week=10+int((met-week9stop)//604800)
 filenames=sorted(set(re.findall(r'href="(lat_photon_weekly_w%03d_[^"]+\.fits)"'%week,index)))
 if len(filenames)!=1:raise RuntimeError((name,week,filenames))
 jobs.append((name,row,old,filenames[0]))
def work(job):
 name,row,old,filename=job;path=cache/filename;url=base+filename
 if not path.exists():
  response=requests.get(url,timeout=(20,90));response.raise_for_status()
  if response.status_code!=200:raise RuntimeError('Partial archive response')
  path.write_bytes(response.content)
 digest=hashlib.sha256(path.read_bytes()).hexdigest()
 met=float(row.trigger_met);lo=met-100;hi=met+old['duration']
 with fits.open(path,checksum=True) as hd:
  hd.verify('exception');events=hd['EVENTS'].data
  assert hd['EVENTS'].header['TSTART']<=lo<hi<=hd['EVENTS'].header['TSTOP']
  center=SkyCoord(float(row.ra)*u.deg,float(row.dec)*u.deg)
  first=np.flatnonzero((events['TIME']>=lo)&(events['TIME']<=hi)&(events['ENERGY']>=100)&(events['ENERGY']<=500000))
  dist=center.separation(SkyCoord(events['RA'][first]*u.deg,events['DEC'][first]*u.deg)).deg
  selected=first[dist<=12.]
  sub=fits.BinTableHDU(data=events[selected],header=hd['EVENTS'].header,name='EVENTS')
  sub.header['TSTART']=lo;sub.header['TSTOP']=hi
  sub.header.add_history('Audit subset: 12 deg, 100-500000 MeV, trigger-100 s to trigger+duration.')
  sub.header.add_history('No claim of IRF, exposure, zenith, quality or background validation.')
  gti=hd['GTI'].data;gkeep=(gti['STOP']>lo)&(gti['START']<hi)
  g=gti[gkeep].copy();g['START']=np.maximum(g['START'],lo);g['STOP']=np.minimum(g['STOP'],hi)
  gh=fits.BinTableHDU(data=g,header=hd['GTI'].header,name='GTI')
  file=out/(name+'_audit_subset.fits');fits.HDUList([hd[0].copy(),sub,gh]).writeto(file,overwrite=True,checksum=True)
  pd.DataFrame({'source_row_zero_based':selected,'time_MET_s':events['TIME'][selected],'dt_trigger_s':events['TIME'][selected]-met,'energy_MeV':events['ENERGY'][selected],'ra_deg':events['RA'][selected],'dec_deg':events['DEC'][selected]}).to_csv(out/(name+'_audit_subset.csv'),index=False)
  inside_gti=np.zeros(len(selected),bool)
  for start,stop in zip(g['START'],g['STOP']):inside_gti|=(events['TIME'][selected]>=start)&(events['TIME'][selected]<stop)
  result={'target':name,'catalog_name':row['name'],'trigger_MET':met,'ra_deg':float(row.ra),'dec_deg':float(row.dec),'source_url':url,'source_sha256':digest,'source_bytes':path.stat().st_size,'source_rows':len(events),'subset_rows_12deg':len(selected),'subset_rows_in_GTI':int(inside_gti.sum()),'subset_sha256':hashlib.sha256(file.read_bytes()).hexdigest(),'GTI_duration_s':float(np.sum(g['STOP']-g['START'])),'status':'ARCHIVE_RECOVERED_NOT_INFERENCE_READY','limitations':['No spacecraft-quality selection','No zenith cut','No explicit event class/type selection','No exposure or IRF computation','No background or intrinsic-delay inference'],'scope':'Not a reproduction of the original 3-degree extraction: catalogue coordinates, time window and archive processing may differ.'}
  print(json.dumps(result),flush=True);return result
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(work,jobs))
(R/'results/weekly_recovery.json').write_text(json.dumps(results,indent=2)+'\n')
print('Recovered',len(results),'targets; no GR/LIV bound inferred.')
