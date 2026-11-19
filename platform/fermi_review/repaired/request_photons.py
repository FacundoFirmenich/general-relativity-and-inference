"""Submit reproducible FSSC queries. Preserve response and do not infer from failure."""
from pathlib import Path
import json,re,requests
import pandas as pd
from astropy.time import Time
R=Path(__file__).resolve().parents[1]
targets=json.loads((R/'raw/archived_targets.json').read_text())
catalog=pd.read_csv(R/'results/NASA_HEASARC_FERMILGRB.csv')
rows=[]
for name,item in targets.items():
 match=catalog[catalog['name'].str.startswith(name[:9])]
 if len(match)==1:
  row=match.iloc[0];met=float(row.trigger_met);ra=float(row.ra);dec=float(row.dec);origin='HEASARC FERMILGRB'
 else:
  met=float((Time(item['trigger_time'],scale='utc')-Time('2001-01-01T00:00:00',scale='utc')).sec);ra=item['ra'];dec=item['dec'];origin='archived notebook target; independent trigger verification pending'
 payload={'coordfield':f'{ra},{dec}','coordsystem':'J2000','shapefield':'12','timefield':f'{met-100},{met+item["duration"]}','timetype':'MET','energyfield':'100,500000','photonOrExtendedOrNone':'Photon','spacecraft':'on','destination':'query'}
 try:
  response=requests.post('https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi',data=payload,timeout=45)
  (R/'raw'/f'{name}_query_response.html').write_text(response.text)
  urls=re.findall(r'https://fermi\.gsfc\.nasa\.gov/[^\s\x22<>]+',response.text)
  item_out={'grb':name,'payload':payload,'target_source':origin,'status':response.status_code,'response_urls':urls,'state':'SERVER_MAINTENANCE' if 'down for maintenance' in response.text else ('QUERY_SUBMITTED' if urls else 'NO_QUERY_RESULT')}
 except Exception as e:item_out={'grb':name,'payload':payload,'target_source':origin,'error':str(e)}
 rows.append(item_out);(R/'raw/new_queries.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(item_out),flush=True)
