from pathlib import Path
import requests,re,json,hashlib,concurrent.futures
R=Path(__file__).resolve().parents[1];out=R/'raw/spacecraft';out.mkdir(exist_ok=True)
base='https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly/spacecraft/'
s=requests.get(base,timeout=40);s.raise_for_status();(out/'index.html').write_text(s.text)
rows=json.loads((R/'results/weekly_recovery.json').read_text())
def get(x):
 week=re.search(r'_w(\d+)_',x['source_url']).group(1)
 names=sorted(set(re.findall(r'href="(lat_spacecraft_weekly_w'+week+r'_[^"]+\.fits)"',s.text)))
 if len(names)!=1:raise ValueError(names)
 name=names[0];p=out/name
 if not p.exists():
  response=requests.get(base+name,timeout=90);response.raise_for_status();p.write_bytes(response.content)
 result={'target':x['target'],'url':base+name,'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()};print(json.dumps(result),flush=True);return result
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:result=list(pool.map(get,rows))
(out/'manifest.json').write_text(json.dumps(result,indent=2)+'\n')
