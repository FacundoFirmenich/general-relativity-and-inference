"""Download bounded, public Gaia/JPL source records and record exact provenance."""
import urllib.request,urllib.parse,json,hashlib,datetime,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def get(url,name,params=None):
    target=ROOT/'raw'/name
    if target.exists(): return
    if params: url+='?'+urllib.parse.urlencode(params)
    request=urllib.request.Request(url,headers={'User-Agent':'Yukawa-observation-audit/1.0 (scientific-reproducibility)'})
    start=time.time()
    try:
        with urllib.request.urlopen(request,timeout=90) as response:
            data=response.read(15_000_001)
            if len(data)>15_000_000:raise ValueError('bounded download exceeded')
            record={'url':url,'status':response.status,'mime':response.headers.get('Content-Type')}
        target.write_bytes(data)
        record.update(path=str(target.relative_to(ROOT)),bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
    except Exception as e:record={'url':url,'error':str(e)}
    record.update(retrieved_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),seconds=time.time()-start)
    with (ROOT/'downloads.jsonl').open('a') as f:f.write(json.dumps(record)+'\n')
    print(json.dumps(record),flush=True)

if __name__=='__main__':
    for number in ['1566','161989','21','1917']:
        if number in ['1566','161989']:
            get('https://ssd-api.jpl.nasa.gov/sb_radar.api',f'radar_{number}.json',
                {'des':number,'type':'R','bp':'C','coords':'true','ref':'true','notes':'true'})
        get('https://gea.esac.esa.int/tap-server/tap/sync',f'gaia_{number}.csv',
            {'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv',
             'QUERY':f'SELECT * FROM gaiafpr.sso_observation WHERE number_mp={number} ORDER BY epoch'})
