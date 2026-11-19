"""Export actual computed records as CSV; no manuscript generation."""
from pathlib import Path
import json,csv
R=Path(__file__).resolve().parents[1];out=R/'tables';out.mkdir(exist_ok=True)
for base in ['results','research/results','research/revision_v17/results','research/revision_v18/results','research/revision_v19/results','research/fractalgo/results']:
 for p in (R/base).glob('*.json'):
  rows=[]
  def walk(x,key=''):
   if isinstance(x,dict):
    for k,v in x.items():walk(v,key+'/'+str(k))
   elif isinstance(x,list):
    for i,v in enumerate(x):walk(v,key+'/'+str(i))
   else:rows.append([key,x])
  walk(json.loads(p.read_text()))
  with (out/('_'.join(p.relative_to(R).parts).replace('.json','.csv'))).open('w',newline='') as f:
   w=csv.writer(f);w.writerow(['field / campo','value / valor']);w.writerows(rows)
print('Tablas regeneradas / Tables regenerated:',len(list(out.glob('*.csv'))))
