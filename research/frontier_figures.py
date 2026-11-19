"""Copy verified frontier figures without silently rerunning analysis."""
from pathlib import Path
import json,hashlib,shutil
R=Path(__file__).resolve().parents[1];F=R/'research/dcae_frontier'
m=json.loads((F/'records/execution.json').read_text());assert m['status']=='completed'
for lang in ['es','en']:
 rel=f'figures/cosmos_{lang}.pdf';p=F/rel
 assert hashlib.sha256(p.read_bytes()).hexdigest()==m['outputs'][rel]
 shutil.copy2(p,R/f'figures/dcae_cosmos_{lang}.pdf')
print('Copied two hash-checked COSMOS figures.')
