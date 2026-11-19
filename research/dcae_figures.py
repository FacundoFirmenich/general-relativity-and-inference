"""Copy the executed DCAE figure artifacts into the publication figure tree."""
from pathlib import Path
import hashlib,json,shutil
R=Path(__file__).resolve().parents[1]
P=R/'research/dcae'
m=json.loads((P/'records/computation_manifest.json').read_text())
assert m['run']['status']=='completed'
recorded={x['path']:x['sha256'] for x in m['outputs']}
for lang in ['en','es']:
    for name in ['diagnostics','repair']:
        rel=f'figures/{name}_{lang.upper()}.pdf';src=P/rel
        assert hashlib.sha256(src.read_bytes()).hexdigest()==recorded[rel],rel
        shutil.copy2(src,R/f'figures/dcae_{name}_{lang}.pdf')
print('Copied four hash-checked DCAE figure artifacts.')
