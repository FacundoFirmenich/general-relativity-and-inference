"""Reject private manuscript paths and embedded manuscript archives."""
from pathlib import Path
import json,re,base64,zipfile,io
R=Path(__file__).resolve().parents[1]
def forbidden(name):
 p=Path(name)
 return 'publication' in p.parts or p.name.startswith(('PAPER_','SUPPLEMENTARY_','LIVING_SI'))
fail=[];payloads=0
for p in R.rglob('*'):
 if not p.is_file() or any(x.startswith('.') for x in p.relative_to(R).parts):continue
 if forbidden(p.relative_to(R)):fail.append(str(p))
 if p.name.startswith('GRAVITAS_') and p.suffix in ['.py','.ipynb']:
  text=p.read_text()
  if p.suffix=='.ipynb':text='\n'.join(''.join(c.get('source',[])) for c in json.loads(text)['cells'])
  m=re.search(r"^PAYLOAD='([^']+)'",text,re.M)
  if m:
   payloads+=1
   with zipfile.ZipFile(io.BytesIO(base64.b64decode(m[1]))) as z:
    fail.extend(x for x in z.namelist() if forbidden(x))
assert not fail,fail
print(json.dumps({'manuscripts_absent':True,'embedded_payloads_checked':payloads}))
