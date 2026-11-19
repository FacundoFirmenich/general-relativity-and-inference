"""Paired SI sensitivity calculation; archived reproduction remains unchanged."""
from pathlib import Path
import importlib.util, json, hashlib, math
ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'analyses/casimir/code/auditar_casimir.py'
spec=importlib.util.spec_from_file_location('casimir_archived',p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
old={k:getattr(m,k) for k in ['C','HBAR','EV','KB','WP','GAMMA']}
# h, e, k_B and c have exact SI definitions; hbar is h/(2*pi), evaluated numerically.
new={'C':299792458.0,'HBAR':6.62607015e-34/(2*math.pi),'EV':1.602176634e-19,'KB':1.380649e-23}
new['WP']=13*new['EV']/new['HBAR'];new['GAMMA']=.05*new['EV']/new['HBAR']
rows=[]
for nm in [10.,20.,40.,80.,110.]:
 for k,v in old.items():setattr(m,k,v)
 archived=m.pressure(nm*1e-9)
 for k,v in new.items():setattr(m,k,v)
 si=m.pressure(nm*1e-9)
 fine=m.pressure(nm*1e-9,ny=768,nu=512)
 rows.append(dict(gap_nm=nm,pressure_archived_constants_Pa=archived,pressure_SI_Pa=si,relative_SI_shift=si/archived-1,relative_quadrature_difference=si/fine-1))
result={'source':'https://www.bipm.org/en/publications/si-brochure','archived_code_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'old':old,'SI':new,'h_exact_J_s':'6.62607015e-34','scope_es':'Sensibilidad de presión Drude T=0, con parámetros materiales fijos; no es un ajuste experimental ni una cota Yukawa. hbar se deriva de h/(2 pi); no tiene una última cifra decimal finita.','scope_en':'Sensitivity of T=0 Drude pressure with fixed material parameters; not an experimental fit or a Yukawa bound. hbar is derived from h/(2 pi); it has no final finite decimal digit.','rows':rows}
assert max(abs(r['relative_quadrature_difference']) for r in rows)<1e-5
(ROOT/'constants_audit/results.json').write_text(json.dumps(result,indent=2)+'\n')
for lang in ['ES','EN']:
 es=lang=='ES'
 title='Auditoría de constantes SI' if es else 'SI constants audit'
 body='# '+title+'\n\n'+result['scope_es' if es else 'scope_en']+'\n\n'
 body+=('Se preserva el código archivado: sus constantes redondeadas permiten reproducir su tabla, pero no constituyen los valores SI actuales. Esta ejecución compara ambos conjuntos sin modificar el material original.\n\n' if es else 'Archived code is preserved: its rounded constants reproduce its table but are not current SI values. This execution compares both sets without modifying the original material.\n\n')
 body+='c = 299792458 m/s; h = 6.62607015e-34 J s; e = 1.602176634e-19 C; k_B = 1.380649e-23 J/K.\n\n'
 body+=('| Separación (nm) | Cambio relativo de presión (%) | Error relativo de cuadratura |\n' if es else '| Gap (nm) | Relative pressure shift (%) | Relative quadrature error |\n')+'|---:|---:|---:|\n'
 for r in rows:body+=f"| {r['gap_nm']:.0f} | {100*r['relative_SI_shift']:.8f} | {r['relative_quadrature_difference']:.3e} |\n"
 body+='\n'+('Reproducir: ' if es else 'Reproduce: ')+'`python constants_audit/audit_si.py`\n\n'+result['source']+'\n'
 (ROOT/f'constants_audit/README_{lang}.md').write_text(body)
print(json.dumps(result,indent=2))
