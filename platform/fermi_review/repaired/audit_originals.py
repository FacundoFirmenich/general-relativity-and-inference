"""Execute numerical checks on archived notebook definitions without running its data claims."""
from pathlib import Path
import json,ast,hashlib
import numpy as np
from catalog_audit import K
R=Path(__file__).resolve().parents[1]
n=json.loads((R/'originals/fermilat.LIV.ipynb').read_text());tree=ast.parse(''.join(n['cells'][0]['source']))
node=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='generate_GRB090510_mock_data')
scope={'np':np};exec(compile(ast.Module(body=[node],type_ignores=[]),'archived_mock','exec'),scope)
E,t,sigma,z=scope[node.name]();truth=2*E**(-.5);residual=t-truth
from scipy.integrate import quad
rows=[]
for z in [.340,.903,1.822,4.35]:
 old=(1+z)*quad(lambda u:1/np.sqrt(.315*(1+u)**3+.685),0,z)[0]
 rows.append({'z':z,'old_over_correct_linear':old/K(z,1),'old_over_correct_quadratic':old/(1.5*K(z,2))})
files=[{'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'cells':len(json.loads(p.read_text())['cells'])} for p in sorted((R/'originals').glob('fermi*.ipynb'))]
result={'originals':files,'mock_execution':{'n':len(E),'E_min_GeV':float(E.min()),'E_max_GeV':float(E.max()),'injected_noise_sigma_s':.1,'declared_likelihood_sigma_s':float(sigma[0]),'realized_noise_RMS_s':float(np.sqrt(np.mean(residual**2))),'scope':'synthetic only; not a NASA photon sample'},'LIV_kernel_check':rows,'beta8_archived_photon_counts':{'GRB090510':142,'GRB080916C':163,'GRB090902B':401,'GRB130427A':924,'GRB221009A':196},'beta8_narrative_photon_counts':{'GRB090510':196,'GRB080916C':240,'GRB090902B':179,'GRB130427A':475,'GRB221009A':255},'beta8_example':{'GRB090510_output_N0':.033,'GRB090510_narrative_N0':6.13e-10,'GRB090510_output_alpha':1.609,'GRB090510_narrative_alpha':2.45},'checks':{'mock_reproduced_50':len(E)==50,'noise_sigma_mismatch':bool(np.all(sigma==.05)),'kernel_mismatch_nonzero_z':all(abs(x['old_over_correct_linear']-1)>.01 for x in rows)},'scope':'Checks on archived code and outputs. Raw photon extraction and physical likelihood are not revalidated.'}
assert all(result['checks'].values());(R/'results/originals_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
