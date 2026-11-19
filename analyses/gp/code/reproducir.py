"""Execute the finite experiment contract and record hashes before and after."""
from pathlib import Path
import sys,subprocess,os,json,hashlib,time,platform,datetime
import numpy,scipy,mpmath,matplotlib
R=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
inputs=sorted(list((R/'sources/cassini').glob('*'))+list((R/'code').glob('*.py'))+
    list((R/'sources').glob('*source.tex'))+[R/'records/CONTRATO.md'])
before=[{'path':str(p.relative_to(R)),'sha256':sha(p)} for p in inputs if p.is_file()]
scripts=['gp_covarianza_cassini.py','gp_emulador.py','gp_realizaciones.py','numeros_y_unidades.py','verificar_calculos.py','figuras.py']
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
t=time.time();start=datetime.datetime.now(datetime.timezone.utc).isoformat();commands=[]
for name in scripts:
    cmd=[sys.executable,'code/'+name];begin=time.time()
    p=subprocess.run(cmd,cwd=R,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    (R/'records'/name.replace('.py','.log')).write_text(p.stdout)
    commands.append({'argv':cmd,'exit_status':p.returncode,'runtime_seconds':time.time()-begin})
    print(name,p.returncode,round(time.time()-begin,3),flush=True)
    if p.returncode:break
for row in before:row['sha256_after']=sha(R/row['path'])
ok=all(c['exit_status']==0 for c in commands) and len(commands)==len(scripts)
outputs=[{'path':str(p.relative_to(R)),'sha256':sha(p),'kind':'result'} for folder in ['results','figures'] for p in sorted((R/folder).glob('*')) if p.is_file()]
checks=json.loads((R/'results/independent_checks.json').read_text()) if ok else {}
manifest={'schema_version':2,'claim_id':'GP-YUKAWA-2026-11-17',
 'repository':{'commit':None,'dirty':True,'reason':'Standalone analysis directory, not a Git checkout; code identity is recorded by SHA256'},
 'command':commands,'environment':{'software':[{'name':name,'version':version} for name,version in [('Python',platform.python_version()),('numpy',numpy.__version__),('scipy',scipy.__version__),('mpmath',mpmath.__version__),('matplotlib',matplotlib.__version__)]],'hardware':platform.platform(),'variables':{k:env[k] for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']}},
 'mathematics':{'assertion_tested':'Joint GP law, emulator boundary fidelity, fixed-kernel injection response, and numerical conversion consistency',
 'coefficient_domain':'Real arithmetic approximated by float64; independent 60-digit check for sparse GP predictions',
 'conventions':'Dimensionless diagnostic x is not graviton mass; optical A is not gamma; physical length uses hbar unless explicitly contrasted with h',
 'inputs':[r['path'] for r in before],
 'bounds':{'diagnostic_x':[0,3],'training_nodes':[4,7,13],'realization_grids':[51,101,201,401,801],'simulations':4096,'cassini_epochs':1696},
 'non_claims':['No recalculated INPOP mass bound','No proof of independence in unpublished GPUR implementation','No calibrated instrumental covariance from sensitivity kernels','No universal invalidity of GPR']},
 'randomness':{'used':True,'generator':'numpy.default_rng (PCG64)','seed':20261117},
 'run':{'started_at':start,'runtime_seconds':time.time()-t,'exit_status':0 if ok else 1,'status':'completed' if ok else 'failed','resource_limits':{'memory_bytes':None,'cpu_seconds':None,'max_cores':None,'max_threads':1}},
 'outputs':outputs,'checks':checks,'result':'All executed assertions passed' if ok else 'Failure recorded',
 'residual_risks':['Self-review, no independent reviewer','Published INPOP fitting nodes and executable implementation were not recovered','Conditional Cassini reconstruction inherited from the archived prior analysis; instrument covariance not calibrated'],
 'input_artifacts':before,'declared_results':['results/gp_cassini_refits.json','results/gp_emulator.json','results/gp_realizations.json','results/compton_consistency.json','results/number_registry.csv','results/independent_checks.json']}
assert all(r['sha256']==r['sha256_after'] for r in before),'Input changed during run'
(R/'records/execution_manifest_v2.json').write_text(json.dumps(manifest,indent=2)+'\n')
if not ok:raise SystemExit(1)
