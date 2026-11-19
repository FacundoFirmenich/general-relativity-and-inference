"""Execute the scoped table reanalysis and numerical regression checks."""
from pathlib import Path
import sys,subprocess,json,hashlib,datetime,platform,os,time
from importlib.metadata import version
R=Path(__file__).resolve().parent/'dcae_frontier'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
inputs=sorted([*R.glob('data/*.csv'),*R.glob('src/*.py'),*R.glob('scripts/*.py'),*R.glob('tests/*.py'),R/'records/source_audit.json'])
before={str(p.relative_to(R)):sha(p) for p in inputs}
start=time.time();stamp=datetime.datetime.now(datetime.timezone.utc).isoformat();commands=[]
env=os.environ.copy();env.update(MPLBACKEND='Agg')
for name,args in [('analysis',['scripts/run_analysis.py']),('tests',['-m','unittest','discover','-s','tests','-v']),('figures',['scripts/plot_integrated.py'])]:
 p=subprocess.run([sys.executable,*args],cwd=R,env=env,capture_output=True,text=True)
 (R/f'records/{name}.log').write_text(p.stdout+p.stderr)
 commands.append({'argv':[sys.executable,*args],'cwd':'research/dcae_frontier','exit_status':p.returncode,'log':f'records/{name}.log'})
 if p.returncode:raise RuntimeError(p.stdout+p.stderr)
after={str(p.relative_to(R)):sha(p) for p in inputs};assert before==after
outputs=[R/'results/REANALYSIS.json',*R.glob('figures/cosmos_*')]
m={'schema':'dcae.frontier.execution.v1','scope':'Published table reanalysis, 16 finite tests and publication figures; no raw-catalogue Monte Carlo or empirical FKP execution','started_at':stamp,'runtime_seconds':time.time()-start,'status':'completed','commands':commands,'environment':{'python':platform.python_version(),'numpy':version('numpy'),'matplotlib':version('matplotlib')},'inputs':before,'inputs_unchanged':before==after,'outputs':{str(p.relative_to(R)):sha(p) for p in outputs},'new_independent_measurements':0,'physical_DCAE':'NOT_ESTIMABLE'}
(R/'records/execution.json').write_text(json.dumps(m,indent=2)+'\n')
print('Frontier: reanalysis, 16 tests, bilingual figures passed; input hashes unchanged.')
