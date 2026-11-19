#!/usr/bin/env python3
"""Single entry point. Numerical modes run without TeX, browser or network."""
from pathlib import Path
import os,sys,subprocess,json,hashlib,time,platform,datetime,argparse
from importlib.metadata import version,PackageNotFoundError
R=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['all','calculations','figures','tables','validate','documents','dcae','frontier','revision','graph','information','fermi','stability','topology'],nargs='?',default='all')
    args=parser.parse_args();mode=args.mode
    (R/'records').mkdir(exist_ok=True)
    (R/'figures').mkdir(exist_ok=True)
    plans={
      'calculations':[['analyses/gp/code/'+s] for s in ['gp_covarianza_cassini.py','gp_emulador.py','gp_realizaciones.py','numeros_y_unidades.py','verificar_calculos.py']]+[
        ['analyses/casimir/code/auditar_casimir.py'],['-m','gravitas_si.audits'],['-m','gravitas_si.inference'],['-m','gravitas_si.sensitivity'],['-m','gravitas_si.contrast_checks']]+[[f'research/{name}.py'] for name in ['gaia_covariance','inference_checks','physical_checks','mond_quadrupole','mond_diagnostics','cassini_odf_audit']],
      'figures':[['-c','from gravitas_si.visuals import main; main("es"); main("en")'],['-m','gravitas_si.record_visuals'],['research/figures.py'],['research/mond_figures.py']],
      'tables':[['tools/export_tables.py']],
      'validate':[['analyses/gp/code/verificar_calculos.py'],['-m','gravitas_si.audits'],['-m','gravitas_si.contrast_checks']],
      'documents':[['tools/private_documents.py']],
      'dcae':[['research/dcae/verify.py'],['research/dcae/validate_manifest.py','research/dcae/records/computation_manifest.json','--root','research/dcae']]}
    plans['revision']=[[f'research/revision_v17/{name}.py'] for name in ['audit','cassini','information']]
    plans['information']=[['research/revision_v18/information_audit.py']]
    plans['stability']=[['research/revision_v19/posterior_stability.py']]
    plans['topology']=[[f'research/fractalgo/{name}.py'] for name in ['test_core','run_cases','figures']]
    plans['fermi']=[[f'platform/fermi_review/{x}'] for x in ['repaired/catalog_audit.py','repaired/audit_originals.py','repaired/audit_hunter.py','photon_followup/run_timing.py']]
    plans['calculations'].extend(plans['stability']+plans['fermi']+plans['topology'][:2])
    plans['figures'].append(['research/fractalgo/figures.py'])
    plans['revision'].extend(plans['information'])
    plans['revision'].append(['platform/audit_gravitools.py'])
    plans['graph']=[['tools/build_epistemic_graph.py'],['tools/render_epistemic_print.py']]
    plans['calculations'].extend(plans['revision'])
    plans['figures'].append(['research/revision_v17/figures.py'])
    plans['figures'].extend(plans['graph'])
    plans['frontier']=[['research/frontier_run.py'],['research/frontier_figures.py']]
    plans['calculations'].extend([['research/dcae/verify.py'],['research/frontier_run.py']])
    plans['figures'].append(['research/dcae_figures.py'])
    plans['dcae'].append(['research/dcae/normalization_audit.py'])
    plans['calculations'].append(['research/dcae/normalization_audit.py'])
    plans['dcae'].append(['research/dcae_figures.py'])
    plans['dcae'].extend(plans['frontier'])
    plans['figures'].append(['research/frontier_figures.py'])
    plan=plans['calculations']+plans['tables']+plans['figures'] if mode=='all' else plans[mode]
    missing=[p for p in ['analyses/gp/sources/cassini/dss25_diseno_erfa_constant.npz','analyses/gp/sources/cassini/contraste_robusto_pesos.npz'] if not (R/p).exists()]
    if missing and mode not in ['documents']:
        raise SystemExit('Faltan entradas numericas / Numerical inputs missing: '+', '.join(missing)+'. Ver / See NUMERICAL_DATA_MANIFEST.json; usar el paquete completo / use the complete package.')
    if mode=='figures':plan=[['research/dcae/verify.py'],['research/frontier_run.py']]+plan
    # Document compilation is explicit: all numerical results do not require a TeX installation.
    files=[]
    for folder in ['gravitas_si','inputs','analyses','docs','publication']:
        files.extend(p for p in (R/folder).rglob('*') if p.is_file() and ('sources' in p.parts or p.suffix in ['.py','.h5','.md','.tex','.css']) and '__pycache__' not in p.parts and 'results' not in p.parts and 'figures' not in p.parts and p.name not in ['revision_supp.md','revision_main.md','information_supp.md','transformations_rendered.md','topology_supp.md','topology_main.md'])
    files += [R/'run.py',R/'requirements.txt']
    files += [p for p in (R/'research/dcae_frontier').rglob('*') if p.is_file() and ('data' in p.parts or p.suffix=='.py') and '__pycache__' not in p.parts and 'original' not in p.parts]
    files += list((R/'research').glob('*.py'))
    files += list((R/'research/revision_v17').glob('*.py'))+list((R/'research/revision_v18').glob('*.py'))+[R/'platform/audit_gravitools.py',R/'tools/build_epistemic_graph.py',R/'tools/render_epistemic_print.py',R/'ui/epistemic/template.html']
    files += list((R/'research/dcae').glob('*.py'))+[R/'research/dcae/CONTRACT.md']
    files += list((R/'research/raw').glob('gaia_*.csv'))+list((R/'research/raw').glob('radar_*.json'))+list((R/'research/raw').glob('sparc_*.mrt'))+list((R/'research/raw').glob('cassini_*'))
    if mode=='revision':
        files=[R/'run.py',R/'gravitas_si/inference.py',R/'analyses/gp/code/gp_covarianza_cassini.py',R/'research/physical_checks.py',R/'research/results/gaia_covariance_transits.csv']+list((R/'research/revision_v17').glob('*.py'))+[R/'analyses/gp/sources/cassini'/name for name in ['dss25_diseno_erfa_constant.npz','contraste_robusto_pesos.npz']]+list((R/'research/revision_v18').glob('*.py'))+[R/'platform/audit_gravitools.py',R/'platform/originals/GraviTools_original.py']
    files += list((R/'research/fractalgo').glob('*.py'))+[R/'research/fractalgo/CONTRACT.md']
    files += list((R/'research/revision_v19').glob('*.py')) + list((R/'platform/fermi_review/repaired').glob('*.py')) + list((R/'platform/fermi_review/photon_followup').glob('*.py'))
    files += [p for sub in ['raw/spacecraft','results/event_subsets','originals'] for p in (R/'platform/fermi_review'/sub).rglob('*') if p.is_file()]
    before=[dict(path=str(p.relative_to(R)),sha256=sha(p)) for p in sorted(set(files))]
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLBACKEND='Agg')
    start=time.time();stamp=datetime.datetime.now(datetime.timezone.utc).isoformat();commands=[]
    for number,cmd in enumerate(plan):
        t=time.time();full=[sys.executable,*cmd];p=subprocess.run(full,cwd=R,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        log=f'records/{mode}_{number:02d}.log';(R/log).write_text(p.stdout)
        commands.append(dict(argv=full,exit_status=p.returncode,log=log,runtime_seconds=time.time()-t))
        print(f'{number+1}/{len(plan)} {cmd[-1]}: exit {p.returncode}, {time.time()-t:.1f}s',flush=True)
        if p.returncode:
            print(p.stdout[-3000:],flush=True);break
    for row in before:row['sha256_after']=sha(R/row['path'])
    versions=[]
    for name in ['numpy','scipy','mpmath','matplotlib','h5py','dynesty','emcee','plotly']:
        try:versions.append(dict(name=name,version=version(name)))
        except PackageNotFoundError:versions.append(dict(name=name,version=None))
    success=all(x['exit_status']==0 for x in commands) and len(commands)==len(plan) and all(x['sha256']==x['sha256_after'] for x in before)
    outputs=[dict(path=str(p.relative_to(R)),sha256=sha(p)) for folder in ['results','figures','analyses/gp/results','analyses/casimir/results','research/results','research/dcae/results','research/dcae/figures','research/dcae/records','research/dcae_frontier/results','research/dcae_frontier/records','research/dcae_frontier/figures','research/revision_v17/results','research/revision_v18/results','ui/epistemic','platform/audit_runs','research/revision_v19/results','platform/fermi_review/results','platform/fermi_review/photon_followup/results','research/fractalgo/results'] for p in sorted((R/folder).rglob('*')) if p.is_file()]
    manifest=dict(schema_version=2,claim_id='GRAVITAS-SI-2026-11-18',mode=mode,
      repository=dict(commit=None,dirty=True,reason='Standalone release; file SHA256 identifies code'),command=commands,
      environment=dict(python=platform.python_version(),platform=platform.platform(),software=versions,variables={k:env[k] for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','MPLBACKEND']}),
      mathematics=dict(assertion_tested='Conditional GP, force/PPN benchmarks, hierarchical optical inference, DESI operator, Casimir normalization and finite synthetic DCAE controls',
        coefficient_domain='float64; independent 60-digit sparse GP checks',conventions='SI; metric signature -+++; A optical amplitude is not gamma or graviton mass',
        bounds=dict(cassini_epochs=1696,gp_realizations=4096,nested_live_points=250,mcmc_chains=4,mcmc_retained_per_chain=12000,sbc_replicates=4000,ppc_draws=256),
        non_claims=['No full INPOP refit','No new graviton bound','No external instrument calibration','No universal GP invalidity','No independent peer review']),
      randomness=dict(used=True,generator='NumPy PCG64 and scrambled SciPy Sobol',seeds='Explicit in each module: 20261117; 20261118/19; 20261200--20261400'),
      run=dict(started_at=stamp,runtime_seconds=time.time()-start,status='completed' if success else 'failed',exit_status=0 if success else 1),
      input_artifacts=before,outputs=outputs,
      checks=dict(inputs_unchanged=all(x['sha256']==x['sha256_after'] for x in before),commands_completed=len(commands)),
      result='Numerical execution passed; temporal posterior predictive discrepancy remains scientifically unresolved' if success else 'Execution failure; see logs',
      residual_risks=['Fixed temporal kernel and empirical baseline weights','Conditional reconstruction inherited','Priors are scenarios without independent calibration','Windows/macOS launchers not executed on those OSs; Colab cloud runtime not externally tested'])
    if mode in ['dcae','frontier']:
        manifest['mathematics']={'assertion_tested':'Scoped DCAE mathematics and/or published COSMOS table reanalysis, FKP floating-point and finite inclusion checks','bounds':{'frontier_unit_tests':16,'dcae_finite_checks':15 if mode=='dcae' else 0,'new_independent_observations':0},'non_claims':['No new physical DCAE measurement','No raw COSMOS recovery rerun','No empirical FKP estimate','No rerun of earlier gravity calculations']}
        manifest['randomness']={'used':mode=='dcae','generator':'PCG64' if mode=='dcae' else None,'seed':20261118 if mode=='dcae' else None,'scope':'Synthetic clipping checks only; published table reanalysis is deterministic'}
        manifest['result']='Scoped checks completed; physical DCAE remains NOT_ESTIMABLE' if success else 'Execution failed; see command logs'
    if mode=='revision':
        manifest['mathematics']={'assertion_tested':'Gaussian contrast distances, four-component conditional Cassini hierarchy with deleted-pass prediction, exact Gaussian decision retention and eccentricity sensitivity conversion','bounds':{'cassini_epochs':1696,'passes':20,'kernel_components':4,'ppc_draws_per_model':1024,'new_independent_observations':0},'non_claims':['No raw telemetry rerun','No all-parameter ephemeris refit','No physical alpha/lambda posterior','No full hierarchical SBC','No external peer review']}
        manifest['randomness']={'used':True,'generator':'NumPy PCG64','seeds':[2026111807,2026111818],'scope':'Cassini posterior predictive simulation and separate information-control calibration; numerical integration deterministic'}
        manifest['result']='Scoped revision executed; local Cassini discrepancies retained' if success else 'Execution failed; see logs'
    if mode in ['fermi','stability']:
        manifest['mathematics']={'assertion_tested':mode+' scoped execution','non_claims':['No new GR confirmation or refutation','No LIV propagation bound','No continuous emulator convergence claim']}
        manifest['result']='Scoped execution completed; physical interpretation remains bounded by documented assumptions' if success else 'Execution failed'
    if mode=='topology':
        manifest['mathematics']={'assertion_tested':'Fixed-graph prior-relative GP excursions and conditional Cassini residual topology','coefficient_domain':'float64','conventions':'strict thresholds; chronological graph; gaps > twice median cadence break edges','bounds':{'GP_draws':1024,'GP_scenarios':36,'Cassini_draws_per_model':512,'Cassini_passes':20,'epochs':1696},'non_claims':['No dimension estimate','No new force bound','No full hierarchical SBC','No raw telemetry refit','No universal superiority of a diagnostic']}
        manifest['randomness']={'used':True,'generator':'NumPy PCG64','seeds':[2026111901,2026111902,2026111903]}
        manifest['result']='Bounded topology cases executed; exploratory local Cassini discrepancy retained' if success else 'Execution failed; see logs'
    (R/f'records/execution_manifest_{mode}_v2.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if mode=='all':(R/'records/execution_manifest_v2.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if not success:raise SystemExit(1)
if __name__=='__main__':main()
