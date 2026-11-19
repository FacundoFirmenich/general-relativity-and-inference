"""Execute the bounded contract in CONTRACT.md; no new physical bounds."""
from pathlib import Path
import sys,json,csv
import numpy as np
from scipy.linalg import cho_factor,cho_solve,eigh
from scipy.special import logsumexp
R=Path(__file__).resolve().parents[2];sys.path.insert(0,str(R))
from research.fractalgo.core import posterior,path_summary,cell_weights,excursion,path_edges
OUT=R/'research/fractalgo/results';OUT.mkdir(exist_ok=True)
SEED=2026111902

def summary(x):
    x=np.asarray(x,dtype=float)
    return dict(mean=float(x.mean()),mcse=float(x.std(ddof=1)/np.sqrt(len(x))),q025=float(np.quantile(x,.025)),q975=float(np.quantile(x,.975)))

def gp_case():
    rng=np.random.default_rng(SEED);nsim=1024;x=np.linspace(0,3,201);nodes=np.array([0.,.6,1.5,3.])
    kernel=lambda a,b:9*np.exp(-.5*((a[:,None]-b[None,:])/.5)**2)
    cross=kernel(x,nodes);C=kernel(x,x)-cross@cho_solve(cho_factor(kernel(nodes,nodes)),cross.T);C=(C+C.T)/2
    ev,U=eigh(C);assert ev.min()>-1e-9
    keep=ev>1e-10;L=U[:,keep]*np.sqrt(ev[keep]);assert np.max(abs(L@L.T-C))<1e-9
    joint=rng.normal(size=(nsim,keep.sum()))@L.T
    independent=rng.normal(size=(nsim,len(x)))*np.sqrt(np.maximum(C.diagonal(),0))
    src=json.loads((R/'analyses/gp/results/gp_emulator.json').read_text());mean=x**4-2*src['z_hat']*x*x
    rows=[];draw_rows=[]
    for count in [51,101,201]:
        ids=np.linspace(0,200,count,dtype=int);xx=x[ids];weights=cell_weights(xx)
        for prior_name in ['uniform_x','exp_minus_x']:
            prior=weights*(1 if prior_name=='uniform_x' else np.exp(-xx));prior/=prior.sum()
            for law,delta in [('joint_GP',joint),('independent_marginals',independent)]:
                q=mean[ids]+delta[:,ids];l=-.5*(q-q.min(axis=1)[:,None]);lr=l-logsumexp(l+np.log(prior),axis=1)[:,None];p=np.exp(lr+np.log(prior))
                assert np.allclose(p.sum(axis=1),1)
                for threshold in [.5,1.,2.]:
                    counts,longest,largest=path_summary(lr,np.log(threshold),mass=p)
                    selected=(p*(lr>np.log(threshold))).sum(axis=1)
                    row=dict(nodes=count,prior=prior_name,law=law,threshold=threshold,components=summary(counts),largest_component_mass=summary(largest),selected_mass=summary(selected));rows.append(row)
                    draw_rows.extend((count,prior_name,law,threshold,k,int(c),float(a),float(b)) for k,(c,a,b) in enumerate(zip(counts,largest,selected)))
    with (OUT/'gp_draws.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['nodes','prior','law','threshold','draw','components','largest_component_mass','selected_mass']);w.writerows(draw_rows)
    np.savez_compressed(OUT/'gp_inputs_and_draws.npz',x=x,mean=mean,covariance=C,joint_delta=joint,independent_delta=independent)
    out={'scope':'Conditional-cost-profile Gaussian uncertainty experiment; not recovered INPOP training data',
         'profile_source':'analyses/gp/results/gp_emulator.json','profile_z_hat':src['z_hat'],'seed':SEED,'draws':nsim,
         'covariance':'Squared-exponential, amplitude 3, length .5; conditioned at 0,.6,1.5,3',
         'mesh':'nested 51/101/201; trapezoidal cell prior masses, not exact-segment continuous integration',
         'rows':rows,'nonclaims':['No new force or graviton bound','No dimension estimate','No continuum convergence proof','No calibrated model-selection threshold']}
    (OUT/'gp.json').write_text(json.dumps(out,indent=2)+'\n')
    print('GP completed:',len(rows),'scenarios',flush=True)

def cassini_case():
    from research.revision_v17.cassini import Model,integrate
    lengths=[300.,1800.,7200.,86400.];models=[];refs=[];checks=[]
    for ell in lengths:
        m=Model(ell);coarse=integrate(m);ref=integrate(m,True);err=abs(coarse['logz']-ref['logz']);assert err<.002
        models.append(m);refs.append(ref);checks.append(err)
    probs=np.exp(np.array([r['logz'] for r in refs])-logsumexp([r['logz'] for r in refs]))
    m=models[0];graphs=[];node_rows=[]
    for ids,_ in m.blocks:
        ids=ids[np.argsort(m.epoch[ids],kind='stable')];dt=np.diff(m.epoch[ids]);assert (dt>0).all()
        gap=2*np.median(dt);connected=dt<=gap
        graphs.append((ids,connected));node_rows.extend(dict(original_row=int(i),epoch=float(m.epoch[i]),pass_id=int(m.pid[i])) for i in ids)
    assert len(node_rows)==1696 and len({r['original_row'] for r in node_rows})==1696
    graph_record={'nodes':node_rows,'edges':[[int(ids[i]),int(ids[i+1])] for ids,c in graphs for i in np.flatnonzero(c)],'gap_rule':'within-pass dt <= 2*median positive cadence; no cross-pass edges'}
    (OUT/'cassini_graph.json').write_text(json.dumps(graph_record,indent=2)+'\n')
    rng=np.random.default_rng(SEED+1);draws=512;thresholds=[0.,1.,2.];rows=[];arrays={}
    for name,prob in [('reference_7200',np.array([0,0,1,0.])),('mixture',probs)]:
        obs=np.empty((draws,len(graphs),len(thresholds),3));rep=np.empty_like(obs)
        choices=rng.choice(4,draws,p=prob)
        for k,j in enumerate(choices):
            m=models[j];ref=refs[j];x=ref['x'][rng.choice(len(ref['x']),p=ref['weights'])];s=np.exp(x[0]);tau=np.exp(x[1])
            latent=m.Q@(100*rng.normal(size=m.p))+tau*np.sqrt(m.ev)*rng.normal(size=m.n)+3*m.t*rng.normal()
            noise=s*rng.normal(size=m.n)
            posterior_latent=m.y-noise-s*s*m.inverse(x,m.y-latent-noise)
            observed=m.unrotate(m.y-posterior_latent)/s
            replicated=m.unrotate(rng.normal(size=m.n))
            for b,(ids,c) in enumerate(graphs):
                for h,t in enumerate(thresholds):
                    for target,values in [(obs,observed),(rep,replicated)]:
                        # Union of positive and negative excursions, with no edge across sign changes.
                        cp,lp,_=path_summary(values[ids],t,c);cn,ln,_=path_summary(-values[ids],t,c)
                        target[k,b,h]=[int(cp[0]+cn[0]),int(max(lp[0],ln[0])),float(np.mean(abs(values[ids])>t))]
        arrays[name+'_observed']=obs;arrays[name+'_replicated']=rep
        for b,(ids,c) in enumerate(graphs):
            for h,t in enumerate(thresholds):
                a=obs[:,b,h];z=rep[:,b,h]
                lower=z[:,0]<=a[:,0];upper=z[:,1]>=a[:,1]
                rows.append(dict(model=name,pass_id=int(m.pid[ids[0]]),epochs=len(ids),threshold=t,broken_edges=int((~c).sum()),
                    components_observed_median=float(np.median(a[:,0])),components_predictive_025=float(np.quantile(z[:,0],.025)),components_predictive_975=float(np.quantile(z[:,0],.975)),
                    longest_observed_median=float(np.median(a[:,1])),longest_predictive_025=float(np.quantile(z[:,1],.025)),longest_predictive_975=float(np.quantile(z[:,1],.975)),
                    coverage_observed_median=float(np.median(a[:,2])),coverage_predictive_median=float(np.median(z[:,2])),
                    components_lower_tail_count=int(lower.sum()),longest_upper_tail_count=int(upper.sum()),draws=draws,
                    components_lower_tail=float(lower.mean()),longest_upper_tail=float(upper.mean())))
        print('Cassini topology completed:',name,flush=True)
    np.savez_compressed(OUT/'cassini_predictive_draws.npz',thresholds=thresholds,**arrays)
    out={'scope':'Exploratory posterior predictive residual topology conditional on reconstructed DSS-25 inputs; no raw telemetry refit',
         'seed':SEED+1,'draws_per_model':draws,'kernel_lengths_seconds':lengths,'kernel_probabilities':probs.tolist(),
         'quadrature_refinement_errors':checks,'rows':rows,
         'limits':['Posterior predictive fractions are not uniformly calibrated p-values','Discrete ties are included in tails','No multiple-testing exclusion rule','Thresholds/graph selected after prior residual diagnosis; exploratory','Same fixed nuisance basis and inherited empirical weights','No full hierarchical SBC or independent instrument calibration']}
    (OUT/'cassini.json').write_text(json.dumps(out,indent=2)+'\n')

def main():gp_case();cassini_case()
if __name__=='__main__':main()
