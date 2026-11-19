"""Conditional block prediction and noise-scale mixture, reconstructed DSS-25 only.
No physical calibration, new telemetry, station comparison or ephemeris refit.
"""
from pathlib import Path
import sys,json,time,datetime
import numpy as np
from scipy.linalg import eigh
from scipy.special import logsumexp
R=Path(__file__).resolve().parents[2];sys.path.insert(0,str(R))
from gravitas_si.inference import ConditionalModel,logprior,BASE,AUNIT
from analyses.gp.code.gp_covarianza_cassini import intkernel

class Model(ConditionalModel):
    def __init__(self,ell):
        src=R/'analyses/gp/sources/cassini';z=np.load(src/'dss25_diseno_erfa_constant.npz');t=np.load(src/'contraste_robusto_pesos.npz')
        self.pid=z['pass_id'];self.w=z['weights'];self.epoch=z['epoch'];xx=z['design']*self.w[:,None]
        Q,s,_=np.linalg.svd(xx/np.linalg.norm(xx,axis=0),full_matrices=False);assert s[-1]>1e-8
        yy=self.w*(z['residual']+t['constant__delta_SPI']);tt=self.w*t['constant__template_per_A']*AUNIT
        self.blocks=[];rows=[];ev=[];self.rotated_groups=[];start=0
        for pno in np.unique(self.pid):
            ids=np.flatnonzero(self.pid==pno);dt=self.epoch[ids,None]-self.epoch[None,ids]
            K=intkernel(dt,ell)*(BASE*self.w[ids,None])*(BASE*self.w[None,ids]);e,U=eigh(K)
            assert e[0]>-1e-7*max(1,e[-1]);ev.extend(np.maximum(e,0));rows.append(U.T@np.c_[yy[ids],tt[ids],Q[ids]])
            self.blocks.append((ids,U));self.rotated_groups.append(np.arange(start,start+len(ids)));start+=len(ids)
        self.ev=np.array(ev);Z=np.vstack(rows);self.y,self.t,self.Q=Z[:,0],Z[:,1],Z[:,2:];self.n,self.p=self.Q.shape;self.observed=yy

def integrate(m,expanded=False):
    u=np.linspace(-.5,.5,151) if expanded else np.linspace(-.25,.25,101)
    v=np.linspace(-17,3,241) if expanded else np.linspace(-13,2,181)
    x=np.array(np.meshgrid(u,v,indexing='ij')).reshape(2,-1).T
    F=np.c_[100*m.Q,3*m.t];Z=np.c_[m.y,F];Gs=[];lds=[]
    for ids in m.rotated_groups:
        d=np.exp(2*x[:,0,None])+np.exp(2*x[:,1,None])*m.ev[ids]
        G=np.einsum('gi,ij,ik->gjk',1/d,Z[ids],Z[ids],optimize=True)
        Gs.append(G);lds.append(np.log(d).sum(axis=1))
    G=sum(Gs);ld=sum(lds)
    def ll(G,ld,n):
        H=G[:,1:,1:]+np.eye(F.shape[1]);b=G[:,1:,0]
        sol=np.linalg.solve(H,b[...,None])[...,0]
        return -.5*(n*np.log(2*np.pi)+ld+np.linalg.slogdet(H)[1]+G[:,0,0]-np.einsum('gi,gi->g',b,sol))
    wt=np.ones((len(u),len(v)));wt[[0,-1]]*=.5;wt[:,[0,-1]]*=.5
    lp=np.array([logprior(xx) for xx in x])+np.log(wt.ravel()*(u[1]-u[0])*(v[1]-v[0]))
    full=ll(G,ld,m.n);lz=logsumexp(full+lp);w=np.exp(full+lp-lz)
    train=[]
    for ids,g,d in zip(m.rotated_groups,Gs,lds):train.append(float(logsumexp(ll(G-g,ld-d,m.n-len(ids))+lp)))
    ix=np.argmax(w);assert abs(full[ix]-m.loglike(x[ix]))<1e-6
    return dict(logz=float(lz),training_logz=train,x=x,weights=w)

def ppc(models,refs,draws=1024):
    rng=np.random.default_rng(2026111807);probs=np.exp(np.array([r['logz'] for r in refs])-logsumexp([r['logz'] for r in refs]));obs=[];rep=[]
    choices=rng.choice(len(models),draws,p=probs)
    def stat(res,m):return np.array([np.mean(res[ids][:-1]*res[ids][1:]) for ids,U in m.blocks])
    for j in choices:
        m=models[j];ref=refs[j];x=ref['x'][rng.choice(len(ref['x']),p=ref['weights'])];s=np.exp(x[0]);tau=np.exp(x[1])
        lat=m.Q@(100*rng.normal(size=m.p))+tau*np.sqrt(m.ev)*rng.normal(size=m.n)+3*m.t*rng.normal()
        noise=s*rng.normal(size=m.n);latpost=m.y-noise-s*s*m.inverse(x,m.y-lat-noise)
        obs.append(stat(m.unrotate(m.y-latpost)/s,m));rep.append(stat(m.unrotate(rng.normal(size=m.n)),m))
    obs=np.array(obs);rep=np.array(rep);counts=np.array([len(ids)-1 for ids,U in models[0].blocks])
    return dict(draws=draws,pass_observed_median=np.median(obs,axis=0).tolist(),pass_replicated_025=np.quantile(rep,.025,axis=0).tolist(),pass_replicated_975=np.quantile(rep,.975,axis=0).tolist(),pass_tail_fractions=np.mean(rep>=obs,axis=0).tolist(),aggregate_tail_fraction=float(np.mean(rep@counts>=obs@counts)),kernel_posterior_probabilities=probs.tolist())

def main():
    start=time.time();ells=[300.,1800.,7200.,86400.];models=[];refs=[];checks=[]
    for ell in ells:
        m=Model(ell);r=integrate(m);wide=integrate(m,True);err=abs(r['logz']-wide['logz']);assert err<.002,(ell,err)
        models.append(m);refs.append(wide);checks.append(err);print('Integrated ell',ell,'logZ',wide['logz'],'refinement',err,flush=True)
    lp=np.array([r['logz'] for r in refs]);tr=np.array([r['training_logz'] for r in refs]);mixlogz=logsumexp(lp)-np.log(len(lp));trainmix=logsumexp(tr,axis=0)-np.log(len(lp))
    ppc_ref=ppc([models[2]],[refs[2]]);ppc_mix=ppc(models,refs)
    m=models[0];rows=[]
    for i,(ids,U) in enumerate(m.blocks):
        utc=lambda e:datetime.datetime.fromtimestamp(float(e),datetime.timezone.utc).isoformat()
        rows.append(dict(pass_id=int(m.pid[ids[0]]),n=len(ids),start_utc=utc(m.epoch[ids].min()),end_utc=utc(m.epoch[ids].max()),reference_log_predictive=lp[2]-tr[2,i],mixture_log_predictive=mixlogz-trainmix[i],difference=(mixlogz-trainmix[i])-(lp[2]-tr[2,i]),reference_ppc_tail=ppc_ref['pass_tail_fractions'][i],mixture_ppc_tail=ppc_mix['pass_tail_fractions'][i]))
    out=R/'research/revision_v17/results';out.mkdir(exist_ok=True)
    result=dict(scope='Conditional on reconstructed DSS-25 hybrid residuals, fixed weights, fixed optical range and full-design nuisance basis; pass deletion is not independent external validation',kernel_lengths_seconds=ells,kernel_prior=[.25]*4,logz_by_kernel=lp.tolist(),quadrature_refinement_errors=checks,reference_ppc=ppc_ref,mixture_ppc=ppc_mix,pass_predictions=rows,sum_log_score_difference=float(sum(r['difference'] for r in rows)),limits=['No inter-station comparison: this design is DSS-25','Weights and reconstruction were obtained before this diagnostic; no raw-data cross validation','Noise-length alternatives chosen after seeing the original failure; exploratory, not preregistered confirmation','Same data-dependent design prior is held fixed in folds','No physical alpha/lambda posterior or all-parameter ephemeris hierarchy'],runtime_seconds=time.time()-start)
    (out/'cassini.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['reference_ppc','mixture_ppc','pass_predictions']},indent=2),flush=True)
    figures(result)

def figures(result):
    rows=result['pass_predictions'];ppc_ref=result['reference_ppc'];ppc_mix=result['mixture_ppc']
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for lang in ['en','es']:
        en=lang=='en';fig,axs=plt.subplots(1,2,figsize=(10,3.5),layout='constrained');x=np.arange(len(rows))
        for p,c,label,shift in [(ppc_ref,'#7653C7','7200 s',-.12),(ppc_mix,'#18A7D1','Mixture' if en else 'Mezcla',.12)]:
            axs[0].plot(x+shift,p['pass_observed_median'],'o',color=c,ms=3,label=label)
        axs[0].fill_between(x,ppc_ref['pass_replicated_025'],ppc_ref['pass_replicated_975'],color='#7653C7',alpha=.12,label='95% predictive reference' if en else 'Referencia predictiva 95%')
        axs[0].set(ylabel='Adjacent residual product' if en else 'Producto residual adyacente',xlabel='Pass index (chronological)' if en else 'Índice de pasada (cronológico)');axs[0].legend(frameon=False,fontsize=7)
        vals=[r['difference'] for r in rows];axs[1].bar(x,vals,color='#9A5498');axs[1].axhline(0,color='#17365D',linewidth=.7);axs[1].set(ylabel='Mixture − reference log score' if en else 'Log-score mezcla − referencia',xlabel='Held-out pass index' if en else 'Índice de pasada excluida')
        for ax in axs:ax.spines[['top','right']].set_visible(False)
        fig.savefig(R/f'figures/cassini_blocks_{lang}.pdf');plt.close(fig)
if __name__=='__main__':main()
