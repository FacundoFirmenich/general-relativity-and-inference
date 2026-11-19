#!/usr/bin/env python3
"""Execute the finite assertions in CONTRACT.md, record provenance, draw figures."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/pae_matplotlib')
from pathlib import Path
import csv, datetime, hashlib, json, math, platform, sys, time
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import ndtr, erfcx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dcae import from_logs, from_scales

ROOT=Path(__file__).resolve().parent
os.chdir(ROOT)
for d in ['results','figures','records']: Path(d).mkdir(exist_ok=True)
SEED=20261118
INPUTS=['dcae.py','verify.py','CONTRACT.md']
RESULTS=['results/checks.json','results/null_bias.csv','figures/diagnostics_ES.pdf','figures/diagnostics_EN.pdf','figures/repair_ES.pdf','figures/repair_EN.pdf']
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
before={p:sha(p) for p in INPUTS}
t0=time.monotonic();started=datetime.datetime.now(datetime.timezone.utc).isoformat()
checks=[]
def check(name,condition,details):
    checks.append({'assertion':name,'passed':bool(condition),'details':details})
    if not condition: raise AssertionError(name)

status='failed';exit_code=1
try:
    xs=np.linspace(0,3,101)
    log_lf=lambda t:1.3*t
    log_lo=lambda t:.7*t+.1*np.sin(t)
    errors=[abs(quad(lambda v:.6-.1*np.cos(v),0,float(t),epsabs=1e-12)[0]-(log_lf(t)-log_lo(t))) for t in xs]
    check('integral versus endpoints: 101 scales',max(errors)<1e-11,{'max_absolute_error':max(errors)})
    q=np.geomspace(1,100,1001)
    err=0.
    for beta in [0.,.5,1.,1.5]:
        direct=np.array([from_scales(x,x**beta,1,1)['index'] for x in q])
        err=max(err,float(np.max(np.abs(direct-np.maximum(0,1-q**(beta-1))))))
    check('power laws: 4004 values',err<1e-12,{'max_absolute_error':err})
    r=from_scales(10,2,1,.5);scaled=from_scales(10000,2000,1000,500)
    check('common unit invariance',abs(r['index']-scaled['index'])<1e-14,{'base':r,'rescaled':scaled})
    check('missing is not zero',from_scales(1,None,1,1)['status']=='NOT_ESTIMABLE',{})
    check('zero outside declared domain',from_scales(1,0,1,1)['status']=='INVALID_INPUT',{})
    check('saturation is flagged',from_logs(1000,0,0,0)['numerical_saturation'],from_logs(1000,0,0,0))
    i01=from_scales(2,1,1,1)['index'];i12=from_scales(1,1,2,1)['index'];i02=from_scales(1,1,1,1)['index']
    check('clipping breaks general composition',i02!=i01+i12-i01*i12,{'I01':i01,'I12':i12,'I02':i02})
    # Identical volume, distinct response locations: h=(1,0), selection=(1,0) or (0,1).
    h=np.array([1.,0.]);masks=[np.array([1.,0.]),np.array([0.,1.])]
    coverage=[float(m.mean()) for m in masks]
    fisher=[float(np.sum(m*h*h)) for m in masks]
    check('same coverage different scalar Fisher',coverage[0]==coverage[1] and fisher==[1.,0.],{'coverage':coverage,'fisher':fisher})
    A=np.array([[1.,0.]]);B=np.array([[0.,1.]])
    ranks=[int(np.linalg.matrix_rank(x)) for x in [A,np.vstack([A,A]),np.vstack([A,B])]]
    check('repetition versus complementary repair',ranks==[1,1,2],{'ranks':ranks,'null_direction':[0,1]})
    # Binary symmetric channel: error .1 followed by error .2 -> .26, exact rational in SI.
    P=np.array([[.9,.1],[.1,.9]]);K=np.array([[.8,.2],[.2,.8]]);PK=P@K
    check('finite garbling example',np.allclose(PK,[[.74,.26],[.26,.74]],rtol=0,atol=1e-15),{'P':P.tolist(),'K':K.tolist(),'PK':PK.tolist(),'equal_prior_Bayes_error':[.1,.26]})
    rng=np.random.Generator(np.random.PCG64(SEED));rows=[]
    for sigma in [.05,.1,.2,.5,1.]:
        analytic=.5-.5*erfcx(sigma/np.sqrt(2))
        numerical=quad(lambda z:-np.expm1(-sigma*z)*np.exp(-z*z/2)/np.sqrt(2*np.pi),0,12,epsabs=1e-12)[0]
        z=rng.normal(0,sigma,200000);iv=-np.expm1(-np.maximum(0,z))
        estimate=float(iv.mean());se=float(iv.std(ddof=1)/np.sqrt(len(iv)))
        rows.append({'sigma':sigma,'analytic_mean':float(analytic),'quadrature_mean':numerical,'MC_mean':estimate,'MC_se':se,'standardized_difference':float((estimate-analytic)/se),'N':len(iv)})
    check('Gaussian clipping mean against quadrature',max(abs(r['analytic_mean']-r['quadrature_mean']) for r in rows)<1e-11,rows)
    check('Monte Carlo consistency, five fixed sigma cases',max(abs(r['standardized_difference']) for r in rows)<6,{'criterion':'six simulation SE; finite stochastic check, not proof'})
    with open('results/null_bias.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    a=np.array([1,-1,-1,1]);C=np.array([[.04,.01],[.01,.09]])
    S=np.block([[C,C],[C,C]])
    check('shared reference cancels exactly',float(a@S@a)==0.,{'variance_same_reference':float(a@S@a),'wrong_independent_reference':float(a@scipy.linalg.block_diag(C,C)@a)})
    ss=np.array([1.,2.,4.,8.]);n=100;sig=2.
    sd=sig*ss**2/np.sqrt(n)
    J=n*ss**(-4)/sig**2
    check('linear Gaussian inverse-square exact CR bound',np.allclose(sd**2,1/J,rtol=1e-14),{'s':ss.tolist(),'sd':sd.tolist(),'n':n,'sigma':sig})
    kap=[2.,.5,3.];local=[.1,.2,.05];e=0.
    for k,d in zip(kap,local):e=k*e+d
    bound=sum(local[i]*np.prod(kap[i+1:]) for i in range(3))
    check('scalar chain attains Lipschitz bound',abs(e-bound)<1e-15,{'error':e,'bound':float(bound)})

    palette=['#263F69','#7B55A4','#18A7BF','#A64F91']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'axes.labelcolor':'#263F69','axes.titlecolor':'#263F69'})
    for lang in ['ES','EN']:
        es=lang=='ES';fig,axes=plt.subplots(2,2,figsize=(10,7.4),layout='constrained')
        ax=axes[0,0]
        for b,c in zip([0,.5,1,1.5],palette):ax.plot(q,np.maximum(0,1-q**(b-1)),color=c,label=rf'$\beta={b}$',ls='--' if b==1.5 else '-')
        ax.set(xscale='log',xlabel=r'$L_F/L_{F0}$',ylabel=r'$I_{\mathrm{DCAE}}$',title='A · Leyes condicionadas' if es else 'A · Conditional power laws');ax.legend(frameon=False,ncol=2)
        ax=axes[0,1];t=np.linspace(0,2*np.pi,300);d=.8*np.sin(t)
        ax.plot(t,d,color=palette[0],label=r'$\Delta$');ax.plot(t,-np.expm1(-np.maximum(0,d)),color=palette[3],label=r'$I$');ax.axhline(0,color='#AAAAAA',lw=.7)
        ax.set(xlabel=r'$t=\log(s/s_0)$',title='B · Crecimiento y recuperación' if es else 'B · Growth and recovery');ax.legend(frameon=False)
        ax=axes[1,0];ax.bar(['A','B'],fisher,color=[palette[1],palette[2]],width=.5)
        ax.set(ylim=(0,1.2),ylabel='Información de Fisher' if es else 'Fisher information',title='C · Misma cobertura = 1/2' if es else 'C · Same coverage = 1/2')
        ax.text(1,.04,'0',ha='center',color=palette[0])
        ax=axes[1,1];grid=np.linspace(0,1,200);ax.plot(grid,.5-.5*erfcx(grid/np.sqrt(2)),color=palette[1],label='Analítica' if es else 'Analytic')
        ax.errorbar([r['sigma'] for r in rows],[r['MC_mean'] for r in rows],yerr=[3*r['MC_se'] for r in rows],fmt='o',color=palette[2],ms=4,label='MC ± 3 SE')
        ax.set(xlabel=r'$\sigma_\Delta$',ylabel=r'$E[\widehat I\mid\Delta=0]$',title='D · Sesgo del recorte bajo la nula' if es else 'D · Clipping bias under the null');ax.legend(frameon=False)
        fig.suptitle('Controles sintéticos · no datos cosmológicos' if es else 'Synthetic controls · not cosmological data',fontsize=12,color=palette[0])
        fig.savefig(f'figures/diagnostics_{lang}.pdf');plt.close(fig)
        # Geometry of local likelihood surfaces; z is chi^2, not physical distance.
        fig=plt.figure(figsize=(10,4.1),layout='constrained')
        v=np.linspace(-2,2,70);u,w=np.meshgrid(v,v)
        for j,z in enumerate([u*u,u*u+w*w]):
            ax=fig.add_subplot(1,2,j+1,projection='3d');ax.plot_surface(u,w,z,color=palette[1 if j==0 else 2],alpha=.8,linewidth=0,rstride=2,cstride=2)
            ax.set(xlabel=r'$\theta_1$',ylabel=r'$\theta_2$',zlabel=r'$\chi^2$',title=('A · Una dirección invisible' if es else 'A · One invisible direction') if j==0 else ('A + B · Reparación complementaria' if es else 'A + B · Complementary repair'))
            ax.view_init(28,-55)
        fig.savefig(f'figures/repair_{lang}.pdf');plt.close(fig)
    status='completed';exit_code=0
finally:
    Path('results/checks.json').write_text(json.dumps({'checks':checks,'all_passed':status=='completed','scope':'synthetic finite assertions in CONTRACT.md'},indent=2))
    after={p:sha(p) for p in INPUTS}
    if before!=after:status='inputs-changed'
    outputs=[{'path':p,'sha256':sha(p),'kind':'result'} for p in RESULTS if Path(p).exists()]
    manifest={'schema_version':2,'claim_id':'PAE-DCAE-0.2-synthetic-contract','repository':{'commit':'unavailable: standalone snapshot identified by hashes','dirty':False},'command':[sys.executable,'verify.py'],'environment':{'software':[{'name':'Python','version':platform.python_version()},{'name':'NumPy','version':np.__version__},{'name':'SciPy','version':scipy.__version__},{'name':'Matplotlib','version':matplotlib.__version__}],'hardware':platform.platform()+'; '+platform.machine()},'mathematics':{'assertion_tested':'Finite statements enumerated in CONTRACT.md','coefficient_domain':'real models evaluated in IEEE float64; exact small examples separately proved in supplement','conventions':'positive scales; dimensionless s; conditional Gaussian cases','inputs':INPUTS,'bounds':{'quadrature_cases':101,'powerlaw_cases':4004,'null_MC_per_sigma':200000,'sigmas':[.05,.1,.2,.5,1.]},'non_claims':['No empirical validation of PAE','No exhaustive survey of models','Finite numerical checks are not universal proofs']},'randomness':{'used':True,'generator':'NumPy PCG64','seed':SEED},'run':{'started_at':started,'runtime_seconds':time.monotonic()-t0,'exit_status':exit_code,'status':status,'resource_limits':{'memory_bytes':None,'cpu_seconds':None,'max_cores':None,'max_threads':None}},'outputs':outputs,'checks':checks,'result':'Implementation and finite assertions verified in stated range' if status=='completed' else 'See failing check','residual_risks':['Self-review, no independent mathematical audit','Float64/library dependence','No real observational data'], 'input_artifacts':[{'path':p,'sha256':before[p],'sha256_after':after[p]} for p in INPUTS],'declared_results':RESULTS}
    Path('records/computation_manifest.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps({'status':status,'checks':len(checks),'seconds':manifest['run']['runtime_seconds'],'null_bias':rows if 'rows' in globals() else []},indent=2))
