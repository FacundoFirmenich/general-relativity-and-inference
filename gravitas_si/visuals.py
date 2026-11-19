"""Accessible blue/violet/white publication and offline interactive figures."""
from pathlib import Path
import json,runpy,shutil,os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from cycler import cycler
from .figure_language import translate_figure
R=Path(__file__).resolve().parents[1];O=R/'figures'
NAVY='#17365D';BLUE='#2455D6';CYAN='#18A7D1';VIOLET='#7653C7';PALE='#DDE8F7';WHITE='#FFFFFF'
CMAP=LinearSegmentedColormap.from_list('gravitas',[NAVY,BLUE,VIOLET,CYAN,WHITE])
def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':NAVY,'axes.labelcolor':NAVY,
        'axes.edgecolor':NAVY,'axes.titlecolor':NAVY,'xtick.color':NAVY,'ytick.color':NAVY,
        'grid.color':PALE,'figure.facecolor':WHITE,'axes.facecolor':WHITE,'savefig.facecolor':WHITE,
        'patch.edgecolor':NAVY,'hatch.color':NAVY,'axes.prop_cycle':cycler(color=[BLUE,VIOLET,CYAN,NAVY]),
        'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight','pdf.fonttype':42})
def save(fig,name):
    translate_figure(fig,os.environ.get('GRAVITAS_FIGLANG','es'))
    fig.savefig(O/(name+'.pdf'));fig.savefig(O/(name+'.svg'));fig.savefig(O/(name+'.png'),dpi=180);plt.close(fig)

def main(lang='es'):
    global O
    os.environ['GRAVITAS_FIGLANG']=lang
    O=R/'figures' if lang=='es' else R/'figures/en'
    O.mkdir(exist_ok=True,parents=True);style()
    # Base plots are regenerated from their data with an explicit, approved palette.
    script=(R/'analyses/gp/code/figuras.py').read_text()
    script=script.replace("blue='#17618b';orange='#c45f24';gray='#363d46'",f"blue='{BLUE}';orange='{VIOLET}';gray='{NAVY}'")
    script=script.replace("'k--'","'--',color=gray")
    exec(compile(script,str(R/'analyses/gp/code/figuras.py'),'exec'),{'__file__':str(R/'analyses/gp/code/figuras.py')})
    base=R/'analyses/gp/figures'
    if lang=='en':base=base/'en'
    for p in base.glob('*'):
        if p.is_file():shutil.copy2(p,O/p.name)
    rr=np.geomspace(.2,5,70);ll=np.geomspace(.08,30,80);x,y=np.meshgrid(rr,ll)
    alpha=.1
    response=lambda r,l:1+alpha*(1+r/l)*np.exp(-r/l)
    z=100*(response(x,y)/response(1,y)-1)
    np.savetxt(R/'results/yukawa_surface.csv',np.c_[np.log10(x.ravel()),np.log10(y.ravel()),z.ravel()],delimiter=',',header='log10_r_over_r0,log10_lambda_over_r0,relative_shape_percent',comments='')
    np.savetxt(R/'results/yukawa_surface_tikz.csv',np.c_[np.log10(x[::3,::3].ravel()),np.log10(y[::3,::3].ravel()),z[::3,::3].ravel()],delimiter=',',header='log10_r_over_r0,log10_lambda_over_r0,relative_shape_percent',comments='')
    fig=plt.figure(figsize=(11.3,4.8),layout='constrained');ax=fig.add_subplot(121,projection='3d')
    ax.plot_surface(np.log10(x),np.log10(y),z,cmap=CMAP,linewidth=0,antialiased=True,shade=False)
    ax.set(xlabel='log₁₀(r/r₀)',ylabel='log₁₀(λ/r₀)',zlabel='Cambio de forma (%)',title='Yukawa después de normalizar GM')
    ax.view_init(elev=26,azim=-130)
    for a in [ax.xaxis,ax.yaxis,ax.zaxis]:
        a.set_pane_color((1,1,1,1));a.line.set_color(NAVY);a.set_tick_params(colors=NAVY)
        a._axinfo['grid']['color']=PALE
    ax2=fig.add_subplot(122)
    for l,col in zip([.2,1,5,30],[VIOLET,BLUE,CYAN,NAVY]):
        ax2.semilogx(rr,100*(response(rr,l)/response(1,l)-1),label=f'λ/r₀ = {l:g}',color=col)
    ax2.axhline(0,color=NAVY,lw=.7);ax2.set(xlabel='r/r₀',ylabel='Cambio de forma (%)',title='Cortes que permiten comprobar la superficie')
    ax2.legend(frameon=False);ax2.grid(alpha=.6);save(fig,'yukawa_3d')
    import plotly.graph_objects as go
    colors=[[0,NAVY],[.25,BLUE],[.5,VIOLET],[.8,CYAN],[1,WHITE]]
    fig=go.Figure(go.Surface(x=np.log10(x),y=np.log10(y),z=z,colorscale=colors,
        colorbar=dict(title='Cambio (%)',tickfont=dict(color=NAVY)),
        hovertemplate='log₁₀(r/r₀)=%{x:.2f}<br>log₁₀(λ/r₀)=%{y:.2f}<br>Cambio=%{z:.3f}%<extra></extra>'))
    axspec=dict(backgroundcolor=WHITE,gridcolor=PALE,zerolinecolor=NAVY,color=NAVY,showbackground=True)
    fig.update_layout(template='none',title='GRAVITAS-SI | Respuesta Yukawa tras normalizar GM',
        font=dict(family='Arial',color=NAVY),paper_bgcolor=WHITE,plot_bgcolor=WHITE,
        scene=dict(xaxis=dict(title='log₁₀(r/r₀)',**axspec),yaxis=dict(title='log₁₀(λ/r₀)',**axspec),
                   zaxis=dict(title='Cambio de forma (%)',**axspec),bgcolor=WHITE,camera=dict(eye=dict(x=-1.6,y=-1.4,z=1.1))),
        modebar=dict(bgcolor=WHITE,color=BLUE,activecolor=VIOLET),margin=dict(l=10,r=10,t=75,b=10))
    html=fig.to_html(full_html=False,include_plotlyjs=True,config={'displaylogo':False,'responsive':True})
    (O/'yukawa_interactivo.html').write_text('<!doctype html><html lang="es"><meta charset="utf-8"><title>GRAVITAS-SI</title><style>body{background:'+WHITE+';color:'+NAVY+';font:17px Arial;margin:2rem}a{color:'+BLUE+'}h1{color:'+NAVY+'}</style><h1>Geometría observable de una corrección Yukawa</h1><p>α = 0,1. Se divide por la respuesta en r₀ para eliminar una normalización de GM. Es un cálculo analítico ilustrativo: no representa una detección ni un ajuste de efemérides.</p>'+html+'<p>Rotar y ampliar con el ratón. Todas las bibliotecas están incluidas; funciona sin conexión. Fuente numérica: <code>results/yukawa_surface.csv</code>.</p><p><a href="https://orcid.org/0009-0002-6578-3811">Facundo Firmenich · ORCID 0009-0002-6578-3811</a></p></html>',encoding='utf-8')
    j=json.loads((R/'results/hierarchical.json').read_text());a=np.load(R/'results/hierarchical_arrays.npz')
    fig,axs=plt.subplots(1,3,figsize=(11.2,3.5),layout='constrained')
    for k,col in enumerate([BLUE,VIOLET]):
        if k==0:
            from scipy.special import ndtr
            xx=np.linspace(-.006,.003,450);density=np.zeros_like(xx)
            from scipy.stats import norm
            select=a['weights']>1e-8
            for mu,var,w in zip(a['mean'][select],a['variance'][select],a['weights'][select]):density+=w*norm.pdf(xx,mu*.001,np.sqrt(var)*.001)
            axs[0].plot(xx*1000,density/1000,color=BLUE);axs[0].fill_between(xx*1000,0,density/1000,color=CYAN,alpha=.15)
    axs[0].axvline(0,color=NAVY,ls='--');axs[0].set(xlabel='A × 10³ (amplitud óptica)',ylabel='Densidad',title='Posterior condicional')
    for ts,col in zip([3e-14,1e-13,3e-13],[BLUE,VIOLET,CYAN]):
        rows=[r for r in j['sensitivity'] if np.isclose(r['tau_prior_scale']/ts,1)]
        axs[1].semilogx([r['A_prior_sd'] for r in rows],[r['log_B10'] for r in rows],'o-',color=col,label=f'escala τ={ts:.0e}')
    axs[1].axhline(0,color=NAVY,ls='--');axs[1].set(xlabel='DE previa de A',ylabel='log B₁₀',title='La evidencia depende del prior');axs[1].legend(fontsize=8,frameon=False)
    axs[2].hist(a['ppc_rep'][:,0],bins=22,density=True,color=CYAN,alpha=.5,label='Réplicas')
    axs[2].hist(a['ppc_obs'][:,0],bins=22,density=True,histtype='step',lw=2,color=VIOLET,label='Residuos observados')
    axs[2].set(xlabel='Producto de residuos adyacentes',ylabel='Densidad',title='Fallo predictivo temporal');axs[2].legend(fontsize=8,frameon=False)
    for ax in axs:ax.grid(alpha=.5)
    save(fig,'jerarquia_diagnosticos')
    d=np.load(R/'results/desi_operator_arrays.npz')
    fig,axs=plt.subplots(1,2,figsize=(10,3.2),layout='constrained')
    axs[0].semilogy(np.arange(1,len(d['singular'])+1),d['singular'],'o',ms=3,color=BLUE)
    axs[0].set(xlabel='Índice del valor singular',ylabel='Valor singular (unidades del producto)',title='DESI: 72 direcciones de respuesta')
    axs[1].plot(d['null'],color=VIOLET,lw=.7);axs[1].set(xlabel='Coordenada del espectro de entrada',ylabel='Vector de norma unidad',title='Una de 975 direcciones no observables')
    for ax in axs:ax.grid(alpha=.5)
    save(fig,'desi_operador')
    import csv
    files=list((R/'analyses/casimir/results').glob('*.csv'))
    tablefile=next(f for f in files if 'pressure' in f.name)
    rows=list(csv.DictReader(tablefile.open()))
    fig,axs=plt.subplots(1,2,figsize=(10,3.2),layout='constrained')
    gap=np.array([float(r['gap_nm']) for r in rows])
    raw=np.array([float(r['table_Drude_Pa']) for r in rows]);correct=np.array([float(r['independent_T0_Drude_Pa']) for r in rows])
    axs[0].loglog(gap,-raw,color=BLUE,label='Tabla Drude');axs[0].loglog(gap,-correct,color=CYAN,ls='--',label='Integral T=0')
    axs[0].loglog(gap,-np.pi*raw,color=VIOLET,label='Tabla × π');axs[0].set(xlabel='Separación (nm)',ylabel='Magnitud de presión (Pa)',title='Normalización física')
    axs[1].plot(gap,correct/raw,color=BLUE);axs[1].axhline(np.pi,color=VIOLET,ls='--');axs[1].set(xlabel='Separación (nm)',ylabel='Cociente de presiones',title='Una escala errónea conserva R²')
    axs[0].legend(frameon=False)
    for ax in axs:ax.grid(alpha=.5)
    save(fig,'casimir_normalizacion')
    # Vector TikZ/PGFPlots source reads the very same calculated table.
    (O/'yukawa_pgfplots.tex').write_text(r'''\begin{tikzpicture}
\begin{axis}[width=.85\linewidth,height=8cm,view={-130}{26},xlabel={$\log_{10}(r/r_0)$},ylabel={$\log_{10}(\lambda/r_0)$},zlabel={Cambio (\%)},axis line style={navy},tick style={navy},colormap={gv}{color(0cm)=(navy);color(1cm)=(blue);color(2cm)=(violet);color(3cm)=(cyan)},mesh/ordering=x varies]
\addplot3[surf,shader=interp,mesh/cols=24] table[x=log10_r_over_r0,y=log10_lambda_over_r0,z=relative_shape_percent,col sep=comma]{results/yukawa_surface_tikz.csv};
\end{axis}\end{tikzpicture}
''')
    (O/'yukawa_tikz_document.tex').write_text(r'''\documentclass[border=6pt]{standalone}
\usepackage{pgfplots}\pgfplotsset{compat=1.18}
\definecolor{navy}{HTML}{17365D}\definecolor{blue}{HTML}{2455D6}
\definecolor{cyan}{HTML}{18A7D1}\definecolor{violet}{HTML}{7653C7}
\begin{document}\color{navy}\input{figures/yukawa_pgfplots.tex}\end{document}
''')
    if lang=='en':
        p=O/'yukawa_pgfplots.tex';p.write_text(p.read_text().replace('Cambio','Change'))
        p=O/'yukawa_tikz_document.tex';p.write_text(p.read_text().replace('figures/yukawa_pgfplots.tex','figures/en/yukawa_pgfplots.tex'))
    print('Seven publication figures, offline 3D HTML, and PGFPlots source regenerated.')
if __name__=='__main__':
    main('es')
    main('en')
