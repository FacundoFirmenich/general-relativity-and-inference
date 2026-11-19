"""Bilingual TikZ views of the same epistemic model used by the offline HTML."""
from pathlib import Path
import json
R=Path(__file__).resolve().parents[1]
def main():
 g=json.loads((R/'ui/epistemic/graph.json').read_text());nodes={n['id']:n for n in g['nodes']}
 for lang in ['es','en']:
  def T(es,en):return es if lang=='es' else en
  def b(v):return v[lang].replace('&',r'\&')
  out=R/'figures'/lang;out.mkdir(exist_ok=True)
  lines=[r'\begin{tikzpicture}[x=1cm,y=1cm,>=Latex,font=\small,text=navy]',r'\path[use as bounding box] (-.1,.6) rectangle (16.4,-14.0);',r'\node[anchor=center,font=\sffamily\bfseries\small,text=violet] at (8.1,.3) {'+T('MAPA EPISTÉMICO · DEL RESULTADO A SU ALCANCE','EPISTEMIC MAP · FROM RESULTS TO THEIR SCOPE')+'};',r'\node[draw=navy!35,fill=navy!3,rounded corners=3pt,text width=15.8cm,align=center,inner sep=8pt] at (8.1,-.7) {\textbf{'+T('Premisa física y operador de observación','Physical premise and observation operator')+r'}\\[3pt]'+T('Plantilla efectiva explícita; calibración de masa; excentricidad; predicciones propias de cada teoría.','Explicit effective template; mass calibration; eccentricity; theory-specific predictions.')+'};']
  for j,t in enumerate([T('CRITERIO','CRITERION'),T('EVIDENCIA Y CONTROL','EVIDENCE AND CHECK'),T('LÍMITE DE LA AFIRMACIÓN','CLAIM BOUNDARY')]):
   lines.append(r'\node[anchor=west,font=\sffamily\bfseries\footnotesize,text=violet] at ('+str(j*5.55)+',-1.7) {'+t+'};')
  rows=[('gp',r'$Db=Db^{\prime}$; $DKD^T=DK^{\prime}D^T$',T('Bandas marginales iguales no fijan la ley posterior.','Equal marginal bands do not fix the posterior law.'),r'$D_{\rm KL}=1.402585$ nats',T('Control sintético de covarianza; soportes singulares explícitos.','Synthetic covariance control; singular supports explicit.'),T('No es un factor de corrección empírico de INPOP.','Not an empirical INPOP correction factor.')),
  ('cassini',r'$y=A s+X\beta+\epsilon$',T('1696 épocas reconstruidas; cuatro covarianzas temporales.','1696 reconstructed epochs; four temporal covariances.'),T('Pasada 12: $0/1024$ excedencias','Pass 12: $0/1024$ exceedances'),T('La mezcla mejora el control agregado; persiste el fallo local.','The mixture improves the aggregate check; local failure persists.'),T(r'Sin reducción primaria independiente ni posterior físico $(\alpha,\lambda)$.',r'No independent primary reduction or physical $(\alpha,\lambda)$ posterior.')),
  ('information',r'$K_{\pi,L}=1-R_{\pi,L}(E)/R_{\pi,L}(\varnothing)$',T('Retención relativa a un objeto, prior y pérdida declarados.','Retention relative to a declared target, prior and loss.'),r'$J_{\rm mix}=0.345509<0.555556$',T('Control de mezcla; Fisher, Rényi y riesgo distinguen pérdidas.','Mixture control; Fisher, Rényi and risk distinguish losses.'),T('No es una probabilidad de verdad ni una barrera universal.','Not a truth probability or a universal barrier.')),
  ('dcae',r'$I_+=1-e^{-\max(0,\Delta)}$',T(r'Conservar $\Delta$ firmado junto al índice recortado.',r'Retain signed $\Delta$ alongside the clipped index.'),r'$\Delta_{\rm soporte}=-0.0092143$' if lang=='es' else r'$\Delta_{\rm support}=-0.0092143$',T('COSMOS: el orden cambia con tolerancia y observable.','COSMOS: ordering changes with tolerance and observable.'),T('Cobertura no identifica. DCAE físico: no estimable con estas tablas.','Coverage does not identify. Physical DCAE: not estimable from these tables.'))]
  for row,(r,f,desc,metric,check,limit) in enumerate(rows):
   y=-3-row*2.55;color=['violet','plum','cyan','navy'][row]
   titles=[b(g['routes'][r]),b(nodes[r+'_check']['title']),T('Alcance condicionado','Conditional scope')]
   bodies=[f+r'\\[5pt]'+desc,metric+r'\\[5pt]'+check,limit]
   for j in range(3):
    lines.append(r'\node[draw='+color+r'!55,fill='+color+r'!3,rounded corners=3pt,text width=4.65cm,minimum height=1.95cm,align=left,inner sep=7pt'+(',dashed' if j==2 else '')+'] ('+r+str(j)+') at ('+str(2.45+j*5.55)+','+str(y)+r') {\textbf{'+titles[j]+r'}\\[6pt]'+bodies[j]+'};')
   for j in range(2):lines.append(r'\draw[->,'+color+r'!60,thick] ('+r+str(j)+'.east)--('+r+str(j+1)+'.west);')
  lines+=[r'\node[draw=violet!40,fill=violet!5,rounded corners=2pt,text width=15.8cm,align=center,inner sep=7pt] at (8.1,-12.8) {\textbf{'+T('Identificabilidad + ley conjunta + adecuación predictiva','Identifiability + joint law + predictive adequacy')+r'}\\[4pt]'+T('Las comprobaciones numéricas no sustituyen la validación instrumental independiente.','Numerical checks do not replace independent instrumental validation.')+'};',r'\node[anchor=west,font=\tiny,text=navy!65] at (0,-13.8) {'+T('Mismo modelo de contenidos que el grafo interactivo. SHA-256: ','Same content model as the interactive graph. SHA-256: ')+g['content_sha256'][:16]+'.};',r'\end{tikzpicture}']
  (out/'epistemic_graph.tex').write_text('\n'.join(lines)+'\n')
  lines=[r'\begin{tikzpicture}[x=1cm,y=1cm,>=Latex,font=\small,text=navy]',r'\path[use as bounding box] (-.1,.6) rectangle (16.4,-7.9);',r'\node[anchor=center,font=\sffamily\bfseries\small,text=violet] at (8.1,.3) {'+T('RELATIVIDAD GENERAL · QUÉ AUTORIZA UNA INFERENCIA','GENERAL RELATIVITY · WHAT WARRANTS AN INFERENCE')+'};',r'\node[anchor=center,align=center,font=\large\bfseries,text width=16cm] at (8.1,-.6) {'+T('Mismos resúmenes no garantizan las mismas conclusiones.','The same summaries need not warrant the same conclusions.')+'};']
  chain=[T('Modelo físico','Physical model'),T('Operador de medida','Measurement operator'),T('Ley conjunta','Joint law'),T('Contraste predictivo','Predictive check')]
  for j,title in enumerate(chain):
   lines.append(r'\node[draw=navy!30,fill=navy!3,rounded corners=2pt,text width=3.5cm,minimum height=.8cm,align=center,inner sep=5pt] (c'+str(j)+') at ('+str(1.9+j*4.15)+r',-1.7) {\textbf{'+title+'}};')
   if j:lines.append(r'\draw[->,navy!60,thick] (c'+str(j-1)+'.east)--(c'+str(j)+'.west);')
  panels=[('violet',T('01 · EMULACIÓN GP','01 · GP EMULATION'),r'$Db,\ DKD^T$',T('Preservar contrastes,\nno sólo bandas marginales.','Preserve contrasts,\nnot only marginal bands.')),
   ('cyan',T('02 · INFORMACIÓN','02 · INFORMATION'),r'$K_{\pi,L}=1-R(E)/R(\varnothing)$',T('Medir la pérdida para\nun objeto y una decisión.','Measure loss for\na target and a decision.')),
   ('plum',T('03 · CASSINI','03 · CASSINI'),T('Pasada 12: $0/1024$','Pass 12: $0/1024$'),T('Conservar el fallo local\naunque mejore el agregado.','Retain local failure\neven if the aggregate improves.'))]
  for j,(color,title,formula,body) in enumerate(panels):
   x=2.5+j*5.6
   lines.append(r'\node[draw='+color+r'!55,fill='+color+r'!4,rounded corners=4pt,text width=4.55cm,minimum height=2.65cm,align=center,inner sep=8pt] at ('+str(x)+r',-4.05) {\textcolor{'+color+r'}{\sffamily\bfseries '+title+r'}\\[12pt]'+formula+r'\\[12pt]'+body.replace('\n',r'\\')+'};')
  lines += [r'\node[draw=violet!35,fill=white,rounded corners=3pt,align=center,text width=15.7cm,inner sep=10pt] at (8.1,-6.45) {\textbf{'+T(r'Convergencia numérica $\neq$ validación física',r'Numerical convergence $\neq$ physical validation')+r'}\\[6pt]'+T('La contribución: criterios verificables, jerarquía explícita y evidencia trazable.','The contribution: testable criteria, an explicit hierarchy and traceable evidence.')+'};',r'\node[anchor=west,font=\footnotesize,text=navy!70] at (0,-7.5) {'+T('Cassini: reconstrucción condicionada. GP e información: controles definidos. DCAE: descriptor, no ley universal.','Cassini: conditional reconstruction. GP and information: specified controls. DCAE: descriptor, not universal law.')+'};',r'\end{tikzpicture}']
  (out/'visual_abstract.tex').write_text('\n'.join(lines)+'\n')
 print('TikZ ES/EN generated from epistemic model / generado desde el modelo epistémico')
if __name__=='__main__':main()
