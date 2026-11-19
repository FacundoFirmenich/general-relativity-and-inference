"""Build bilingual self-contained notebook and visual software payloads."""
from pathlib import Path
import base64,hashlib,io,json,zipfile,nbformat
R=Path(__file__).resolve().parents[1]
def selected():
 keep=[]
 for p in R.rglob('*'):
  rel=p.relative_to(R)
  if not p.is_file() or any((x.startswith('.') and x not in ['.streamlit','.github']) or x=='__pycache__' for x in rel.parts):continue
  if rel.parts[0] in ['.streamlit','.github'] and str(rel) not in ['.streamlit/config.toml','.github/workflows/reproduce.yml']:continue
  if p.name.startswith('GRAVITAS_') or p.suffix in ['.aux','.out']:continue
  if p.suffix=='.log' and not (rel.parts[0]=='records' and p.name.startswith('all_')):continue
  if any(x in rel.parts for x in ['publication','tables','figures']):continue
  if 'epistemic_export' in p.name or p.name.startswith('epistemic_') and rel.parts[0]=='records':continue
  if rel.parts[0]=='records' and ('standalone' in p.name or 'epistemic_export' in p.name):continue
  keep.append(p)
 return sorted(keep)
EXTRACT='''from pathlib import Path
import base64,hashlib,zipfile,io,tempfile,subprocess,sys

def extract():
    raw=base64.b64decode(PAYLOAD)
    if hashlib.sha256(raw).hexdigest()!=PAYLOAD_SHA256:raise RuntimeError('Carga alterada / Payload integrity failure')
    work=Path(tempfile.mkdtemp(prefix='GRAVITAS_'))
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        for item in archive.infolist():
            if not (work/item.filename).resolve().is_relative_to(work):raise RuntimeError('Ruta insegura / Unsafe archive path')
        archive.extractall(work)
    return work
'''
def main():
 b=io.BytesIO()
 with zipfile.ZipFile(b,'w',zipfile.ZIP_DEFLATED,compresslevel=9) as z:
  for p in selected():
   info=zipfile.ZipInfo(str(p.relative_to(R)),date_time=(2026,11,18,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=(0o755 if p.suffix in ['.bash','.command'] else 0o644)<<16;z.writestr(info,p.read_bytes())
 payload=b.getvalue();enc=base64.b64encode(payload).decode();digest=hashlib.sha256(payload).hexdigest();prefix='PAYLOAD_SHA256='+repr(digest)+'\nPAYLOAD='+repr(enc)+'\n'
 standalone='#!/usr/bin/env python3\n'+prefix+EXTRACT+'''\nif __name__=='__main__':
    work=extract()
    print('Carpeta / Workspace:',work,flush=True)
    mode=sys.argv[1] if len(sys.argv)>1 else 'ui'
    if mode!='extract':raise SystemExit(subprocess.call([sys.executable,str(work/'launch_gravitas.py'),mode],cwd=work))
'''
 (R/'GRAVITAS_STANDALONE.py').write_text(standalone)
 for lang in ['es','en']:
  def T(es,en):return es if lang=='es' else en
  nb=nbformat.v4.new_notebook();nb.metadata.update(kernelspec={'name':'python3','display_name':'Python 3','language':'python'},language_info={'name':'python'},colab={'name':f'GRAVITAS_{lang.upper()}.ipynb','provenance':[]})
  nb.cells=[nbformat.v4.new_markdown_cell(T('# Relatividad general y métodos de inferencia\n## GRAVITAS · cuaderno autónomo\n\nFacundo Firmenich · [ORCID](https://orcid.org/0009-0002-6578-3811) · f.firmenich@cedesur.org\n\nEjecutar todas extrae una carga verificada, instala dependencias y reproduce cálculos, tablas y figuras. No exige GitHub ni Drive; la primera instalación necesita Internet. Los manuscritos y su información suplementaria se entregan por separado y no están incrustados en este cuaderno.\n\nEl grafo ofrece tour, play, pausa, fuentes descargables y límites de alcance. Los cálculos Cassini son condicionales a entradas reconstruidas. No hay nueva cota física ni validación instrumental independiente. iPhone/iPad acceden mediante Colab en el navegador, no mediante un ejecutable nativo.','# General Relativity and Inference Methods\n## GRAVITAS · self-contained notebook\n\nFacundo Firmenich · [ORCID](https://orcid.org/0009-0002-6578-3811) · f.firmenich@cedesur.org\n\nRun all extracts a verified payload, installs dependencies and reproduces calculations, tables and figures. Neither GitHub nor Drive is required; first installation needs Internet. Manuscripts and supplementary information are delivered separately and are not embedded in this notebook.\n\nThe graph offers a tour, play, pause, downloadable sources and claim boundaries. Cassini calculations are conditional on reconstructed inputs. No new physical bound or independent instrumental validation is claimed. iPhone/iPad access uses Colab in a browser, not a native executable.'))]
  nb.cells[0].source+='\n\n'+T('Software archivado v1.11: ','Archived software v1.11: ')+'[10.5281/zenodo.24381485](https://doi.org/10.5281/zenodo.24381485) · '+T('DOI general: ','Concept DOI: ')+'[10.5281/zenodo.24381484](https://doi.org/10.5281/zenodo.24381484).\n'
  setup=prefix+EXTRACT+'''\nWORK=extract()
from importlib.metadata import version,PackageNotFoundError
needed=[]
for line in (WORK/'requirements.txt').read_text().splitlines()+['ipywidgets==8.1.7']:
    name,wanted=line.split('==')
    try:
        if version(name)!=wanted:needed.append(line)
    except PackageNotFoundError:needed.append(line)
if needed:subprocess.check_call([sys.executable,'-m','pip','install',*needed])
print('Workspace / Carpeta:',WORK)
print('SHA-256:',PAYLOAD_SHA256)
'''
  nb.cells.append(nbformat.v4.new_code_cell(setup))
  nb.cells.append(nbformat.v4.new_markdown_cell(T('## 1 · Ejecutar cálculos, tablas y figuras\nIncluye los controles de información, Cassini, Gaia, QUMOND, DCAE y la generación del grafo. Los fallos predictivos científicos se conservan aunque termine correctamente el programa.','## 1 · Run calculations, tables and figures\nIncludes information controls, Cassini, Gaia, QUMOND, DCAE and graph generation. Scientific predictive failures are retained even when the program exits successfully.')))
  nb.cells.append(nbformat.v4.new_code_cell("subprocess.run([sys.executable,str(WORK/'run.py'),'all'],cwd=WORK,check=True)"))
  nb.cells.append(nbformat.v4.new_markdown_cell(T('## 2 · Explorar la investigación','## 2 · Explore the investigation')))
  nb.cells.append(nbformat.v4.new_code_cell("from IPython.display import HTML, display\nimport html\ngraph=(WORK/'ui/epistemic/epistemic_graph_"+lang+".html').read_text()\ndisplay(HTML('<iframe style=\"width:100%;height:1100px;border:0\" srcdoc=\"'+html.escape(graph,quote=True)+'\"></iframe>'))"))
  nb.cells.append(nbformat.v4.new_markdown_cell(T('## 3 · Reproducir por bloques\nElegí el modo y pulsá el botón. No se ejecuta código arbitrario desde el selector.','## 3 · Reproduce by block\nChoose a mode and press the button. The selector does not execute arbitrary code.')))
  labels=[(T('FractalGo: topología','FractalGo: topology'),'topology'),(T('Fermi-LAT','Fermi-LAT'),'fermi'),(T('Estabilidad posterior','Posterior stability'),'stability'),(T('Sólo figuras','Figures only'),'figures'),(T('Sólo cálculos','Calculations only'),'calculations'),(T('Tablas','Tables'),'tables'),(T('Información','Information'),'information'),(T('DCAE / COSMOS','DCAE / COSMOS'),'frontier'),(T('Grafo y TikZ','Graph and TikZ'),'graph')]
  code="import ipywidgets as widgets\nselector=widgets.Dropdown(options="+repr(labels)+")\nbutton=widgets.Button(description="+repr(T('Ejecutar','Run'))+")\noutput=widgets.Output()\ndef execute_block(_):\n    button.disabled=True\n    try:\n        with output:\n            output.clear_output()\n            result=subprocess.run([sys.executable,str(WORK/'run.py'),selector.value],cwd=WORK,capture_output=True,text=True)\n            print(result.stdout)\n            print(result.stderr)\n            result.check_returncode()\n    finally:button.disabled=False\nbutton.on_click(execute_block)\ndisplay(widgets.VBox([selector,button,output]))"
  nb.cells.append(nbformat.v4.new_code_cell(code))
  nb.cells.append(nbformat.v4.new_markdown_cell(T('## 4 · Guardar y descargar\nEl ZIP conserva resultados, fuentes y registros. Activá DESCARGAR en Colab para descargarlo.','## 4 · Save and download\nThe ZIP retains results, sources and records. Enable DOWNLOAD in Colab to download it.')))
  nb.cells.append(nbformat.v4.new_code_cell("import shutil\narchive=shutil.make_archive(str(WORK.parent/(WORK.name+'_outputs')),'zip',WORK)\nprint(archive)\nDOWNLOAD=False\nif DOWNLOAD:\n    from google.colab import files\n    files.download(archive)"))
  nbformat.write(nb,R/f'GRAVITAS_{lang.upper()}.ipynb')
 (R/'records/standalone_payload_v1_11.json').write_text(json.dumps({'sha256':digest,'bytes':len(payload),'files':len(selected()),'scope_es':'Carga autónoma: cálculos, interfaz y grafo; sin archivos privados de lectura.','scope_en':'Self-contained payload: calculations, UI and graph; no private reading archive.'},indent=2)+'\n')
 print(json.dumps({'payload_bytes':len(payload),'files':len(selected()),'sha256':digest}))
if __name__=='__main__':main()
