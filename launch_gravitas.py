#!/usr/bin/env python3
"""Double-click launchers call this installer and local-only dashboard launcher."""
from pathlib import Path
import sys, os, subprocess, venv, socket, time, threading, webbrowser, urllib.request
ROOT=Path(__file__).resolve().parent

def environment():
    if sys.version_info<(3,11):raise SystemExit('Se requiere Python 3.11 o posterior / Python 3.11 or newer is required. Alternativa / Alternative: Colab standalone.')
    directory=ROOT/'.gravitas-env'
    python=directory/('Scripts/python.exe' if os.name=='nt' else 'bin/python')
    if not python.exists():
        print('Preparando Gravitas / Preparing Gravitas…',flush=True)
        venv.EnvBuilder(with_pip=True).create(directory)
    # A receipt binds installed dependencies to their exact requested versions.
    import hashlib
    files=[ROOT/'requirements.txt',ROOT/'requirements-ui.txt']
    digest=hashlib.sha256(b''.join(p.read_bytes() for p in files)).hexdigest()
    marker=directory/'requirements.sha256'
    if not marker.exists() or marker.read_text()!=digest:
        print('Primera instalación: se necesita internet / First installation needs internet.',flush=True)
        subprocess.run([str(python),'-m','pip','install','-r',str(files[0]),'-r',str(files[1])],check=True)
        marker.write_text(digest)
    return python

def main():
    python=environment()
    mode=sys.argv[1] if len(sys.argv)>1 else 'ui'
    if mode!='ui':
        if mode not in ['all','calculations','figures','tables','validate','documents','dcae','frontier','revision','graph','information','fermi','stability','topology']:
            raise SystemExit('Modo desconocido / Unknown mode')
        raise SystemExit(subprocess.call([str(python),str(ROOT/'run.py'),mode],cwd=ROOT))
    if not any((ROOT/'figures').glob('*.png')):
        subprocess.run([str(python),str(ROOT/'run.py'),'figures'],cwd=ROOT,check=True)
    with socket.socket() as sock:
        sock.bind(('127.0.0.1',0));port=sock.getsockname()[1]
    url=f'http://127.0.0.1:{port}'
    def open_when_ready():
        for attempt in range(150):
            try:
                with urllib.request.urlopen(url+'/_stcore/health',timeout=1) as r:
                    if r.status==200:webbrowser.open(url);return
            except (OSError,ValueError):pass
            time.sleep(.3)
        print('Abrir / Open:',url,flush=True)
    threading.Thread(target=open_when_ready,daemon=True).start()
    print('Gravitas:',url,'\nCerrar esta ventana detiene la aplicación / Close this window to stop.',flush=True)
    try:
        subprocess.run([str(python),'-m','streamlit','run',str(ROOT/'ui/app.py'),'--server.address=127.0.0.1',f'--server.port={port}','--server.headless=true','--browser.gatherUsageStats=false','--theme.base=light','--theme.primaryColor=#795d96'],cwd=ROOT,check=True)
    except KeyboardInterrupt:pass

if __name__=='__main__':main()
