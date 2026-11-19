"""Allowlisted local execution with an exclusive lock and per-run records."""
from pathlib import Path
import subprocess, sys, os, json, datetime, hashlib, time
ROOT=Path(__file__).resolve().parents[1]
COMMANDS={'graph':['tools/build_epistemic_graph.py'],'platform':['platform/audit_gravitools.py'],
          'information':['research/revision_v18/information_audit.py'],
          **{k:['run.py',k] for k in ['all','calculations','figures','tables','documents','revision','fermi','stability','topology']}}

def execute(mode,report=lambda line:None):
    if mode not in COMMANDS:raise ValueError('Unknown computation mode')
    directory=ROOT/'records/ui_runs';directory.mkdir(parents=True,exist_ok=True)
    lock=directory/'running.lock'
    try:fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    except FileExistsError:raise RuntimeError('Another run holds records/ui_runs/running.lock. Wait for it to finish; after a crash inspect the PID before removing a stale lock.')
    os.write(fd,str(os.getpid()).encode());os.close(fd)
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')
    log=directory/(stamp+'_'+mode+'.log');start=time.monotonic()
    argv=[sys.executable,'-u',*COMMANDS[mode]]
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MPLBACKEND='Agg')
    try:
        with log.open('w') as stream,subprocess.Popen(argv,cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True) as process:
            for line in process.stdout:stream.write(line);stream.flush();report(line)
            status=process.wait()
        result=dict(mode=mode,argv=argv,started_utc=stamp,runtime_seconds=time.monotonic()-start,exit_status=status,status='completed' if status==0 else 'failed',log=str(log.relative_to(ROOT)),log_sha256=hashlib.sha256(log.read_bytes()).hexdigest())
        log.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
        return result
    finally:lock.unlink(missing_ok=True)
