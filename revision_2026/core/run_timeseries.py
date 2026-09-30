"""Run upstream LiCSBAS steps 11--16 with fixed settings and resumable logs."""
from pathlib import Path
import subprocess,json,os,time,argparse
from common import ROOT,OUT,dump

DEST=OUT/'timeseries';LOG=DEST/'logs';LOG.mkdir(parents=True,exist_ok=True)
PY=ROOT/'.venv_licsbas/Scripts/python.exe'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--branch',choices=['full','train','both'],default='both');args=ap.parse_args()
    env=os.environ.copy();env.update(PYTHONUTF8='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MPLBACKEND='Agg')
    for branch in (['full','train'] if args.branch=='both' else [args.branch]):
        inp=DEST/('GEOC_'+branch);out=DEST/('TS_'+branch)
        if not (inp/'slc.mli.par').exists():raise ValueError('Inputs not prepared')
        steps=[('11_check_unw',['-d',str(inp),'-t',str(out)],'info/11ifg_stats.txt'),('12_loop_closure',['-d',str(inp),'-t',str(out),'--n_para','1'],'info/12ref.txt'),('13_sb_inv',['-d',str(inp),'-t',str(out),'--n_para','1','--mem_size','512'],'cum.h5'),('14_vel_std',['-t',str(out),'--mem_size','512'],'results/vstd'),('15_mask_ts',['-t',str(out),'--noautoadjust'],'results/mask'),('16_filt_ts',['-t',str(out),'--n_para','1'],'cum_filt.h5')]
        for name,arguments,expected in steps:
            success=LOG/(branch+'_'+name+'.success.json')
            if success.exists() and (out/expected).exists():continue
            command=[str(PY),str(ROOT/'licsbas_portable.py'),'LiCSBAS'+name+'.py',*arguments]
            print('START',branch,name,flush=True);started=time.time()
            with (LOG/(branch+'_'+name+'.log')).open('w',encoding='utf-8') as log:
                result=subprocess.run(command,cwd=DEST,env=env,stdout=log,stderr=subprocess.STDOUT)
            if result.returncode or not (out/expected).exists():raise RuntimeError(f'{branch} {name} failed: return={result.returncode}, expected={expected}; inspect log')
            dump(success,dict(command=command,returncode=result.returncode,seconds=time.time()-started,expected_output=str(out/expected)))
            print('COMPLETE',branch,name,round(time.time()-started,1),'s',flush=True)
    print('LICSBAS PROCESSING COMPLETE',flush=True)

if __name__=='__main__':main()
