"""Recompute held-out residuals from raw phases and the fixed training solution."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import h5py

ROOT=Path(__file__).resolve().parent.parent
STUDY=ROOT.parent/'01_当前返修工作区/InSAR_revision_20260926/strengthening'
TS=STUDY/'results/timeseries'
OUT=ROOT/'inputs/training_quality'
OUT.mkdir(parents=True,exist_ok=True)
HASHES={}

def raw(path,shape,dtype='float32'):
    HASHES[path.relative_to(STUDY).as_posix()]=hashlib.sha256(path.read_bytes()).hexdigest()
    return np.fromfile(path,dtype=dtype).reshape(shape)

def main():
    geom=json.loads((TS/'geometry.json').read_text(encoding='utf-8'))
    shape=tuple(geom['shape'])
    hold=json.loads((TS/'holdout_design.json').read_text(encoding='utf-8'))
    assert not set(hold['training'])&set(hold['held_out'])
    source=TS/'TS_train/cum.h5'
    HASHES['results/timeseries/TS_train/cum.h5']=hashlib.sha256(source.read_bytes()).hexdigest()
    with h5py.File(source,'r') as f:
        dates=[str(int(d)) for d in f['imdates'][:]];cum=f['cum'][:]
        ref=f['refarea'][()];ref=ref.decode() if isinstance(ref,bytes) else str(ref)
    x1,x2,y1,y2=map(int,re.split('[:/]',ref))
    metrics={key:raw(TS/'TS_train/results'/key,shape) for key in
             ['coh_avg','n_unw','vstd','stc','n_gap','maxTlen','n_ifg_noloop','n_loop_err','resid_rms','mask','vel']}
    old=np.load(TS/'matched_quality.npz')
    coeff=-299792458/5405000000/(4*np.pi)*1000
    squared=[];valids=[];pair_ids=[]
    for pair in hold['held_out']:
        aa,bb=pair.split('_')
        if aa not in dates or bb not in dates: raise RuntimeError('Held-out date absent from training')
        phase=raw(TS/'GEOC_full'/pair/(pair+'.unw'),shape)
        good=np.isfinite(phase)&(phase!=0)
        reference=np.where(good[y1:y2,x1:x2],phase[y1:y2,x1:x2],np.nan)
        if not np.isfinite(reference).any(): raise RuntimeError('Held-out reference unavailable')
        obs=(phase-float(np.nanmean(reference)))*coeff
        predicted=cum[dates.index(bb)]-cum[dates.index(aa)]
        valid=good&np.isfinite(predicted)
        ss=np.zeros(shape,dtype=np.float32)
        ss[valid]=(obs[valid]-predicted[valid])**2
        squared.append(ss);valids.append(valid);pair_ids.append(pair)
    ss=np.stack(squared);vs=np.stack(valids)
    n=vs.sum(axis=0)
    rmse=np.sqrt(np.divide(ss.sum(axis=0,dtype=np.float64),n,out=np.full(shape,np.nan),where=n>0))
    assert np.array_equal(n,old['heldout_count'])
    difference=float(np.nanmax(np.abs(rmse-old['heldout_rmse_mm'])))
    assert difference<2e-5, difference
    np.savez_compressed(OUT/'quality_curve_inputs.npz',**{'train_'+k:v for k,v in metrics.items()},
                        population=old['population'],area_m2=old['area_m2'],
                        training_support=old['training_support'],heldout_pair_squared_residual=ss,
                        heldout_pair_valid=vs,heldout_pair_ids=np.array(pair_ids))
    metadata={'source':'Fixed 393-pair training inversion and 44 held-out edges in connected 437-pair stack',
              'dates':dates,'heldout_pairs':pair_ids,'training_pairs':hold['training'],
              'training_reference':ref,'los_mm_per_radian':coeff,'geometry':geom,
              'residual_agreement_with_frozen_S3_max_abs_mm':difference,'source_sha256':HASHES,
              'definition':'Squared LOS displacement prediction residuals; no external physical ground truth',
              'checks':{'training_heldout_edges_disjoint':True,'all_heldout_dates_in_training':True,
                        'all_heldout_references_evaluable':True,'counts_match_frozen_S3':True}}
    (OUT/'source_provenance.json').write_bytes((json.dumps(metadata,indent=2)+'\n').encode())
    print('Prepared quality inputs',len(pair_ids),shape,'max agreement difference',difference,flush=True)

if __name__=='__main__': main()
