"""Recompute every held-out edge and training-only strata from retained inputs."""
from pathlib import Path
from datetime import datetime
import json,re,hashlib,subprocess
import numpy as np,pandas as pd,h5py,networkx as nx,rasterio
from affine import Affine

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;P=ROOT/'results/timeseries';OUT=HERE/'evidence'
def raw(path,shape,dtype='float32'):return np.fromfile(path,dtype=dtype).reshape(shape)
def main():
    OUT.mkdir(exist_ok=True);h=json.loads((P/'holdout_design.json').read_text());geom=json.loads((P/'geometry.json').read_text());shape=tuple(geom['shape'])
    training=set(h['training']);held=set(h['held_out']);allpairs=training|held
    assert len(training)==393 and len(held)==44 and not(training&held)
    assert {p.name for p in (P/'GEOC_train').iterdir() if p.is_dir()}==training
    graph=nx.Graph();graph.add_edges_from(p.split('_') for p in sorted(allpairs));chosen=[]
    for pair in np.random.default_rng(h['seed']).permutation(sorted(allpairs)):
        a,b=pair.split('_');graph.remove_edge(a,b)
        if nx.is_connected(graph):chosen.append(pair)
        else:graph.add_edge(a,b)
        if len(chosen)==44:break
    assert set(chosen)==held
    passes=np.zeros(shape,dtype='int32');input_hashes=[]
    for pair in sorted(training):
        for ext in ['unw','cc']:
            a=P/'GEOC_train'/pair/(pair+'.'+ext);b=P/'GEOC_full'/pair/(pair+'.'+ext)
            digest=hashlib.sha256(a.read_bytes()).hexdigest();assert digest==hashlib.sha256(b.read_bytes()).hexdigest()
            input_hashes.append(dict(pair=pair,kind=ext,sha256=digest))
        phase=raw(P/'GEOC_train'/pair/(pair+'.unw'),shape);cc=raw(P/'GEOC_train'/pair/(pair+'.cc'),shape,'uint8')
        passes+=(np.isfinite(phase)&(phase!=0)&(cc.astype(float)/255>=.3))
    support=passes/393
    kept=np.load(P/'matched_quality.npz');assert np.max(np.abs(support-kept['training_support']))<1e-7
    with h5py.File(P/'TS_train/cum.h5','r') as f:
        dates=[str(int(d)) for d in f['imdates'][:]];cum=f['cum'][:].astype(float);velocity=f['vel'][:].astype(float)
        ref=f['refarea'][()];ref=ref.decode() if isinstance(ref,bytes) else str(ref)
    x0,x1,y0,y1=map(int,re.split('[:/]',ref))
    train_mask=raw(P/'TS_train/results/mask',shape);train_vel=raw(P/'TS_train/results/vel',shape)
    accepted=np.isfinite(train_mask)&(train_mask>0)&np.isfinite(train_vel)
    assert np.array_equal(accepted,kept['training_valid'])
    full_mask=raw(P/'TS_full/results/mask',shape);full_vel=raw(P/'TS_full/results/vel',shape)
    full_valid=np.isfinite(full_mask)&(full_mask>0)&np.isfinite(full_vel)
    assert np.array_equal(full_valid,kept['full_valid'])
    names=['time_series','linear_velocity','zero_displacement'];squares={k:np.zeros(shape) for k in names};counts=np.zeros(shape,dtype='int32');records=[]
    factor=-299792458/(5405000000*4*np.pi)*1000
    for pair in sorted(held):
        a,b=pair.split('_');assert a in dates and b in dates
        phase=raw(P/'GEOC_full'/pair/(pair+'.unw'),shape).astype(float);valid=np.isfinite(phase)&(phase!=0)
        rp=phase[y0:y1,x0:x1];rv=valid[y0:y1,x0:x1];assert rv.any()
        observed=(phase-rp[rv].mean())*factor
        predicted=cum[dates.index(b)]-cum[dates.index(a)]
        common=valid&np.isfinite(predicted)&np.isfinite(velocity);counts+=common
        years=(datetime.strptime(b,'%Y%m%d')-datetime.strptime(a,'%Y%m%d')).days/365.25
        row=dict(pair=pair,observations=int(common.sum()))
        for name,pred in [('time_series',predicted),('linear_velocity',velocity*years),('zero_displacement',np.zeros(shape))]:
            error=observed[common]-pred[common];squares[name][common]+=error**2;row[name+'_rmse']=float(np.sqrt(np.mean(error**2)))
        records.append(row)
    groups={'all_common':counts>0,'training_mask_pass':(counts>0)&accepted,'training_mask_fail':(counts>0)&~accepted}
    for lo,hi in zip([0,.2,.4,.6,.8],[.2,.4,.6,.8,1.0000001]):groups[f'support_{lo:.1f}_{min(hi,1):.1f}']=(counts>0)&(support>=lo)&(support<hi)
    rows=[];old=pd.read_csv(P/'heldout_prediction_comparison.csv')
    for group,sel in groups.items():
        n=int(counts[sel].sum())
        for model,ss in squares.items():
            rmse=float(np.sqrt(ss[sel].sum()/n));r=old[(old.group==group)&(old.model==model)].iloc[0]
            assert n==r.observations and int(sel.sum())==r.pixels
            assert abs(rmse-r.pooled_rmse_mm)<1e-4,(group,model,rmse,r.pooled_rmse_mm)
            rows.append(dict(group=group,model=model,pixels=int(sel.sum()),observations=n,rmse_mm=rmse,retained_rmse_mm=r.pooled_rmse_mm,difference=rmse-r.pooled_rmse_mm))
    cohort=[]
    for branch,pairs in [('full',allpairs),('train',training)]:
        bad=set()
        for n in ['11bad_ifg.txt','12bad_ifg.txt']:
            bad.update(re.findall(r'\b\d{8}_\d{8}\b',(P/f'TS_{branch}/info'/n).read_text()))
        assert bad<=pairs;cohort.append(dict(branch=branch,inputs=len(pairs),rejected=len(bad),inversion_pairs=len(pairs-bad)))
        assert len(pairs-bad)==(431 if branch=='full' else 382)
        for stage in ['11_check_unw','12_loop_closure','13_sb_inv','14_vel_std','15_mask_ts','16_filt_ts']:
            rec=json.loads((P/'logs'/f'{branch}_{stage}.success.json').read_text());assert rec['returncode']==0
            assert all('GEOC_full' not in a for a in rec['command']) if branch=='train' else True
        assert (P/f'TS_{branch}/cum.h5').exists() and (P/f'TS_{branch}/cum_filt.h5').exists()
    git=subprocess.run(['git','diff','--name-only'],cwd=ROOT/'external/LiCSBAS2',capture_output=True,text=True,check=True);assert not git.stdout.strip()
    head=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT/'external/LiCSBAS2',capture_output=True,text=True,check=True).stdout.strip();assert head==geom['upstream_commit']
    pd.DataFrame(rows).to_csv(OUT/'holdout_independent_groups.csv',index=False);pd.DataFrame(records).to_csv(OUT/'holdout_independent_pairs.csv',index=False)
    pd.DataFrame(input_hashes).to_csv(OUT/'training_input_hashes.csv',index=False)
    report=dict(cohort=cohort,holdout_design_recreated=True,all_44_edges_recomputed=True,training_inputs_identical_to_full_copies=True,
        no_heldout_edges_in_training=True,training_support_recomputed=True,reference_from_training_h5=ref,upstream_commit=head,tracked_upstream_source_unchanged=True,
        max_group_rmse_difference_mm=max(abs(r['difference']) for r in rows),tolerance_mm=1e-4,
        tolerance_reason='Independent float64 phase/reference arithmetic compared with original float32 intermediate subtraction; much smaller than 0.001 mm reported rounding',
        pooling='Each evaluable pixel-edge residual receives equal weight; all models use identical common observations; strata from training inputs only',
        limitation='Acquisition dates and atmospheric components are shared; internal edge reconstruction, not independent physical accuracy or new-date prediction')
    (OUT/'holdout_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
