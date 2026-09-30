"""Audit retained held-out edges, full workflow completion and prediction arithmetic."""
import json,re,subprocess
import numpy as np,pandas as pd,h5py
from common import ROOT,OUT,dump

def main():
    p=OUT/'timeseries';h=json.loads((p/'holdout_design.json').read_text());g=json.loads((p/'geometry.json').read_text());shape=tuple(g['shape'])
    assert len(h['training'])==393 and len(h['held_out'])==44 and not(set(h['training'])&set(h['held_out']))
    actual={x.name for x in (p/'GEOC_train').iterdir() if x.is_dir()};assert actual==set(h['training'])
    logs=[]
    for branch in ['full','train']:
        for step in ['11_check_unw','12_loop_closure','13_sb_inv','14_vel_std','15_mask_ts','16_filt_ts']:
            q=p/'logs'/(branch+'_'+step+'.success.json');r=json.loads(q.read_text());assert r['returncode']==0;logs.append(str(q.relative_to(ROOT)))
    z=np.load(p/'matched_quality.npz');summary=json.loads((p/'matched_summary.json').read_text());conf=pd.read_csv(p/'support_finalmask_confusion.csv')
    for _,c in conf.groupby(['support','threshold']):
        assert int(c.pixels.sum())==np.prod(shape);assert np.isclose(c.population.sum(),z['population'].sum(),rtol=1e-12)
    bins=pd.read_csv(p/'heldout_support_bins.csv');assert int(bins.residual_observations.sum())==int(z['heldout_count'].sum())
    prediction=pd.read_csv(p/'heldout_prediction_comparison.csv')
    assert prediction.groupby('group').observations.nunique().eq(1).all()
    assert prediction.groupby('group').pixels.nunique().eq(1).all()
    valid=np.isfinite(z['heldout_rmse_mm']);reconstructed=float(np.sqrt(np.sum(z['heldout_rmse_mm'][valid]**2*z['heldout_count'][valid])/z['heldout_count'].sum()))
    assert np.isclose(reconstructed,summary['holdout']['all_predictions']['pooled_rmse_mm'],rtol=1e-12)
    table=pd.read_csv(p/'heldout_pairs.csv');errors=[]
    with h5py.File(p/'TS_train/cum.h5','r') as f:
        dates=f['imdates'][()].astype(str).tolist();ref=f['refarea'][()].decode();xx,xe,yy,ye=map(int,re.split('[:/]',ref));cum=f['cum'][:];velocity=f['vel'][:]
    for _,r in table[table.status.eq('evaluated')].iloc[::11].iterrows():
        a,b=r.pair.split('_');phase=np.fromfile(p/'GEOC_full'/r.pair/(r.pair+'.unw'),dtype=np.float32).reshape(shape);good=np.isfinite(phase)&phase.__ne__(0)
        reference=np.nanmean(np.where(good[yy:ye,xx:xe],phase[yy:ye,xx:xe],np.nan))
        obs=(phase-float(reference))*(-299792458/(5405000000*4*np.pi)*1000)
        pred=cum[dates.index(b)]-cum[dates.index(a)];ok=good&np.isfinite(pred)
        err=abs(float(np.sqrt(np.mean((obs[ok]-pred[ok])**2)))-r.rmse_mm);assert err<1e-5;errors.append(err)
        from datetime import datetime
        dt=(datetime.strptime(b,'%Y%m%d')-datetime.strptime(a,'%Y%m%d')).days/365.25
        common=ok&np.isfinite(velocity)
        for model,prediction_field in [('time_series',pred),('linear_velocity',velocity*dt),('zero_displacement',np.zeros(shape))]:
            error=obs[common].astype('float64')-prediction_field[common].astype('float64')
            assert abs(float(np.sqrt(np.mean(error**2)))-r[model+'_common_rmse_mm'])<1e-5
    git=subprocess.run(['git','diff','--name-only'],cwd=ROOT/'external/LiCSBAS2',capture_output=True,text=True,check=True)
    assert not git.stdout.strip(), 'Tracked upstream scientific source changed'
    head=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT/'external/LiCSBAS2',capture_output=True,text=True,check=True).stdout.strip()
    assert head=='2627e50a1c703becd4a0b18754b05c1cf8c564d5'
    dump(p/'audit.json',dict(all_twelve_stages_complete=True,heldout_input_leakage='No held-out edge occurs in training input directories',confusion_population_conservation='pass',heldout_bin_observation_conservation='pass',pooled_residual_reconstruction='pass',independently_recomputed_pair_rmse_max_difference_mm=max(errors),upstream_tracked_source_diff=git.stdout,logs=logs,scope='Implementation and design audit, not external physical accuracy'))
    print('MATCHED AUDIT PASS',flush=True)

if __name__=='__main__':main()
