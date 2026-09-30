"""E2 same-chain mask association and held-out-pair predictive residuals.

Run with .venv_licsbas Python for h5py. Neither residuals nor bootstrap spread
are independent ground truth for vertical motion or permanently missing pixels.
"""
import json,re
from datetime import datetime
from pathlib import Path
import numpy as np,pandas as pd,h5py
from affine import Affine
from scipy.stats import spearmanr
from common import ROOT,OUT,dump,edges,pixel_areas,separable_resample

DEST=OUT/'timeseries'

def readraw(path,shape):return np.fromfile(path,dtype='float32').reshape(shape)

def supports(folder,pairs,shape):
    total=np.zeros(shape,dtype='float64');validtotal=np.zeros(shape,dtype='float64')
    for pair in pairs:
        unw=readraw(folder/pair/(pair+'.unw'),shape);cc=np.fromfile(folder/pair/(pair+'.cc'),dtype='uint8').reshape(shape)
        v=np.isfinite(unw)&(unw!=0);validtotal+=v;total+=v&(cc>=77)
    return (total/len(pairs)).astype('float32'),(validtotal/len(pairs)).astype('float32')

def main():
    geom=json.loads((DEST/'geometry.json').read_text(encoding='utf-8'));shape=tuple(geom['shape']);tr=Affine(*geom['transform']);x,y=edges(tr,shape)
    hold=json.loads((DEST/'holdout_design.json').read_text(encoding='utf-8'))
    allpairs=sorted(hold['training']+hold['held_out'])
    matched,unwfull=supports(DEST/'GEOC_full',allpairs,shape)
    train_support,unwtrain=supports(DEST/'GEOC_train',hold['training'],shape)
    g=np.load(OUT/'fine/grid.npz');old=np.load(OUT/'fine/mean_support.npz')['q30']
    legacy=separable_resample(old.astype('float64'),g['x_edges'],g['y_edges'],x,y)[::-1].copy()
    population=separable_resample(g['population'],g['x_edges'],g['y_edges'],x,y,normalization='source')[::-1].copy()
    area=pixel_areas(x,y)[::-1].copy()
    metrics={}
    for branch in ['full','train']:
        path=DEST/('TS_'+branch)/'results'
        metrics[branch]={k:readraw(path/k,shape) for k in ['mask','vel','vstd','stc','n_gap','maxTlen','n_loop_err','n_ifg_noloop','resid_rms','n_unw','coh_avg']}
    full_valid=np.isfinite(metrics['full']['mask'])&(metrics['full']['mask']>0)&np.isfinite(metrics['full']['vel'])
    train_valid=np.isfinite(metrics['train']['mask'])&(metrics['train']['mask']>0)&np.isfinite(metrics['train']['vel'])
    confusion=[];bins=[]
    for name,s in [('legacy28',legacy),('matched437',matched)]:
        for threshold in [.2,.3,.5,.7,.8]:
            passed=s>=threshold
            for pred in [False,True]:
                for actual in [False,True]:
                    sel=(passed==pred)&(full_valid==actual)
                    confusion.append(dict(support=name,threshold=threshold,predicted_supported=pred,final_valid=actual,pixels=int(sel.sum()),population=float(population[sel].sum()),area_m2=float(area[sel].sum())))
        for lo,hi in zip([0,.2,.4,.6,.8],[.2,.4,.6,.8,1.0000001]):
            sel=(s>=lo)&(s<hi);n=int(sel.sum())
            if not n:continue
            row=dict(support=name,bin_lower=lo,bin_upper=min(hi,1),pixels=n,population=float(population[sel].sum()),final_valid_pixel_fraction=float(full_valid[sel].mean()),final_valid_population_fraction=float(population[sel&full_valid].sum()/population[sel].sum()) if population[sel].sum()>0 else np.nan)
            for metric in ['vstd','resid_rms','stc','n_gap','maxTlen']:
                vals=metrics['full'][metric][sel];vals=vals[np.isfinite(vals)];row[metric+'_median']=float(np.median(vals)) if len(vals) else np.nan
            bins.append(row)
    pd.DataFrame(confusion).to_csv(DEST/'support_finalmask_confusion.csv',index=False);pd.DataFrame(bins).to_csv(DEST/'support_quality_bins.csv',index=False)
    residual_ss=np.zeros(shape);residual_n=np.zeros(shape,dtype='int32');pair_rows=[]
    model_ss={name:np.zeros(shape) for name in ['time_series','linear_velocity','zero_displacement']};model_n=np.zeros(shape,dtype='int32')
    with h5py.File(DEST/'TS_train/cum.h5','r') as f:
        dates=[str(int(d)) for d in f['imdates'][:]];cum=f['cum'][:];training_velocity=f['vel'][:]
        ref=f['refarea'][()];ref=ref.decode() if isinstance(ref,bytes) else str(ref)
    x1,x2,y1,y2=map(int,re.split('[:/]',ref))
    coefficient=-299792458/5405000000/(4*np.pi)*1000
    for pair in hold['held_out']:
        aa,bb=pair.split('_');row=dict(pair=pair)
        if aa not in dates or bb not in dates:
            row.update(status='excluded_missing_training_date');pair_rows.append(row);continue
        phase=readraw(DEST/'GEOC_full'/pair/(pair+'.unw'),shape);good=np.isfinite(phase)&(phase!=0)
        refphase=np.where(good[y1:y2,x1:x2],phase[y1:y2,x1:x2],np.nan)
        if not np.isfinite(refphase).any():row.update(status='excluded_invalid_training_reference');pair_rows.append(row);continue
        observed=(phase-float(np.nanmean(refphase)))*coefficient
        predicted=cum[dates.index(bb)]-cum[dates.index(aa)]
        residual=observed-predicted;valid=good&np.isfinite(predicted)
        residual_ss[valid]+=residual[valid].astype('float64')**2;residual_n[valid]+=1
        dt=(datetime.strptime(bb,'%Y%m%d')-datetime.strptime(aa,'%Y%m%d')).days/365.25
        linear=training_velocity*dt;common=valid&np.isfinite(linear);model_n[common]+=1
        for name,prediction in [('time_series',predicted),('linear_velocity',linear),('zero_displacement',np.zeros(shape))]:
            error=observed[common].astype('float64')-prediction[common].astype('float64');model_ss[name][common]+=error**2
            row[name+'_common_rmse_mm']=float(np.sqrt(np.mean(error**2))) if common.any() else np.nan
        row['common_prediction_pixels']=int(common.sum())
        row.update(status='evaluated',pixels=int(valid.sum()),rmse_mm=float(np.sqrt(np.mean(residual[valid]**2))),median_residual_mm=float(np.median(residual[valid])))
        accepted=valid&train_valid
        row['valid_mask_pixels']=int(accepted.sum());row['valid_mask_rmse_mm']=float(np.sqrt(np.mean(residual[accepted]**2))) if accepted.any() else np.nan
        pair_rows.append(row)
    rmse=np.sqrt(np.divide(residual_ss,residual_n,out=np.full(shape,np.nan),where=residual_n>0))
    heldbins=[]
    for lo,hi in zip([0,.2,.4,.6,.8],[.2,.4,.6,.8,1.0000001]):
        sel=(train_support>=lo)&(train_support<hi)&(residual_n>0)
        if not sel.any():continue
        heldbins.append(dict(support='training_pairs_only',lower=lo,upper=min(hi,1),pixels=int(sel.sum()),residual_observations=int(residual_n[sel].sum()),pooled_rmse_mm=float(np.sqrt(residual_ss[sel].sum()/residual_n[sel].sum())),median_pixel_rmse_mm=float(np.median(rmse[sel])),population=float(population[sel].sum())))
    pd.DataFrame(pair_rows).to_csv(DEST/'heldout_pairs.csv',index=False);pd.DataFrame(heldbins).to_csv(DEST/'heldout_support_bins.csv',index=False)
    comparison=[]
    groups=[('all_common',model_n>0),('training_mask_pass',(model_n>0)&train_valid),('training_mask_fail',(model_n>0)&~train_valid)]
    groups += [(f'support_{lo:.1f}_{min(hi,1):.1f}',(model_n>0)&(train_support>=lo)&(train_support<hi)) for lo,hi in zip([0,.2,.4,.6,.8],[.2,.4,.6,.8,1.0000001])]
    for group,sel in groups:
        n=int(model_n[sel].sum())
        for name,ss in model_ss.items():comparison.append(dict(group=group,model=name,pixels=int(sel.sum()),observations=n,pooled_rmse_mm=float(np.sqrt(ss[sel].sum()/n)) if n else np.nan))
    pd.DataFrame(comparison).to_csv(DEST/'heldout_prediction_comparison.csv',index=False)
    correlations=[]
    for support_name,s in [('legacy28',legacy),('matched437',matched),('training_only',train_support)]:
        for metric,v in [('full_vstd',metrics['full']['vstd']),('full_resid_rms',metrics['full']['resid_rms']),('heldout_rmse',rmse)]:
            valid=np.isfinite(v)&np.isfinite(s)
            if valid.sum()>2 and np.std(s[valid])>0 and np.std(v[valid])>0:correlations.append(dict(support=support_name,metric=metric,pixels=int(valid.sum()),spearman_r=float(spearmanr(s[valid],v[valid]).statistic)))
    pd.DataFrame(correlations).to_csv(DEST/'descriptive_correlations.csv',index=False)
    with h5py.File(DEST/'TS_full/cum.h5','r') as f:full_velocity=f['vel'][:]
    np.savez_compressed(DEST/'matched_quality.npz',x_edges=x,y_edges=y,population=population,area_m2=area,legacy_support=legacy,matched_support=matched,training_support=train_support,full_valid=full_valid,training_valid=train_valid,heldout_rmse_mm=rmse,heldout_count=residual_n,full_vel_los=full_velocity,full_vstd=metrics['full']['vstd'],full_residual=metrics['full']['resid_rms'])
    evaluated=[r for r in pair_rows if r['status']=='evaluated'];subset={}
    for name,sel in [('all_predictions',residual_n>0),('training_mask_pass',(residual_n>0)&train_valid),('training_mask_fail',(residual_n>0)&~train_valid)]:
        n=int(residual_n[sel].sum());subset[name]=dict(pixels=int(sel.sum()),observations=n,pooled_rmse_mm=float(np.sqrt(residual_ss[sel].sum()/n)) if n else None)
    summary=dict(shape=list(shape),dates_in_training=len(dates),planned_input_pairs=len(allpairs),heldout_pairs=len(hold['held_out']),evaluated_heldout_pairs=len(evaluated),excluded_heldout_pairs=[r for r in pair_rows if r['status']!='evaluated'],population_total=float(population.sum()),final_valid_population=float(population[full_valid].sum()),final_valid_population_fraction=float(population[full_valid].sum()/population.sum()),final_valid_pixel_fraction=float(full_valid.mean()),mean_matched_support=float(matched.mean()),mean_legacy_support=float(legacy.mean()),training_reference=ref,holdout=subset,interpretation='Mask association is descriptive and partially shares coherence information with the processing algorithm. Held-out edges were excluded from training; training-only support and training mask are used for predictive strata. LOS displacement residuals are not external velocity/vertical-motion accuracy.',row_order='north-to-south arrays, south-to-north y_edges')
    dump(DEST/'matched_summary.json',summary)
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
