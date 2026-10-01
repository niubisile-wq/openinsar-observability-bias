"""Prespecified training-only post-filters, evaluated on fixed held-out edges."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd

def analyse(inputs,out):
    out.mkdir(parents=True,exist_ok=True)
    z=np.load(inputs)
    p=z['population'].astype(float);a=z['area_m2'].astype(float)
    ss=z['heldout_pair_squared_residual'].astype(float);valid=z['heldout_pair_valid'].astype(bool)
    ids=z['heldout_pair_ids'].astype(str)
    base=np.isfinite(z['train_vel'])
    default=base&np.isfinite(z['train_mask'])&(z['train_mask']>0)
    assert len(ids)==44 and np.all(p>=0) and np.all(a>0)
    rules=[('finite_training_velocity','none',np.nan,base),('fixed_training_default_mask','default',np.nan,default)]
    for t in [0.05,0.2,0.3,0.4,0.5,0.6,0.7]:
        rules.append((f'coherence_min_{t:g}','coh_avg_min',t,default&np.isfinite(z['train_coh_avg'])&(z['train_coh_avg']>=t)))
    for t in [0.5,1.0,1.5,2.0]:
        rules.append((f'training_residual_max_{t:g}','resid_rms_max_mm',t,default&np.isfinite(z['train_resid_rms'])&(z['train_resid_rms']<=t)))
    for t in [0.5,1.0,2.0,5.0,100.0]:
        rules.append((f'velocity_std_max_{t:g}','vstd_max_mm_yr',t,default&np.isfinite(z['train_vstd'])&(z['train_vstd']<=t)))
    summaries=[];pair_rows=[]
    counts=valid.sum(axis=0);sums=ss.sum(axis=0);has=counts>0
    pixel_rmse=np.sqrt(np.divide(sums,counts,out=np.full(p.shape,np.nan),where=has))
    for name,axis,t,keep in rules:
        accepted=keep&has
        n=int(counts[accepted].sum());weight=float((p[accepted]*counts[accepted]).sum())
        if n==0: raise RuntimeError('Predefined rule has no evaluable pixels: '+name)
        pairs=[]
        for k,pair in enumerate(ids):
            ok=keep&valid[k];ni=int(ok.sum())
            value=float(np.sqrt(ss[k][ok].sum()/ni)) if ni else np.nan
            pairs.append(value)
            pair_rows.append({'rule':name,'pair':pair,'valid_pixels':ni,'domain_population_covered_pct':100*p[ok].sum()/p.sum(),
                              'rmse_mm':value,'population_weighted_rmse_mm':float(np.sqrt((p[ok]*ss[k][ok]).sum()/p[ok].sum())) if p[ok].sum()>0 else np.nan})
        summaries.append({'rule':name,'training_metric':axis,'threshold':t,
                          'accepted_pixels':int(keep.sum()),'accepted_domain_area_pct':100*a[keep].sum()/a.sum(),
                          'accepted_domain_population_pct':100*p[keep].sum()/p.sum(),
                          'evaluable_accepted_pixels':int(accepted.sum()),'evaluable_domain_population_pct':100*p[accepted].sum()/p.sum(),
                          'accepted_allocation_without_heldout_residual_units':float(p[keep&~has].sum()),
                          'heldout_observations':n,'pooled_rmse_mm':float(np.sqrt(sums[accepted].sum()/n)),
                          'population_weighted_pooled_rmse_mm':float(np.sqrt((p[accepted]*sums[accepted]).sum()/weight)) if weight>0 else np.nan,
                          'pixel_rmse_median_mm':float(np.median(pixel_rmse[accepted])),
                          'pixel_rmse_P95_mm':float(np.quantile(pixel_rmse[accepted],.95)),
                          'pair_rmse_median_mm':float(np.nanmedian(pairs)),
                          'pair_rmse_min_mm':float(np.nanmin(pairs)),'pair_rmse_max_mm':float(np.nanmax(pairs))})
    frame=pd.DataFrame(summaries)
    for axis,ascending in [('coh_avg_min',True),('resid_rms_max_mm',False),('vstd_max_mm_yr',False)]:
        sub=frame[frame.training_metric.eq(axis)].sort_values('threshold',ascending=ascending)
        assert np.all(np.diff(sub.accepted_domain_population_pct)<=1e-10)
    frame.to_csv(out/'quality_coverage_error.csv',index=False)
    pd.DataFrame(pair_rows).to_csv(out/'quality_coverage_error_by_pair.csv',index=False)
    protocol={'input_sha256':hashlib.sha256(inputs.read_bytes()).hexdigest(),
              'training_only_filters':True,'post_filter_of_fixed_training_solution':True,
              'heldout_edges':44,'prespecified_rules':len(rules),'population_denominator':float(p.sum()),
              'area_denominator_m2':float(a.sum()),'common_per_rule_evaluation':'Same frozen pair residuals and valid-prediction fields; report changed spatial coverage explicitly',
              'uncertainty':'Pair ranges describe shared-date, spatially dependent internal prediction errors; no iid confidence interval.',
              'limits':['Held-out edge residuals are internal LOS displacement prediction errors, not external mm/yr velocity or vertical-motion accuracy.',
                        'Filters cannot validate locations with no predictions or missing observations.',
                        'Rules were fixed before this curve was read; no optimal threshold is selected using held-out error.',
                        'No new inversion is fitted for any post-filter. The historic default training mask is retained as the base.'],
              'checks':{'training_filter_nested_coverage':True,'all_rules_have_evaluable_observations':True}}
    (out/'quality_coverage_error_protocol.json').write_bytes((json.dumps(protocol,indent=2)+'\n').encode())
    print(frame[['rule','accepted_domain_population_pct','pooled_rmse_mm','population_weighted_pooled_rmse_mm']].to_string(index=False),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--inputs',type=Path);parser.add_argument('--output',type=Path)
    args=parser.parse_args();root=Path(__file__).resolve().parent.parent
    analyse(args.inputs or root/'inputs/training_quality/quality_curve_inputs.npz',args.output or root/'results/training_quality')
