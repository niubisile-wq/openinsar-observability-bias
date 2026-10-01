"""Controlled missingness on exact official Census block/rate-cell intersections."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd

SEED=20261001
REPEATS=30

def analyse(inputs,out):
    out.mkdir(parents=True,exist_ok=True)
    z=np.load(inputs)
    rate=z['rates_mm_per_year'].ravel().astype(float);valid=z['original_valid'].ravel().astype(bool)
    ci=z['intersection_pixel'];bi=z['intersection_block'];area=z['intersection_area_m2'];weight=z['allocation_weight']
    cp=z['center_pixel'];cb=z['center_block'];nb=len(z['block_geoid'])
    sx=z['pixel_x_m'].ravel();sy=z['pixel_y_m'].ravel()
    wtotal=np.bincount(bi,weights=weight,minlength=nb);full=float(wtotal.sum());known=wtotal>0
    v_ids=np.flatnonzero(valid)
    references={t:np.bincount(bi,weights=weight*(rate[ci]<=t),minlength=nb) for t in [-3.,-5.,-10.]}
    rng=np.random.default_rng(SEED);rows=[];baseline=[]
    def score_observed(obs,rep,pattern,fraction):
        observed=obs[ci]
        ow=np.bincount(bi,weights=weight*observed,minlength=nb)
        oa=np.bincount(bi,weights=area*observed,minlength=nb)
        med=np.full(nb,np.nan)
        present=obs[cp]
        if present.any():
            grouped=pd.Series(rate[cp[present]]).groupby(cb[present]).median()
            med[grouped.index.to_numpy(int)]=grouped.to_numpy()
        natives=known&(ow>0);medok=known&np.isfinite(med)
        common=natives&medok
        if not common.any(): raise RuntimeError('No common computable blocks')
        for threshold,truth in references.items():
            oc=np.bincount(bi,weights=weight*observed*(rate[ci]<=threshold),minlength=nb)
            ac=np.bincount(bi,weights=area*observed*(rate[ci]<=threshold),minlength=nb)
            native=wtotal*np.divide(oc,ow,out=np.full(nb,np.nan),where=ow>0)
            uniform=wtotal*np.divide(ac,oa,out=np.full(nb,np.nan),where=oa>0)
            median=np.where(np.isfinite(med),wtotal*(med<=threshold),np.nan)
            hidden=float(full-ow.sum());lower=float(oc.sum());upper=lower+hidden;truth_total=float(truth.sum())
            assert lower-1e-6<=truth_total<=upper+1e-6
            for method,estimate,defined in [('native_observed_allocation_fraction',native,natives),
                                          ('area_observed_class_fraction',uniform,natives),
                                          ('published_pixel_center_median_class',median,medok)]:
                selected=defined&known
                err=estimate[selected]-truth[selected]
                common_err=estimate[common]-truth[common]
                denominator=float(wtotal[selected].sum());common_p=float(wtotal[common].sum())
                rows.append({'replicate':rep,'mask_geometry':pattern,'removed_valid_pixel_fraction':fraction,
                             'threshold_mm_per_year':threshold,'method':method,
                             'known_footprint_allocation_units':full,'known_footprint_blocks':int(known.sum()),
                             'reference_class_allocation_units':truth_total,'reference_class_share_pct':100*truth_total/full,
                             'computable_blocks':int(selected.sum()),'unestimated_blocks':int((known&~defined).sum()),
                             'unestimated_allocation_units':float(wtotal[known&~defined].sum()),
                             'unestimated_allocation_share_pct':100*wtotal[known&~defined].sum()/full,
                             'reference_class_in_computable_units':float(truth[selected].sum()),
                             'estimated_class_in_computable_units':float(estimate[selected].sum()),
                             'conditional_computable_signed_error_pp':100*err.sum()/denominator,
                             'conditional_computable_local_L1_pp':100*np.abs(err).sum()/denominator,
                             'common_computable_blocks':int(common.sum()),'common_computable_allocation_units':common_p,
                             'common_signed_error_pp':100*common_err.sum()/common_p,
                             'common_local_L1_pp':100*np.abs(common_err).sum()/common_p,
                             'common_signed_error_normalized_to_full_pp':100*common_err.sum()/full,
                             'hidden_allocation_share_pct':100*hidden/full,
                             'whole_known_domain_lower_bound_pct':100*lower/full,
                             'whole_known_domain_upper_bound_pct':100*upper/full,
                             'whole_known_domain_bound_contains_reference':True})
    # Preserve baseline method errors and missing-center availability separately.
    score_observed(valid,-1,'unmasked',0.)
    for rep in range(REPEATS):
        scores={'random_cells':rng.random(len(rate))}
        for pattern,centers in [('one_contiguous_hole',1),('four_contiguous_holes',4)]:
            cx=rng.uniform(sx[v_ids].min(),sx[v_ids].max(),centers)
            cy=rng.uniform(sy[v_ids].min(),sy[v_ids].max(),centers)
            scores[pattern]=np.min((sx[:,None]-cx)**2+(sy[:,None]-cy)**2,axis=1)
        for pattern,score in scores.items():
            order=v_ids[np.argsort(score[v_ids],kind='stable')]
            for fraction in [.1,.3,.5]:
                obs=valid.copy();obs[order[:int(round(len(order)*fraction))]]=False
                score_observed(obs,rep,pattern,fraction)
        if rep%10==0: print('Controlled real-block replicate',rep+1,flush=True)
    frame=pd.DataFrame(rows)
    frame[frame.replicate.eq(-1)].to_csv(out/'real_block_unmasked_baseline.csv',index=False)
    frame=frame[frame.replicate.ge(0)]
    frame.to_csv(out/'real_block_masking_replicates.csv',index=False)
    keys=['mask_geometry','removed_valid_pixel_fraction','threshold_mm_per_year','method']
    summary=frame.groupby(keys,as_index=False).agg(
        repetitions=('replicate','count'),
        mean_unestimated_allocation_pct=('unestimated_allocation_share_pct','mean'),
        min_unestimated_allocation_pct=('unestimated_allocation_share_pct','min'),
        max_unestimated_allocation_pct=('unestimated_allocation_share_pct','max'),
        mean_common_signed_error_pp=('common_signed_error_pp','mean'),
        min_common_signed_error_pp=('common_signed_error_pp','min'),
        max_common_signed_error_pp=('common_signed_error_pp','max'),
        mean_common_local_L1_pp=('common_local_L1_pp','mean'),
        mean_conditional_computable_signed_error_pp=('conditional_computable_signed_error_pp','mean'),
        mean_whole_domain_bound_width_pp=('hidden_allocation_share_pct','mean'),
        mean_lower_bound_pct=('whole_known_domain_lower_bound_pct','mean'),
        mean_upper_bound_pct=('whole_known_domain_upper_bound_pct','mean'))
    summary.to_csv(out/'real_block_masking_summary.csv',index=False)
    protocol={'input_sha256':hashlib.sha256(inputs.read_bytes()).hexdigest(),'seed':SEED,'repetitions':REPEATS,
              'mask_fraction_definition':'Count of originally valid rate pixels, not a preset fraction of population allocation',
              'masks':'Random cell order; nearest distance to one or four random centers in EPSG:3310; nested 10/30/50% deletions',
              'whole_domain_bounds':'Known retained class allocation <= full controlled known-domain class allocation <= retained class allocation + all withheld allocation',
              'fair_comparison':'Each comparator scored against the same unmasked native class allocation on common computable real blocks. Undefined units also reported per method.',
              'native_fraction':'Observed allocation-weighted class proportion extrapolated only inside computable blocks; it is a conditional estimator, not truth for hidden cells.',
              'limits':['All references are from the original known DWR footprint; no reference for genuinely missing deformation.',
                        'The pixel-center median is discontinuous and may be undefined in small irregular blocks even before masking.',
                        'Population is a Census-total-normalized GHSL allocation, not household locations.',
                        'The deterministically correct bounds can be wide; narrow conditional errors do not establish identification.',
                        'No confidence intervals or independence assumptions across spatial repetitions.'],
              'checks':{'all_2430_controlled_method_rows':len(frame)==2430,'all_bounds_contain_reference':bool(frame.whole_known_domain_bound_contains_reference.all()),
                        'no_missing_comparator_values_replaced_by_zero':True}}
    assert protocol['checks']['all_2430_controlled_method_rows']
    (out/'real_block_masking_protocol.json').write_bytes((json.dumps(protocol,indent=2)+'\n').encode())
    print(summary[summary.removed_valid_pixel_fraction.eq(.3)&summary.threshold_mm_per_year.eq(-5.)].to_string(index=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path);p.add_argument('--output',type=Path);args=p.parse_args()
    root=Path(__file__).resolve().parent.parent
    analyse(args.inputs or root/'inputs/fresno_real_blocks/real_block_masking_inputs.npz',args.output or root/'results/real_block_masking')
