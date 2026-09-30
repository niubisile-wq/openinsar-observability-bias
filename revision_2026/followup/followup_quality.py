import os
from pathlib import Path
from collections import defaultdict
import os, json, hashlib, time
import numpy as np
import pandas as pd
import rasterio
from affine import Affine
from rasterio.warp import transform, reproject, Resampling
from scipy.ndimage import gaussian_filter

BASE=Path(os.environ['INSAR_BASE'])
PROJECT=Path(os.environ['INSAR_PROJECT'])
SRC=PROJECT/'results/timeseries/matched_quality.npz'
WORLDPOP=PROJECT/'external/worldpop_thailand_2020_UNadj.tif'
GEOM=PROJECT/'results/timeseries/geometry.json'
OUT=BASE/'analysis_results'
OUT.mkdir(exist_ok=True)
SEED=20260930
REPS=30
SCALES_M=(250,500,1000,2000)
THRESHOLDS=(-2.0,-5.0,-10.0)
BLOCKS=(5,10,20)
ORIGIN_FRACS=(0.0,0.5)

def load_inputs():
    with np.load(SRC) as z:
        d={k:z[k].copy() for k in ['x_edges','y_edges','population','area_m2','matched_support','full_valid','full_vel_los','full_vstd']}
    return d

def coordinates(d):
    x=(d['x_edges'][:-1]+d['x_edges'][1:])/2
    y=(d['y_edges'][:-1]+d['y_edges'][1:])/2
    xx,yy=np.meshgrid(x,y)
    east,north=transform('EPSG:4326','EPSG:32647',xx.ravel().tolist(),yy.ravel().tolist())
    return np.asarray(east).reshape(xx.shape),np.asarray(north).reshape(xx.shape)

def labels_from_metric(x,y,size,ox,oy):
    ix=np.floor((x-ox)/size).astype(np.int64)
    iy=np.floor((y-oy)/size).astype(np.int64)
    pairs=np.column_stack((iy.ravel(),ix.ravel()))
    _,inv=np.unique(pairs,axis=0,return_inverse=True)
    return inv.astype(np.int32),int(inv.max()+1)

def group_sum(label,values,n):
    return np.bincount(label,weights=np.asarray(values).ravel(),minlength=n)

def weighted_median_los(label,v,observed,n):
    ser=pd.Series(np.where(observed.ravel(),v.ravel(),np.nan))
    med=ser.groupby(label,sort=False).median()
    out=np.full(n,np.nan)
    out[med.index.to_numpy(dtype=np.int64)]=med.to_numpy(dtype=float)
    return out

def reproject_worldpop(d):
    frozen=PROJECT/'results/timeseries/worldpop_matched.npz'
    if frozen.is_file() and not WORLDPOP.is_file():
        with np.load(frozen) as z:dst=z['worldpop_average_to_matched_grid'].copy()
        if dst.shape!=d['population'].shape:raise ValueError('Frozen WorldPop matched-grid shape mismatch')
        return dst
    geom=json.loads(GEOM.read_text(encoding='utf-8'))
    dst=np.full(d['population'].shape,np.nan,dtype=np.float64)
    tx=Affine(*geom['transform'])
    with rasterio.open(WORLDPOP) as src:
        reproject(source=rasterio.band(src,1),destination=dst,
                  src_transform=src.transform,src_crs=src.crs,src_nodata=src.nodata,
                  dst_transform=tx,dst_crs=geom['crs'],dst_nodata=np.nan,
                  resampling=Resampling.average,init_dest_nodata=True)
    return dst

def main():
    start=time.time(); d=load_inputs()
    p0=d['population'].astype(float); a=d['area_m2'].astype(float)
    smean=d['matched_support'].astype(float); smask=d['full_valid'].astype(float)
    v=d['full_vel_los'].astype(float); std=d['full_vstd'].astype(float)
    valid=d['full_valid'].astype(bool)&np.isfinite(v)&np.isfinite(std)&np.isfinite(p0)&(p0>=0)&np.isfinite(a)&(a>0)
    x,y=coordinates(d)
    wp_raw=reproject_worldpop(d)
    common=np.isfinite(wp_raw)&(wp_raw>=0)&np.isfinite(p0)&(p0>=0)&np.isfinite(a)&(a>0)
    gh_total=float(p0[common].sum())
    wp_total=float(wp_raw[common].sum())
    if gh_total<=0 or wp_total<=0: raise ValueError('Nonpositive common-domain population')
    pops={'GHSL':p0.copy(),'WorldPop_normalized':wp_raw*(gh_total/wp_total)}
    for arr in pops.values(): arr[~common]=0.0
    total_pop={name:float(q[common].sum()) for name,q in pops.items()}
    pop_model_rows=[]
    for name,q in pops.items():
        pop_model_rows.append(dict(model=name,common_pixels=int(common.sum()),common_GHSL_allocation_units=gh_total,
            native_total_before_normalization=(gh_total if name=='GHSL' else wp_total),
            normalized_total=float(q[common].sum()),worldpop_resampling='Rasterio average to the 0.001-degree matched grid; normalized to GHSL on common coverage',
            independence_note='Alternative model sensitivity only; not independent population truth. Census/ancillary source overlap and resampling effects remain.'))
    pd.DataFrame(pop_model_rows).to_csv(OUT/'population_model_common_domain.csv',index=False)

    # Factorial allocation sensitivity and covariance bound. Pixel centres are
    # assigned to metric cells; report the approximation and compare the 1-km
    # GHSL cases with the independent direct-intersection audit already on file.
    fact=[]; bounds=[]
    supports={'437_pair_mean_input':smean,'final_quality_mask':smask}
    for size in SCALES_M:
      for fy in ORIGIN_FRACS:
       for fx in ORIGIN_FRACS:
        label,n=labels_from_metric(x,y,size,fx*size,fy*size)
        area_g=group_sum(label,a,n)
        for model,q in pops.items():
         qtot=total_pop[model]
         pop_g=group_sum(label,q,n)
         density=q/a
         dbar=np.divide(group_sum(label,density*a,n),area_g,out=np.zeros(n),where=area_g>0)
         varrho=np.divide(group_sum(label,a.ravel()*(density.ravel()-dbar[label])**2,n),area_g,out=np.zeros(n),where=area_g>0)
         sigrho=np.sqrt(np.maximum(varrho,0))
         for support_name,ss in supports.items():
          asum=group_sum(label,a*ss,n)
          sbar=np.divide(asum,area_g,out=np.zeros(n),where=area_g>0)
          native_supported=float(np.sum(q*ss))
          uniform_supported=float(np.sum(pop_g*sbar))
          native_def=(qtot-native_supported)/qtot
          uniform_def=(qtot-uniform_supported)/qtot
          delta=uniform_def-native_def
          cov_num=float(native_supported-uniform_supported)
          cov_resid=delta-cov_num/qtot
          sigmasq=np.divide(group_sum(label,a.ravel()*(ss.ravel()-sbar[label])**2,n),area_g,out=np.zeros(n),where=area_g>0)
          local_bound=float(np.sum(area_g*sigrho*np.sqrt(np.maximum(sigmasq,0)))/qtot)
          local_signed=np.abs(pop_g*sbar-group_sum(label,q*ss,n))/qtot
          local_error_abs=float(local_signed.sum())
          top=max(1,int(np.ceil(n*.1)))
          top_share=float(np.sort(local_signed)[-top:].sum()/local_error_abs) if local_error_abs>0 else 0.0
          fact.append(dict(scale_m=size,origin_x_fraction=fx,origin_y_fraction=fy,population_model=model,support=support_name,
             domain_population=qtot,native_deficit_pct=100*native_def,uniform_deficit_pct=100*uniform_def,
             uniform_minus_native_pp=100*delta,absolute_local_error_share_pct=100*local_error_abs,
             top_decile_of_cells_share_of_absolute_error_pct=100*top_share,
             covariance_identity_pp=100*cov_num/qtot,covariance_identity_residual_pp=100*cov_resid,
             cauchy_schwarz_upper_bound_pp=100*local_bound,
             bound_to_absolute_effect_ratio=(local_bound/abs(delta) if abs(delta)>1e-15 else np.nan),
             integration='UTM projected pixel-centre assignment; partial source-domain pixels retained by source-cell area',
             truth_status='Conditional on each population model and support field; no household truth'))
          # Signed cell-level aggregation effect components for diagnostics.
          signed=(pop_g*sbar-group_sum(label,q*ss,n))/qtot
          bounds.append(dict(scale_m=size,origin_x_fraction=fx,origin_y_fraction=fy,population_model=model,support=support_name,
             cells=n,positive_effect_cell_share=float(np.mean(signed>0)),negative_effect_cell_share=float(np.mean(signed<0)),
             max_positive_cell_effect_pp=float(max(0,np.max(signed))*100),max_negative_cell_effect_pp=float(min(0,np.min(signed))*100),
             mean_absolute_cell_effect_pp=float(np.mean(np.abs(signed))*100)))
    pd.DataFrame(fact).to_csv(OUT/'allocation_covariance_bounds.csv',index=False)
    pd.DataFrame(bounds).to_csv(OUT/'allocation_local_effects.csv',index=False)

    # Controlled mask experiment: compare the published within-block median
    # velocity assignment with the population-weighted within-block prevalence.
    # Ground truth is known only because artificial masking is imposed on pixels
    # that already had valid fitted LOS velocity; no conclusion about original
    # invalid locations is drawn.
    rng=np.random.default_rng(SEED)
    mechanisms=('random','spatial_cluster','high_velocity_spread','low_velocity_spread')
    case_rows=[]
    ids=np.flatnonzero(valid)
    for rep in range(REPS):
      for mechanism in mechanisms:
        if mechanism=='random': score=rng.standard_normal(p0.shape)
        elif mechanism=='spatial_cluster': score=gaussian_filter(rng.standard_normal(p0.shape),8)
        elif mechanism=='high_velocity_spread': score=std+0.25*gaussian_filter(rng.standard_normal(p0.shape),8)
        else: score=-std+0.25*gaussian_filter(rng.standard_normal(p0.shape),8)
        order=ids[np.argsort(score.ravel()[ids],kind='stable')]
        nmiss=int(round(len(ids)*.30)); missing=np.zeros(p0.size,dtype=bool); missing[order[-nmiss:]]=True; missing=missing.reshape(p0.shape)
        observed=valid&~missing
        for k in BLOCKS:
          for fy in ORIGIN_FRACS:
           for fx in ORIGIN_FRACS:
            label,n=labels_from_metric(np.indices(p0.shape)[1],np.indices(p0.shape)[0],k,fx*k,fy*k)
            # The function is generic in its coordinates; here units are source pixels.
            obs_pop={name:group_sum(label,q*observed,n) for name,q in pops.items()}
            miss_pop={name:group_sum(label,q*missing,n) for name,q in pops.items()}
            median=weighted_median_los(label,v,observed,n)
            for threshold in THRESHOLDS:
              target=v<=threshold
              obs_class={name:group_sum(label,q*observed*target,n) for name,q in pops.items()}
              true_hidden={name:float(np.sum(q*missing*target)) for name,q in pops.items()}
              median_flag=np.isfinite(median)&(median<=threshold)
              for name,q in pops.items():
                prevalence=np.divide(obs_class[name],obs_pop[name],out=np.zeros(n),where=obs_pop[name]>0)
                estimates={'within_block_population_prevalence':prevalence*miss_pop[name],
                           'unweighted_within_block_median_velocity':median_flag*miss_pop[name]}
                denominator=float(np.sum(q[valid]))
                for method,est in estimates.items():
                  err=float(est.sum()-true_hidden[name])
                  abs_cell=float(np.abs(est-(group_sum(label,q*missing*target,n))).sum())
                  case_rows.append(dict(replicate=rep,mechanism=mechanism,missing_fraction=.30,block_pixels=k,
                    origin_x_fraction=fx,origin_y_fraction=fy,population_model=name,los_threshold_mm_yr=threshold,
                    method=method,true_artificially_hidden_class_population=true_hidden[name],
                    estimated_artificially_hidden_class_population=float(est.sum()),signed_bias=err,
                    signed_bias_share_of_observed_domain_population_pct=100*err/denominator,
                    absolute_block_error_share_of_observed_domain_population_pct=100*abs_cell/denominator,
                    blocks_without_training_observation=int(np.sum((obs_pop[name]<=0)&(miss_pop[name]>0))),
                    data_domain='pixels passing original LiCSBAS final full-valid mask only',
                    validity='Semisynthetic conditional truth; artificial mask only; not original invalid pixels or household counts'))
      if (rep+1)%10==0: print(f'mask replicates {rep+1}/{REPS}',flush=True)
    cf=pd.DataFrame(case_rows)
    cf.to_csv(OUT/'block_median_baseline_replicates.csv',index=False)
    keys=['mechanism','block_pixels','origin_x_fraction','origin_y_fraction','population_model','los_threshold_mm_yr','method']
    summary=cf.groupby(keys,dropna=False).agg(replicates=('replicate','count'),true_hidden_population=('true_artificially_hidden_class_population','mean'),
       estimated_hidden_population=('estimated_artificially_hidden_class_population','mean'),mean_signed_bias=('signed_bias','mean'),
       sd_signed_bias=('signed_bias','std'),q025_signed_bias=('signed_bias',lambda z:z.quantile(.025)),
       q975_signed_bias=('signed_bias',lambda z:z.quantile(.975)),
       mean_absolute_block_error_share_pct=('absolute_block_error_share_of_observed_domain_population_pct','mean'),
       max_blocks_without_training=('blocks_without_training_observation','max')).reset_index()
    summary.to_csv(OUT/'block_median_baseline_summary.csv',index=False)

    worldpop_used=WORLDPOP if WORLDPOP.is_file() else PROJECT/'results/timeseries/worldpop_matched.npz'
    meta=dict(source=str(SRC),source_sha256=hashlib.sha256(SRC.read_bytes()).hexdigest(),
      worldpop_source=str(worldpop_used),worldpop_sha256=hashlib.sha256(worldpop_used.read_bytes()).hexdigest(),
      seed=SEED,replicates=REPS,projected_scales_m=list(SCALES_M),block_source_pixel_sizes=list(BLOCKS),
      mask_mechanisms=list(mechanisms),mask_fraction=.30,thresholds_mm_yr=list(THRESHOLDS),origins='zero and half-unit in both dimensions',
      population_models='GHSL native; WorldPop average-resampled to processing grid and normalized to GHSL total over common cells; no independent-truth claim',
      support_estimand='Difference in fixed-domain population-weighted mean deficit between native population weighting and area-uniform support within projected units',
      geometry_approximation='Projected pixel-centre grouping for 250-2000 m sensitivity. Compare 1-km GHSL effect estimates with prior direct polygon intersection audit; differences are not to be pooled.',
      covariance_diagnostic='Compute exact identity for grouped centre-assigned data and Cauchy-Schwarz bound. Mathematical bound uses model-specific population density and area-weighted support spread.',
      block_baseline='Unweighted median LOS velocity within visible pixels in artificial square reporting units; assigns all artificially hidden population in a unit to the median category. This is not actual Census block geography.',
      artificial_truth='Velocity known in all originally final-valid pixels before simulated removal. Original failed pixels excluded.',
      interval_interpretation='2.5th-97.5th percentiles across seeded artificial masks are design sensitivity ranges, not confidence intervals.',
      limitations=['No independent household or building occupancy truth','No deformation truth in original invalid pixels','WorldPop and GHSL lineage overlap may occur','Projected centre-bin approximation','Square blocks approximate but do not reproduce real Census polygons'])
    (OUT/'followup_method.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(dict(allocation_factorial_rows=len(fact),local_effect_rows=len(bounds),mask_baseline_replicate_rows=len(cf),mask_summary_rows=len(summary),seconds=round(time.time()-start,1),out=str(OUT)),ensure_ascii=False))
if __name__=='__main__': main()

