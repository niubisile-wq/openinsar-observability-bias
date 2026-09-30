import os
"""Controlled artificial masking of observed LiCSBAS LOS and population.

This estimates model-based population counts in velocity classes on the
observed-valid domain only. It says nothing about motion at genuinely invalid
pixels or physical hazard exposure.
"""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

HERE = Path(os.environ.get('INSAR_PAIRED_ROOT', Path(__file__).resolve().parent))
STUDY=Path(os.environ['INSAR_STRENGTHENING_ROOT'])
SRC=STUDY/'results'/'timeseries'/'matched_quality.npz'
OUT=HERE/'actual_los_masking'
SEED=20260930
REPS=100
KS=(5,10,20)       # nominally ~0.5, 1, 2 km at this latitude
FRACS=(.10,.30)
THRESHOLDS=(-2.,-5.,-10.)  # mm/yr, negative LOS (away from satellite)

def main():
    OUT.mkdir(exist_ok=True)
    with np.load(SRC) as z:
        p=z['population'].astype(float); v=z['full_vel_los'].astype(float)
        valid=z['full_valid'].astype(bool); std=z['full_vstd'].astype(float)
    valid &= np.isfinite(v)&np.isfinite(std)&(p>=0)
    rng=np.random.default_rng(SEED); rows=[]
    # Four artificial missingness mechanisms, each applied only to pixels with
    # observed valid LOS. Spatial fields are continuous across the full grid.
    fields={
        'random': lambda: rng.standard_normal(p.shape),
        'spatial_cluster': lambda: gaussian_filter(rng.standard_normal(p.shape),8),
        # Add a seeded spatial perturbation so repeated masks vary while still
        # preferentially removing the highest/lowest uncertainty pixels.
        'high_velocity_spread': lambda: std + 0.25*gaussian_filter(rng.standard_normal(p.shape),8),
        'low_velocity_spread': lambda: -std + 0.25*gaussian_filter(rng.standard_normal(p.shape),8),
    }
    for rep in range(REPS):
      for mechanism,make_score in fields.items():
        score=make_score(); ids=np.flatnonzero(valid); order=ids[np.argsort(score.ravel()[ids],kind='stable')]
        for frac in FRACS:
          missing=np.zeros(p.size,bool); missing[order[-int(round(len(order)*frac)):]]=True; missing=missing.reshape(p.shape)
          observed=valid&~missing
          for threshold in THRESHOLDS:
            target=v<=threshold
            true=float((p*missing*target).sum()); hiddenp=float((p*missing).sum())
            for k in KS:
              h,w=p.shape
              for dy,dx in [(0,0),(k//2,0),(0,k//2),(k//2,k//2)]:
                  yy,xx=np.indices((h,w)); gy=(yy+dy)//k; gx=(xx+dx)//k
                  labels=(gy*(int(gx.max())+1)+gx).ravel(); n=int(labels.max()+1)
                  def agg(a): return np.bincount(labels,weights=a.ravel(),minlength=n)
                  pp=p; mm=missing; oo=observed; vv=target
                  misspop=agg(pp*mm); truemiss=agg(pp*mm*vv)
                  opop=agg(pp*oo); otarget=agg(pp*oo*vv)
                  # Within-block observed prevalence, extrapolated only to the
                  # artificially masked population in that same block.
                  pred=np.divide(otarget,opop,out=np.zeros(n),where=opop>0)*misspop
                  # Nearest observed pixel(s) to each actual block center,
                  # vectorized; ties are averaged to avoid arbitrary choices.
                  lab2=labels.reshape(h,w); rr,cc=np.indices((h,w)); lab=labels
                  rmin=np.full(n,h,dtype=int); rmax=np.full(n,-1,dtype=int)
                  cmin=np.full(n,w,dtype=int); cmax=np.full(n,-1,dtype=int)
                  np.minimum.at(rmin,lab,rr.ravel()); np.maximum.at(rmax,lab,rr.ravel())
                  np.minimum.at(cmin,lab,cc.ravel()); np.maximum.at(cmax,lab,cc.ravel())
                  cy=(rmin+rmax)/2; cx=(cmin+cmax)/2
                  d2=(rr.ravel()-cy[lab])**2+(cc.ravel()-cx[lab])**2
                  mind=np.full(n,np.inf); np.minimum.at(mind,lab,np.where(oo.ravel(),d2,np.inf))
                  tie=oo.ravel()&(d2==mind[lab]); ntie=np.bincount(lab,weights=tie.astype(float),minlength=n)
                  nclass=np.bincount(lab,weights=(tie&vv.ravel()).astype(float),minlength=n)
                  near=np.divide(nclass,ntie,out=np.zeros(n),where=ntie>0)*misspop
                  for method,estimate in [('within_block_observed_population_prevalence',pred),('nearest_observed_pixel_to_block_center',near)]:
                    err=estimate-truemiss
                    rows.append(dict(replicate=rep,mechanism=mechanism,missing_fraction=frac,los_threshold_mm_yr=threshold,
                        block_pixels=k,origin_y_pixels=dy,origin_x_pixels=dx,method=method,observed_valid_pixels=int(valid.sum()),artificially_masked_population=hiddenp,
                        true_masked_class_population=true,estimated_masked_class_population=float(estimate.sum()),
                        signed_error=float(err.sum()),absolute_error_population=float(abs(err).sum()),
                        absolute_block_error_fraction_observed_population=float(abs(err).sum()/p[valid].sum())))
      if (rep+1)%20==0: print('replicates',rep+1,flush=True)
    df=pd.DataFrame(rows); df.to_csv(OUT/'replicates.csv',index=False)
    keys=['mechanism','missing_fraction','los_threshold_mm_yr','block_pixels','origin_y_pixels','origin_x_pixels','method']
    summary=df.groupby(keys).agg(replicates=('replicate','count'),true_population=('true_masked_class_population','mean'),
      mean_estimate=('estimated_masked_class_population','mean'),mean_signed_error=('signed_error','mean'),
      sd_signed_error=('signed_error','std'),q025_signed_error=('signed_error',lambda x:x.quantile(.025)),
      q975_signed_error=('signed_error',lambda x:x.quantile(.975)),
      mean_absolute_block_error_share=('absolute_block_error_fraction_observed_population','mean')).reset_index()
    summary.to_csv(OUT/'summary.csv',index=False)
    meta=dict(source=str(SRC),source_sha256=hashlib.sha256(SRC.read_bytes()).hexdigest(),seed=SEED,replicates=REPS,
      actual_los_velocity_convention='LiCSBAS: positive toward satellite; negative away; relative LOS mm/yr',
      thresholds_mm_yr=list(THRESHOLDS),missing_fractions=list(FRACS),mechanisms=list(fields),block_sizes_pixels=list(KS),
      origins='four zero/half-block offsets in each axis; partial raster-edge blocks retained',
      valid_domain='Final full-valid mask only; artificial masking imposed only on otherwise valid pixels',
      spread_definition='LiCSBAS bootstrap velocity spread diagnostic; not calibrated physical uncertainty',
      estimand='Population-weighted pixels with observed LOS <= threshold, artificially hidden inside originally valid domain',
      interpretation='Conditional, semisynthetic test of estimating LOS-class population within the observed-valid domain. Not vertical velocity, hazard, a household count, or validation of genuinely unobserved pixels.',
      known_limitation='The original validity mask and LOS estimates may share processing dependencies; this experiment does not infer values in invalid pixels.')
    (OUT/'method.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False),encoding='utf-8')
    print(summary.to_string(index=False))
if __name__=='__main__': main()
