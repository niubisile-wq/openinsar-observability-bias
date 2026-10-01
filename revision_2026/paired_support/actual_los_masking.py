"""S14: hidden-class estimators with explicit availability and bounds."""
import os
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from masking_estimators import geometry,retained_statistics,evaluate,summarize

HERE=Path(os.environ.get('INSAR_PAIRED_ROOT',Path(__file__).resolve().parent))
STUDY=Path(os.environ['INSAR_STRENGTHENING_ROOT'])
SRC=STUDY/'results/timeseries/matched_quality.npz'
OUT=HERE/'actual_los_masking'
SEED=20260930
REPS=100
KS=(5,10,20)
FRACS=(.10,.30)
THRESHOLDS=(-2.,-5.,-10.)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    with np.load(SRC) as z:
        p=z['population'].astype(float);v=z['full_vel_los'].astype(float)
        valid=z['full_valid'].astype(bool);std=z['full_vstd'].astype(float)
    valid &= np.isfinite(v)&np.isfinite(std)&(p>=0)
    geometries={(k,dy,dx):geometry(p,valid,k,dy,dx) for k in KS
        for dy,dx in ((0,0),(k//2,0),(0,k//2),(k//2,k//2))}
    rng=np.random.default_rng(SEED)
    fields={
        'random':lambda:rng.standard_normal(p.shape),
        'spatial_cluster':lambda:gaussian_filter(rng.standard_normal(p.shape),8),
        'high_velocity_spread':lambda:std+.25*gaussian_filter(rng.standard_normal(p.shape),8),
        'low_velocity_spread':lambda:-std+.25*gaussian_filter(rng.standard_normal(p.shape),8)}
    ids=np.flatnonzero(valid);rows=[]
    for rep in range(REPS):
        for mechanism,make_score in fields.items():
            score=make_score();order=ids[np.argsort(score.ravel()[ids],kind='stable')]
            for frac in FRACS:
                missing=np.zeros(p.size,bool);missing[order[-int(round(len(order)*frac)):]]=True
                missing=missing.reshape(p.shape)
                for (k,dy,dx),g in geometries.items():
                    stats=retained_statistics(g,v,missing)
                    for threshold in THRESHOLDS:
                        for row in evaluate(g,v,missing,threshold,stats):
                            rows.append(dict(replicate=rep,mechanism=mechanism,missing_fraction=frac,
                                los_threshold_mm_yr=threshold,block_pixels=k,origin_y_pixels=dy,origin_x_pixels=dx,**row))
        if (rep+1)%10==0:print(f'S14 masks {rep+1}/{REPS}',flush=True)
    df=pd.DataFrame(rows)
    keys=['mechanism','missing_fraction','los_threshold_mm_yr','block_pixels','origin_y_pixels','origin_x_pixels','method']
    df.to_csv(OUT/'replicates.csv',index=False)
    summarize(df,keys).to_csv(OUT/'summary.csv',index=False)
    selected=df[(df.missing_fraction==.3)&(df.los_threshold_mm_yr==-5)&(df.block_pixels==10)]
    paper=summarize(selected,['mechanism','method'])
    paper.to_csv(OUT/'paper_table_minus5_1km_30pct.csv',index=False)
    meta=dict(source_sha256=hashlib.sha256(SRC.read_bytes()).hexdigest(),seed=SEED,replicates=REPS,
        thresholds_mm_yr=list(THRESHOLDS),missing_fractions=list(FRACS),block_sizes_pixels=list(KS),
        row_count=len(df),estimand='Hidden LOS-class allocation inside an artificially masked originally valid domain.',
        undefined_rule='No retained population for prevalence or no retained pixel for centre: undefined, never zero-imputed.',
        common_domain='Both estimators computable; errors normalized to originally valid population in those units.',
        unestimated_denominator='Total allocation in the original final-valid domain.',
        relative_error='Mean of per-case relative errors; zero reference cases undefined and explicitly counted.',
        bounds='Retained full-domain class allocation to retained class allocation plus all hidden allocation.',
        change_from_v051='Removes zero fallback in empty units; masks, seed, rates and original valid domain unchanged.',
        limits='Conditional allocation benchmark, not physical accuracy or evidence in genuinely missing locations.')
    (OUT/'method.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print(paper.to_string(index=False),flush=True)

if __name__=='__main__':main()
