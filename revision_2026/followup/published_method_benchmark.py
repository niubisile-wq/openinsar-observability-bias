"""S22: median-rule adaptation for a common hidden-class allocation target."""
from pathlib import Path
import hashlib,json,os,sys
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'paired_support'))
from masking_estimators import geometry,retained_statistics,evaluate,summarize

STUDY=Path(os.environ['INSAR_STRENGTHENING_ROOT'])
OUT=Path(os.environ.get('INSAR_PAIRED_ROOT',Path(__file__).resolve().parent))/'published_method_benchmark'
SRC=STUDY/'results/timeseries/matched_quality.npz'
SEED=20260930
REPLICATES=100
MASK_FRACTION=.30
BLOCK_PIXELS=10
THRESHOLD=-5.

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    with np.load(SRC) as z:
        pop=z['population'].astype(float);rate=z['full_vel_los'].astype(float)
        valid=z['full_valid'].astype(bool);spread=z['full_vstd'].astype(float)
    valid &= np.isfinite(rate)&np.isfinite(spread)&(pop>=0)
    designs={(dy,dx):geometry(pop,valid,BLOCK_PIXELS,dy,dx) for dy,dx in ((0,0),(5,0),(0,5),(5,5))}
    rng=np.random.default_rng(SEED)
    patterns={
        'random':lambda:rng.standard_normal(pop.shape),
        'spatial_cluster':lambda:gaussian_filter(rng.standard_normal(pop.shape),8),
        'high_velocity_spread':lambda:spread+.25*gaussian_filter(rng.standard_normal(pop.shape),8),
        'low_velocity_spread':lambda:-spread+.25*gaussian_filter(rng.standard_normal(pop.shape),8)}
    ids=np.flatnonzero(valid);rows=[]
    for rep in range(REPLICATES):
        for mechanism,make_score in patterns.items():
            score=make_score();order=ids[np.argsort(score.ravel()[ids],kind='stable')]
            missing=np.zeros(pop.size,bool);missing[order[-int(round(len(order)*MASK_FRACTION)):]]=True
            missing=missing.reshape(pop.shape)
            for (dy,dx),g in designs.items():
                stats=retained_statistics(g,rate,missing,include_median=True)
                for row in evaluate(g,rate,missing,THRESHOLD,stats,include_median=True):
                    rows.append(dict(replicate=rep,mechanism=mechanism,missing_fraction=MASK_FRACTION,
                        los_threshold_mm_yr=THRESHOLD,block_pixels=BLOCK_PIXELS,origin_y_pixels=dy,origin_x_pixels=dx,**row))
        if (rep+1)%20==0:print(f'S22 masks {rep+1}/{REPLICATES}',flush=True)
    df=pd.DataFrame(rows);summary=summarize(df,['mechanism','method'])
    df.to_csv(OUT/'replicates.csv',index=False);summary.to_csv(OUT/'summary.csv',index=False)
    meta=dict(source_sha256=hashlib.sha256(SRC.read_bytes()).hexdigest(),seed=SEED,replicates=REPLICATES,
        missing_fraction=MASK_FRACTION,threshold_mm_per_year=THRESHOLD,block_size_pixels=BLOCK_PIXELS,
        origins=[[0,0],[5,0],[0,5],[5,5]],row_count=len(df),
        target='Population-model allocation in the artificially hidden LOS class, not all residents of median-class blocks.',
        median_adaptation='Multiply hidden allocation by observed unit-median class; adapts rather than duplicates whole-block assignment.',
        undefined_rule='Missing required retained data produces undefined units, not zero class allocation.',
        scoring='Identical common-computable units for all three methods; signed and local L1 errors use their original valid population.',
        availability='Unestimated allocation and worst-case bounds use the original full valid-domain denominator.',
        change_from_v051='No-zero-fallback correction; seed, masks, population, LOS and grouping designs retained.',
        limitations=['Model reference is not household truth.','Original valid-domain LOS only.','Mask ranges are sensitivity, not confidence intervals.'])
    (OUT/'method.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print(summary.to_string(index=False),flush=True)

if __name__=='__main__':main()
