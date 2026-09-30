"""E3 semisynthetic masking with known exposure and synthetic hazard truth.

The fine allocation is a controlled reference, never independent population
validation. Synthetic velocity fields do not validate real unobserved motion.
"""
import argparse
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.stats import rankdata
from common import OUT,CFG,dump

DEST=OUT/'masking';DEST.mkdir(exist_ok=True)

def grouping(shape,k,shift):
    yy,xx=np.indices(shape)
    gx=(xx+shift)//k;gy=(yy+shift)//k
    labels=(gy*(gx.max()+1)+gx).ravel();n=int(labels.max()+1)
    def agg(a):return np.bincount(labels,weights=a.ravel(),minlength=n)
    # Geometric center sample, including clipped edge cells: use the center of
    # their actual included footprint, so no out-of-domain extrapolation.
    xc=np.rint(agg(xx)/agg(np.ones(shape))).astype(int)
    yc=np.rint(agg(yy)/agg(np.ones(shape))).astype(int)
    return labels,agg,yc,xc

def standard(a):return (a-a.mean())/(a.std()+1e-12)

def main():
    global DEST
    ap=argparse.ArgumentParser();ap.add_argument('--region',choices=['chao','dwr'],default='chao');args=ap.parse_args()
    if args.region=='chao':
        z=np.load(OUT/'fine/grid.npz');xx=(z['x_edges'][:-1]+z['x_edges'][1:])/2;yy=(z['y_edges'][:-1]+z['y_edges'][1:])/2
        ix=np.where((xx>=100.35)&(xx<100.75))[0];iy=np.where((yy>=13.55)&(yy<13.95))[0]
        p=z['population'][np.ix_(iy,ix)];a=z['area_m2'][np.ix_(iy,ix)];b=z['builtup_m2'][np.ix_(iy,ix)]
        d=z['domain'][np.ix_(iy,ix)];rectangle=[100.35,13.55,100.75,13.95]
    else:
        DEST=OUT/'dwr/masking';DEST.mkdir(exist_ok=True)
        z=np.load(OUT/'dwr/grid.npz');xx=(z['x_edges'][:-1]+z['x_edges'][1:])/2;yy=(z['y_edges'][:-1]+z['y_edges'][1:])/2
        ix=np.arange(len(xx));iy=np.arange(len(yy));p=z['population']*z['support'];a=z['area_m2']*z['support'];b=z['builtup_m2']*z['support'];d=a>0;rectangle=[-119.9,36.6,-119.7,36.85]
    # The masking experiment has its own fixed rectangular domain. Retain GHSL
    # cells even when the unrelated archived velocity footprint has a hole.
    active=a>0;density=np.log1p(np.divide(p,a,out=np.zeros_like(p),where=active)*1e6);density=(density-density[active].mean())/(density[active].std()+1e-12);rng=np.random.default_rng(CFG['seed'])
    independent=standard(gaussian_filter(rng.normal(size=p.shape),15))
    # Rates in mm/yr are constructed experimental variables with known truth.
    hazard={'density_positive':density>0,'density_negative':density<0,'independent':independent>0}
    velocities={k:np.where(v,-10.,0.) for k,v in hazard.items()}
    np.savez_compressed(DEST/'reference_fields.npz',population=p,builtup_m2=b,area_m2=a,legacy_vlm_domain=d,lon=xx[ix],lat=yy[iy],**velocities)
    groups={}
    for k in [5,10,20]:
        for shift in [0,k//2]:
            lab,agg,yc,xc=grouping(p.shape,k,shift)
            groups[(k,shift)]=(agg,yc,xc,agg(p),agg(a),agg(b))
    rows=[];hazrows=[]
    for replicate in range(CFG['synthetic_replicates']):
        white=standard(rng.normal(size=p.shape));spatial=standard(gaussian_filter(rng.normal(size=p.shape),8))
        scores={'random':white,'spatial':spatial,'density_positive':density+0.5*spatial,'density_negative':-density+0.5*spatial}
        for mechanism,score in scores.items():
            activeids=np.flatnonzero(active);order=activeids[np.argsort(score.ravel()[activeids],kind='stable')]
            for fraction in CFG['synthetic_missing_fractions']:
                m=np.zeros(p.size,dtype=bool);m[order[-int(round(len(order)*fraction)):]]=True;m=m.reshape(p.shape)
                reference=float((p*m).sum());total=float(p.sum());areafrac=float((a*m).sum()/a.sum())
                for (k,shift),(agg,yc,xc,P,A,B) in groups.items():
                    cellref=agg(p*m);cellarea=np.divide(agg(a*m),A,out=np.zeros_like(A),where=A>0)
                    buildingfrac=np.divide(agg(b*m),B,out=cellarea.copy(),where=B>0)
                    for method,estimate in [('area_uniform',P*cellarea),('building_weighted',P*buildingfrac),('cell_center',P*(m[yc,xc]|~active[yc,xc]))]:
                        err=estimate-cellref
                        rows.append(dict(replicate=replicate,mechanism=mechanism,missing_pixel_fraction=fraction,missing_area_fraction=areafrac,cell_size_fine_pixels=k,origin_fine_pixels=shift,method=method,reference=reference,total_population=total,estimate=float(estimate.sum()),signed_error=float(err.sum()),relative_error=float(err.sum()/reference),absolute_cell_error_fraction_total=float(np.abs(err).sum()/total),rmse_cells=float(np.sqrt(np.mean(err**2)))))
                for name,h in hazard.items():
                    truth=float((p*h).sum());observed=float((p*h*~m).sum());hidden=truth-observed
                    # The interval is an algebraic identification bound, not a
                    # fitted estimator or a nominal-confidence uncertainty band.
                    if hidden<-1e-7 or hidden>reference+1e-7:raise ValueError('Identification bound')
                    hazrows.append(dict(replicate=replicate,mechanism=mechanism,missing_pixel_fraction=fraction,hazard=name,truth=truth,observed=observed,hidden=hidden,unsupported_exposure=reference,hidden_share_of_unsupported=hidden/reference,bound_lower=observed,bound_upper=observed+reference))
        if (replicate+1)%10==0:print('MASKING',replicate+1,flush=True)
    df=pd.DataFrame(rows);df.to_csv(DEST/'replicates.csv',index=False)
    group=['mechanism','missing_pixel_fraction','cell_size_fine_pixels','origin_fine_pixels','method']
    summary=df.groupby(group).agg(replicates=('replicate','count'),mean_relative_error=('relative_error','mean'),sd_relative_error=('relative_error','std'),q025_relative_error=('relative_error',lambda x:x.quantile(.025)),q975_relative_error=('relative_error',lambda x:x.quantile(.975)),mean_cell_absolute_error_fraction=('absolute_cell_error_fraction_total','mean')).reset_index()
    summary.to_csv(DEST/'summary.csv',index=False)
    pd.DataFrame(hazrows).to_csv(DEST/'synthetic_hazard_bounds.csv',index=False)
    dump(DEST/'method.json',dict(region=args.region,seed=CFG['seed'],replicates=CFG['synthetic_replicates'],rectangle=rectangle,shape=list(p.shape),domain='Fixed rectangle for Chao; published observed DWR footprint within fixed rectangle for Fresno',missing_fractions=CFG['synthetic_missing_fractions'],mechanisms=list(scores),reference='GHSL 2020 native 3ss population model, held fixed; known introduced missingness; DWR includes fractional published support before introducing masks',grids='5, 10, 20 native pixels with zero/half-pixel-count origins; angular grids, not exact metric grids',center_rule='Fine-pixel introduced mask at included rectangle cell center; center outside observed domain has zero support',building_weighting='Redistribute each coarse-cell population in proportion to GHSL built area, fallback to included area when built area is zero',hazard_truth='Synthetic two-valued -10 or 0 mm/yr; density positive/negative or independent spatial field',interpretation='Monte Carlo intervals describe introduced-mask variation conditional on fixed allocation models, not uncertainty in real exposure or motion',no_physical_validation=True))
    print('E3 COMPLETE',flush=True)

if __name__=='__main__':main()
