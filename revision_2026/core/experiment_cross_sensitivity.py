"""Cross population allocation with pair sampling; preserve extensive moments.

Run from strengthening. Common domain uses only GHSL/WorldPop 2020 coverage.
All four combinations share the GHSL common-domain population total. Origins
are sensitivity cases, not independent statistical replicates.
"""
from pathlib import Path
import hashlib,json,gc
import numpy as np
import pandas as pd
from affine import Affine
from rasterio.warp import reproject,Resampling
from common import OUT,ROOT
from experiment_allocation import geometry,partitions

R=ROOT
D=OUT/'cross_sensitivity'
D.mkdir(exist_ok=True)
ORIGINS=[(0,0),(500,0),(0,500),(500,500)]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    paths=[OUT/'fine/grid.npz',OUT/'population/population_models.npz',
           OUT/'fine/mean_support.npz',OUT/'sampling/mean_224.npz']
    with np.load(paths[0]) as z:g={k:z[k] for k in z.files}
    with np.load(paths[1]) as z:
        models={k:z[k] for k in ['GHSL_2020_3ss','WorldPop_2020_3ss']}
        domain=g['domain'].copy()
        for k in models:domain &= z[k+'_coverage']>1-1e-5
    fields={}
    for key,path in zip(['legacy28','expanded224'],paths[2:]):
        with np.load(path) as z:fields[key]=z['q30'].astype('float64')
    total=float(g['population'][domain].sum())
    original_totals={k:float(v[domain].sum()) for k,v in models.items()}
    models={k:np.where(domain,v,0.)*(total/original_totals[k]) for k,v in models.items()}
    area=np.where(domain,g['area_m2'],0.)
    st,dt50,x50,y50=geometry(g)
    rows=[];conservation=[]
    for resolution in [50,25]:
        x=x50[0]+np.arange(round((x50[-1]-x50[0])/resolution)+1)*resolution
        y=y50[0]+np.arange(round((y50[-1]-y50[0])/resolution)+1)*resolution
        dt=Affine(resolution,0,x[0],0,-resolution,y[-1])
        groups=[partitions(x,y,1000,sx,sy)[2] for sx,sy in ORIGINS]
        def transfer(a,label):
            dst=np.zeros((len(y)-1,len(x)-1),dtype='float64')
            reproject(a[::-1].copy(),dst,src_transform=st,src_crs='EPSG:4326',
                dst_transform=dt,dst_crs='EPSG:32647',resampling=Resampling.sum,
                src_nodata=None,dst_nodata=0,num_threads=2)
            expected=float(a.sum());actual=float(dst.sum())
            # The 25-m run diagnoses GDAL nonlinear-reprojection drift without
            # renormalizing it. Production 50-m sums retain the strict criterion.
            tolerance=max(1e-5,abs(expected)*(1e-8 if resolution==50 else 1e-4))
            if abs(actual-expected)>tolerance:raise ValueError(('conservation',label,expected,actual))
            coarse=[fn(dst[::-1]) for fn in groups]
            for k,b in enumerate(coarse):
                if abs(float(b.sum())-actual)>max(1e-5,abs(actual)*1e-8):raise ValueError(('aggregation',label,k))
            conservation.append(dict(resolution_m=resolution,quantity=label,source_sum=expected,
                projected_sum=actual,relative_drift=(actual-expected)/expected,
                max_coarse_difference=max(abs(float(b.sum())-actual) for b in coarse)))
            del dst;gc.collect()
            return coarse
        A=transfer(area,'area')
        S={sample:transfer(area*f,sample+'_area_support') for sample,f in fields.items()}
        for model,p in models.items():
            P=transfer(p,model+'_population')
            for sample,f in fields.items():
                H=transfer(p*f,model+'_'+sample+'_population_support')
                native=float(np.sum(p*(1-f)))
                for i,(sx,sy) in enumerate(ORIGINS):
                    mean=np.divide(S[sample][i],A[i],out=np.zeros_like(A[i]),where=A[i]>0)
                    native_coarse=P[i]-H[i];uniform=P[i]*(1-mean)
                    difference=uniform-native_coarse
                    covariance=H[i]-P[i]*mean
                    if not np.allclose(difference,covariance,rtol=1e-8,atol=1e-7):raise ValueError('covariance')
                    if resolution==50 and not np.isclose(native_coarse.sum(),native,rtol=1e-8,atol=1e-5):raise ValueError('native target')
                    estimate=float(uniform.sum());delta=estimate-native
                    rows.append(dict(integration_m=resolution,model=model,sample=sample,coherence_cutoff=.3,
                        cell_size_m=1000,origin_x_m=sx,origin_y_m=sy,common_population=total,
                        native_deficit=native,uniform_deficit=estimate,signed_difference=delta,
                        difference_pct_native=100*delta/native,difference_pp_total=100*delta/total,
                        native_deficit_pct_total=100*native/total,uniform_deficit_pct_total=100*estimate/total,
                        projected_native_deficit=float(native_coarse.sum()),covariance_integral=float(covariance.sum())))
                print('DONE',resolution,model,sample,flush=True)
        pd.DataFrame(rows).to_csv(D/'completed_resolutions.csv',index=False)
    df=pd.DataFrame(rows)
    df.to_csv(D/'all_resolutions.csv',index=False)
    primary=df[df.integration_m.eq(50)].copy();primary.to_csv(D/'cross_sensitivity.csv',index=False)
    keys=['model','sample','origin_x_m','origin_y_m']
    cmp=primary.merge(df[df.integration_m.eq(25)],on=keys,suffixes=('_50','_25'))
    cmp['difference_pp_total_25_minus_50']=cmp.difference_pp_total_25-cmp.difference_pp_total_50
    cmp.to_csv(D/'integration_convergence.csv',index=False)
    pd.DataFrame(conservation).to_csv(D/'conservation.csv',index=False)
    maxpp=float(cmp.difference_pp_total_25_minus_50.abs().max())
    if maxpp>.01:raise ValueError(('integration discrepancy exceeds 0.01 percentage points',maxpp))
    summary=[]
    for (model,sample),z in primary.groupby(['model','sample'],sort=False):
        summary.append(dict(model=model,sample=sample,native_deficit=float(z.native_deficit.iloc[0]),
            uniform_min=float(z.uniform_deficit.min()),uniform_max=float(z.uniform_deficit.max()),
            delta_min=float(z.signed_difference.min()),delta_max=float(z.signed_difference.max()),
            pct_min=float(z.difference_pct_native.min()),pct_max=float(z.difference_pct_native.max()),
            pp_min=float(z.difference_pp_total.min()),pp_max=float(z.difference_pp_total.max())))
    audit=dict(inputs=[dict(path=str(p.relative_to(R)),sha256=sha(p)) for p in paths],
        primary_rows=len(primary),common_cells=int(domain.sum()),common_population=total,
        excluded_ghsl_population=float(g['population'][g['domain']&~domain].sum()),
        original_common_domain_totals=original_totals,normalization_factors={k:total/v for k,v in original_totals.items()},
        max_integration_difference_pp_total=maxpp,all_signed_differences_positive=bool((primary.signed_difference>0).all()),
        projection_diagnostic='50-m production conserved at 1e-8 relative; 25-m raw diagnostic allows at most 1e-4 relative moment drift, reports it without renormalization, and must change uniform-minus-native by <0.01 percentage points of common population',
        common_domain='Inherited primary footprint intersected with full GHSL/WorldPop 2020 coverage only; not the four-layer epoch intersection',
        sampling='Legacy 28 and expanded 224 are different pair selections; this comparison jointly changes sample size and composition, not sample size alone',
        transfer='Spherical fine-pixel support; conservative UTM47N transfer of A, A*s, P and P*s separately; aggregate extensive cross moment, never multiply separate projected averages',
        interpretation='Conditional model-allocation comparison; four origins are not independent replicates; no population truth or physical missing-motion validation',
        summary=summary)
    (D/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(audit,indent=2),flush=True)

if __name__=='__main__':main()
