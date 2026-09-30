"""Cross-check complete experiment inventories and recompute selected mechanisms.

Uses saved reference arrays and raw DWR polygons, not production calculation
helpers. This is an internal computational audit, not external validation.
"""
from pathlib import Path
import json,math
import numpy as np,pandas as pd,networkx as nx
from scipy.ndimage import gaussian_filter
from scipy.stats import spearmanr
from shapely.geometry import Polygon,box
from shapely import union_all

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;RES=ROOT/'results';OUT=HERE/'evidence'

def mask_check(region):
    folder=RES/region;df=pd.read_csv(folder/'replicates.csv');summary=pd.read_csv(folder/'summary.csv')
    keys=['mechanism','missing_pixel_fraction','cell_size_fine_pixels','origin_fine_pixels','method']
    counts=df.groupby(keys).replicate.agg(['count','nunique']);assert len(df)==21600 and (counts['count']==100).all() and (counts['nunique']==100).all()
    assert set(df.mechanism)=={'random','spatial','density_positive','density_negative'}
    assert set(df.missing_pixel_fraction)=={.1,.3,.5}
    assert set(df.cell_size_fine_pixels)=={5,10,20} and set(df.method)=={'area_uniform','building_weighted','cell_center'}
    assert np.max(np.abs(df.estimate-df.reference-df.signed_error))<1e-6
    assert np.max(np.abs(df.signed_error/df.reference-df.relative_error))<1e-9
    stats=df.groupby(keys).relative_error.agg(['mean','std']);joined=stats.join(summary.set_index(keys))
    assert np.max(np.abs(joined['mean']-joined.mean_relative_error))<1e-9
    # Recreate all RNG draws; independently calculate four selected repetitions.
    z=np.load(folder/'reference_fields.npz');p=z['population'];a=z['area_m2'];b=z['builtup_m2'];active=a>0
    density=np.log1p(np.divide(p,a,out=np.zeros_like(p),where=active)*1e6);density=(density-density[active].mean())/(density[active].std()+1e-12)
    rng=np.random.default_rng(20260927);rng.normal(size=p.shape) # synthetic independent hazard draw
    yy,xx=np.indices(p.shape);labels=((yy//10)*(math.ceil(p.shape[1]/10))+xx//10).ravel()
    n=int(labels.max()+1)
    agg=lambda v:np.bincount(labels,weights=v.ravel(),minlength=n)
    P,A,B=agg(p),agg(a),agg(b);ones=agg(np.ones(p.shape));cx=np.rint(agg(xx)/ones).astype(int);cy=np.rint(agg(yy)/ones).astype(int)
    rows=[]
    norm=lambda v:(v-v.mean())/(v.std()+1e-12)
    for rep in range(100):
        white=norm(rng.normal(size=p.shape));spatial=norm(gaussian_filter(rng.normal(size=p.shape),8))
        if rep not in {0,1,49,99}:continue
        scores={'random':white,'spatial':spatial,'density_positive':density+.5*spatial,'density_negative':-density+.5*spatial}
        ids=np.flatnonzero(active)
        for mech,score in scores.items():
            ordered=ids[np.argsort(score.ravel()[ids],kind='stable')]
            for fraction in [.1,.3,.5]:
                mask=np.zeros(p.size,dtype=bool);mask[ordered[-round(len(ordered)*fraction):]]=True;mask=mask.reshape(p.shape)
                ref=float(p[mask].sum());area_frac=np.divide(agg(a*mask),A,out=np.zeros_like(A),where=A>0)
                building_frac=np.divide(agg(b*mask),B,out=area_frac.copy(),where=B>0)
                values={'area_uniform':float(P@area_frac),'building_weighted':float(P@building_frac),
                    'cell_center':float(P@(mask[cy,cx]|~active[cy,cx]))}
                for method,estimate in values.items():
                    old=df[(df.replicate==rep)&(df.mechanism==mech)&(df.missing_pixel_fraction==fraction)&
                        (df.cell_size_fine_pixels==10)&(df.origin_fine_pixels==0)&(df.method==method)].iloc[0]
                    assert abs(ref-old.reference)<1e-6 and abs(estimate-old.estimate)<1e-6
                    rows.append(dict(region=region,replicate=rep,mechanism=mech,fraction=fraction,method=method,reference=ref,estimate=estimate,difference=estimate-old.estimate))
    bounds=pd.read_csv(folder/'synthetic_hazard_bounds.csv');assert (bounds.truth>=bounds.bound_lower-1e-6).all() and (bounds.truth<=bounds.bound_upper+1e-6).all()
    pd.DataFrame(rows).to_csv(OUT/(region.replace('/','_')+'_independent.csv'),index=False)
    main=df[(df.missing_pixel_fraction==.3)&(df.cell_size_fine_pixels==10)&(df.origin_fine_pixels==0)]
    synopsis=main.groupby(['mechanism','method']).agg(reference_mean=('reference','mean'),estimate_mean=('estimate','mean'),
       relative_error_mean=('relative_error','mean'),local_absolute_fraction_mean=('absolute_cell_error_fraction_total','mean')).reset_index()
    synopsis['region']=region;synopsis.to_csv(OUT/(region.replace('/','_')+'_interpretation.csv'),index=False)
    return dict(rows=len(df),complete_factorial_groups=len(counts),repeats_each=100,independent_recomputed_rows=len(rows),
                independent_scope='Repetitions 0,1,49,99; all four mechanisms and three fractions; 10-pixel zero-origin cells; three methods',
                synthetic_bounds=len(bounds),maximum_recomputed_difference=max(abs(r['difference']) for r in rows))

def population_check(g,fields):
    z=np.load(RES/'population/population_models.npz');tab=pd.read_csv(RES/'population/model_epoch_comparison.csv');rows=[]
    names=[n for n in z.files if not n.endswith('_coverage')]
    for dn,dk in [('geographic','domain'),('inherited_strong','screen')]:
        joint=g[dk].copy()
        for n in names:joint&=z[n+'_coverage']>1-1e-5
        total_ref=float(g['population'][joint].sum())
        for n in names:
            for comparison,d in [('joint_valid',joint),('own_valid',g[dk]&(z[n+'_coverage']>1-1e-5))]:
                p=z[n][d];total=float(p.sum());excluded=float(g['population'][g[dk]&~d].sum())
                for sk,s in fields.items():
                    u=float(p@(1-s[d]));r=tab[(tab.domain==dn)&(tab.comparison==comparison)&(tab.model==n)&(tab.support==sk)].iloc[0]
                    assert abs(total-r.population)<1e-5 and abs(u-r.unsupported)<1e-5
                    assert abs(excluded-r.excluded_ghsl_population)<1e-5 and abs(u/total*total_ref-r.unsupported_common_total)<1e-5
                    rows.append(dict(domain=dn,model=n,comparison=comparison,support=sk,population=total,deficit=u,share=u/total,common_total=total_ref,excluded_ghsl=excluded))
    cover=np.load(RES/'population/landcover_fractions.npz');sumcover=sum(cover[k].astype(float) for k in cover.files)
    coverage_error=float(np.max(np.abs(sumcover[g['domain']]-1)));assert coverage_error<2e-5
    land=pd.read_csv(RES/'population/landcover_strata.csv')
    for _,r in land.iterrows():
        weights=g[r.exposure]*g['domain']*cover[str(r.class_code)]
        assert abs(float(weights.sum())-r.total)<max(1e-5,abs(r.total)*1e-12)
        assert abs(float((weights*(1-fields[r.support])).sum())-r.unsupported)<max(1e-5,abs(r.unsupported)*1e-12)
    pd.DataFrame(rows).to_csv(OUT/'population_denominators_recomputed.csv',index=False)
    return dict(model_rows=len(rows),landcover_rows=len(land),maximum_landcover_fraction_closure_error=coverage_error)

def sampling_check():
    p=RES/'sampling';expanded=pd.read_csv(p/'expanded_pair_exposure.csv');nested=pd.read_csv(p/'nested_summary.csv')
    design=pd.read_csv(ROOT/'external/licsar_expanded/experiment_input_plan.csv');plan=design[design.sampling].copy() if 'sampling' in design else design[design.sampling_rank>0].copy()
    records=json.loads((ROOT/'external/licsar_expanded/sampling_manifest.json').read_text());ids={r['pair'] for r in records}
    assert len(ids)==224 and all(r['status']=='ok' for r in records)
    assert set(expanded.pair)==ids
    pool=expanded[(expanded.support=='q30')&(expanded.exposure=='population')].set_index('pair')
    assert pool.groupby(['year','quarter']).size().eq(8).all()
    errors=[]
    for _,r in nested.iterrows():
        v=expanded[(expanded.support==r.support)&(expanded.exposure==r.exposure)&(expanded.sampling_rank<=r.pairs/28)]
        assert len(v)==r.pairs
        # Stored nested support fields are rounded once to float32. Bound the
        # weighted effect by a full float32 ULP of support times TOTAL exposure,
        # not the small deficit denominator (which can approach zero).
        diff=float(v.unsupported.mean()-r.unsupported)
        assert abs(diff)<float(r.total)*2**-24+1e-6,(r.pairs,r.support,r.exposure,diff)
        errors.append(dict(pairs=int(r.pairs),support=r.support,exposure=r.exposure,difference=diff,
                           fraction_total=diff/r.total,tolerance=r.total*2**-24+1e-6))
    draws=pd.read_csv(p/'conditional_sampling_replicates.csv');rng=np.random.default_rng(20260927)
    strata=[q.index.to_numpy() for _,q in pool.groupby(['year','quarter'])]
    for rep in range(100):
        perm=[rng.permutation(s) for s in strata]
        for n in [1,2,4,8]:
            selection=np.concatenate([s[:n] for s in perm]);v=float(pool.loc[selection].unsupported.mean())
            retained=draws[(draws.replicate==rep)&(draws.pairs==28*n)].unsupported.iloc[0];assert abs(v-retained)<1e-7
    legacy=pd.read_csv(RES/'fine/pair_exposure.csv');legacy=legacy[legacy.domain=='geographic']
    networks={}
    for name,ids0 in [('legacy28',legacy.pair.unique()),('expanded224',list(ids)),('catalogue2016_2022',design.pair.unique())]:
        graph=nx.Graph();graph.add_edges_from(p.split('_') for p in ids0);networks[name]=dict(pairs=graph.number_of_edges(),dates=len(graph),components=nx.number_connected_components(graph),cycle_rank=graph.number_of_edges()-len(graph)+nx.number_connected_components(graph))
    saved=json.loads((p/'network.json').read_text())
    for name,n in networks.items():
        for key,v in n.items():assert saved[name][key]==v
    relations=[dict(sample='legacy28',pairs=28,purpose='Historical sample and primary mean support'),
       dict(sample='nested28',pairs=28,purpose='Date/baseline stratified subset of expanded pool'),
       dict(sample='expanded224',pairs=224,purpose='Conditional sampling/temporal sensitivity; disconnected date graph'),
       dict(sample='matched437',pairs=437,purpose='Separate 2019-2020 connected processing ROI; 44 edges held out')]
    pd.DataFrame(relations).to_csv(OUT/'sample_relationships.csv',index=False)
    pd.DataFrame(errors).to_csv(OUT/'nested_float32_rounding.csv',index=False)
    return dict(nested_rows=len(nested),maximum_nested_difference_fraction_total=max(abs(r['fraction_total']) for r in errors),
                nested_tolerance='One float32 support ULP (2^-24) times total exposure plus 1e-6; retained mean fields are rounded to float32',
                conditional_draws_recomputed=len(draws),networks=networks)

def dwr_check():
    z=np.load(RES/'dwr/grid.npz');features=json.loads((ROOT/'external/dwr/roi_features.json').read_text())['features']
    assert len(features)==len({f['attributes']['CODE'] for f in features})==38738
    roi=box(math.radians(-119.9),math.sin(math.radians(36.6)),math.radians(-119.7),math.sin(math.radians(36.85)))
    polys=[]
    for f in features:
        rings=f['geometry']['rings'];assert len(rings)==1
        a=np.array(rings[0],dtype=float);a[:,0]=np.deg2rad(a[:,0]);a[:,1]=np.sin(np.deg2rad(a[:,1]))
        poly=Polygon(a);assert poly.is_valid;polys.append(poly)
    united=union_all(polys).intersection(roi);radius=6371008.8**2
    direct_area=united.area*radius;retained_area=float((z['area_m2']*z['support']).sum())
    assert abs(direct_area-retained_area)<1.
    indices=np.random.default_rng(72026).choice(z['population'].size,128,replace=False);rr,cc=np.unravel_index(indices,z['population'].shape);diff=[]
    x,y=z['x_edges'],z['y_edges']
    for r,c in zip(rr,cc):
        pixel=box(math.radians(x[c]),math.sin(math.radians(y[r])),math.radians(x[c+1]),math.sin(math.radians(y[r+1]))).intersection(roi)
        s=united.intersection(pixel).area/pixel.area if pixel.area else 0;diff.append(s-float(z['support'][r,c]))
    assert max(map(abs,diff))<1e-6
    tab=pd.read_csv(RES/'dwr/summary.csv')
    for _,r in tab.iterrows():
        p=z[r.exposure];assert abs(float(p.sum())-r.total)<1e-5
        assert abs(float((p*(1-z['support'])).sum())-r.unsupported)<1e-5
    return dict(unique_polygons=38738,global_union_area_m2=direct_area,retained_union_area_m2=retained_area,
        area_difference_m2=direct_area-retained_area,independent_pixel_intersections=128,max_support_difference=max(map(abs,diff)),
        method='Global union of original polygon rings in longitude-radians/sine-latitude; independent random pixel intersections',
        tolerance='1 m2 global area and 1e-6 pixel fraction, allowing GEOS floating-point union order')

def main():
    OUT.mkdir(exist_ok=True);report={}
    for region in ['masking','dwr/masking']:
        report[region]=mask_check(region);print(region,report[region],flush=True)
    with np.load(RES/'fine/grid.npz') as z:g={k:z[k] for k in z.files}
    with np.load(RES/'fine/mean_support.npz') as z:fields={k:z[k].astype(float) for k in z.files}
    report['population']=population_check(g,fields);print('population checked',flush=True)
    report['sampling']=sampling_check();print('sampling checked',flush=True)
    report['dwr']=dwr_check();print('DWR checked',flush=True)
    (OUT/'sensitivity_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
