"""E5 alternate population model, fixed-resolution epoch and cover strata."""
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import from_bounds,Window
from common import OUT,ROOT,PROJECT,dump,sha256,edges,separable_resample,overlaps

DEST=OUT/'population';DEST.mkdir(exist_ok=True)
CLASSES={10:'Tree cover',20:'Shrubland',30:'Grassland',40:'Cropland',50:'Built-up',60:'Bare/sparse vegetation',70:'Snow/ice',80:'Permanent water',90:'Herbaceous wetland',95:'Mangroves',100:'Moss/lichen'}

def crop(src,tx,ty):
    w=from_bounds(tx[0],ty[0],tx[-1],ty[-1],src.transform)
    c0=max(0,int(np.floor(w.col_off)));r0=max(0,int(np.floor(w.row_off)))
    c1=min(src.width,int(np.ceil(w.col_off+w.width)));r1=min(src.height,int(np.ceil(w.row_off+w.height)))
    return Window(c0,r0,c1-c0,r1-r0)

def population(path,tx,ty):
    with rasterio.open(path) as ds:
        w=crop(ds,tx,ty);v=ds.read(1,window=w,masked=True)
        sx,sy=edges(ds.window_transform(w),v.shape)
        valid=~np.ma.getmaskarray(v)[::-1]&np.isfinite(v.data[::-1])&(v.data[::-1]>=0)
        a=np.where(valid,v.data[::-1],0.).astype('float64')
    aligned=separable_resample(a,sx,sy,tx,ty,normalization='source')
    coverage=separable_resample(valid.astype('float32'),sx,sy,tx,ty)
    return aligned,coverage

def cover(tx,ty):
    path=ROOT/'external/worldcover2021_N12E099.tif';cache=DEST/'landcover_fractions.npz'
    if cache.exists():
        with np.load(cache) as z:return {int(k):z[k] for k in z.files}
    result={k:np.zeros((len(ty)-1,len(tx)-1),dtype='float32') for k in CLASSES}
    validtotal=np.zeros_like(result[10])
    with rasterio.open(path) as ds:
        w=crop(ds,tx,ty)
        # Stream native 10 m rows; exact spherical area overlaps, no modal
        # assignment of an entire population pixel to a minority land class.
        for r0 in range(int(w.row_off),int(w.row_off+w.height),512):
            n=min(512,int(w.row_off+w.height)-r0)
            stripe=Window(w.col_off,r0,w.width,n);a=ds.read(1,window=stripe)[::-1]
            sx,sy=edges(ds.window_transform(stripe),a.shape)
            wx=overlaps(np.deg2rad(sx),np.deg2rad(tx),'target')
            wy=overlaps(np.sin(np.deg2rad(sy)),np.sin(np.deg2rad(ty)),'target')
            for k in np.unique(a):
                if int(k) not in result:continue
                result[int(k)]+=(wx@(wy@(a==k).astype('float32')).T).T.astype('float32')
            print('LANDCOVER rows',r0-int(w.row_off)+n,'/',int(w.height),flush=True)
    np.savez_compressed(cache,**{str(k):v for k,v in result.items()})
    return result

def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--frozen-models',action='store_true',help='Recompute E5 from the released aligned population/land-cover arrays.')
    frozen=parser.parse_args().frozen_models
    with np.load(OUT/'fine/grid.npz') as z:g={k:z[k] for k in z.files}
    with np.load(OUT/'fine/mean_support.npz') as z:support={k:z[k].astype('float64') for k in z.files}
    tx,ty=g['x_edges'],g['y_edges'];rows=[];curves=[];manifest=[]
    models={'GHSL_2020_3ss':(g['population'],np.ones(g['domain'].shape))}
    sources={'WorldPop_2020_3ss':ROOT/'external/worldpop_thailand_2020_UNadj.tif','GHSL_2015_30ss':PROJECT/'analysis/GHS_POP_2015_study_crop.tif','GHSL_2020_30ss':PROJECT/'analysis/GHS_POP_2020_study_crop.tif'}
    if frozen:
        with np.load(DEST/'population_models.npz') as z:
            for name in sources:models[name]=(z[name].copy(),z[name+'_coverage'].copy())
        manifest.append(dict(input='Frozen aligned population_models.npz',sha256=sha256(DEST/'population_models.npz')))
    else:
        for name,path in sources.items():
            models[name]=population(path,tx,ty);manifest.append(dict(model=name,path=str(path),sha256=sha256(path)))
        np.savez_compressed(DEST/'population_models.npz',**{k:v[0] for k,v in models.items()},**{k+'_coverage':v[1] for k,v in models.items()})
    for domain,dk in [('geographic','domain'),('inherited_strong','screen')]:
        # Compare population models on their jointly valid support; report the
        # excluded GHSL population rather than filling WorldPop nodata as zero.
        joint=g[dk].copy()
        for v,c in models.values():joint&=c>1-1e-5
        reference_total=float(g['population'][joint].sum())
        for name,(p,c) in models.items():
            for comparison,d in [('own_valid',g[dk]&(c>1-1e-5)),('joint_valid',joint)]:
                weight=p*d;total=float(weight.sum())
                for sk,s in support.items():
                    u=float((weight*(1-s)).sum())
                    rows.append(dict(domain=domain,comparison=comparison,model=name,support=sk,cells=int(d.sum()),population=total,unsupported=u,unsupported_share=u/total,common_total=reference_total,unsupported_common_total=u/total*reference_total,excluded_ghsl_population=float(g['population'][g[dk]&~d].sum())))
                    if comparison=='joint_valid':
                        for t in np.linspace(0,1,101):curves.append(dict(domain=domain,model=name,support=sk,threshold=t,below_share=float(weight[s<t].sum()/total)))
    pd.DataFrame(rows).to_csv(DEST/'model_epoch_comparison.csv',index=False)
    pd.DataFrame(curves).to_csv(DEST/'model_curves.csv',index=False)
    fractions=cover(tx,ty);totalcover=sum(fractions.values())
    if np.max(np.abs(totalcover[g['domain']]-1))>2e-5:raise ValueError('Incomplete land-cover classification')
    landrows=[]
    for cls,f in fractions.items():
        for exposure in ['population','builtup_m2','area_m2']:
            w=g[exposure]*g['domain']*f;total=float(w.sum())
            if total<=0:continue
            for sk,s in support.items():landrows.append(dict(class_code=cls,class_name=CLASSES[cls],exposure=exposure,support=sk,total=total,unsupported=float((w*(1-s)).sum()),weighted_support=float((w*s).sum()/total)))
    pd.DataFrame(landrows).to_csv(DEST/'landcover_strata.csv',index=False)
    cover_path=DEST/'landcover_fractions.npz' if frozen else ROOT/'external/worldcover2021_N12E099.tif'
    manifest.append(dict(model='ESA_WorldCover_2021_v200',path=str(cover_path),sha256=sha256(cover_path)))
    dump(DEST/'method.json',dict(inputs=manifest,population_transfer='Exact spherical source-area fractions onto GHSL fine grid; no interpolation of counts',model_comparison='GHSL and WorldPop 2020 at native 3ss, including native totals and common-total normalization',epoch_comparison='GHSL 2015 versus 2020 both at 30ss; these are modelled epoch differences, not measured population migration',joint_coverage='All four models must cover each included fine cell to within 1e-5',landcover='ESA WorldCover 2021 v200, exact native-class area fractions within GHSL 3ss pixels',landcover_limitation='Population and built area are assumed uniform within each GHSL pixel for splitting among cover classes; cover stratification does not resolve building occupancy',dependencies='GHSL population and built layers share modelling inputs; WorldPop is an alternative model, not census ground truth'))
    print('E5 COMPLETE',flush=True)

if __name__=='__main__':main()
