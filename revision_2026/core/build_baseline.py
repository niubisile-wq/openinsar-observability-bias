"""E0 common fine-grid baseline with conservative spherical intersections."""
from pathlib import Path
import json,time
import numpy as np,pandas as pd,rasterio
from rasterio.merge import merge
from affine import Affine
from common import ROOT,CFG,PROJECT,BACKUP,OUT,dump,sha256,edges,pixel_areas,sample_array,exact_joint_support

FINE=OUT/'fine';FINE.mkdir(exist_ok=True)
PAIR=FINE/'pairs';PAIR.mkdir(exist_ok=True)
R=BACKUP/'revision_scientific_reports_v1'

def make_grid():
    cells=pd.read_csv(PROJECT/'analysis/reconstructed_cells.csv')
    vlm=PROJECT/'archive/nature/03_exposure_closure/chao_phraya/gridVLM/chaoPhraya_vlm.tif'
    with rasterio.open(vlm) as src:
        vt=src.transform;vv=src.read(1)
    west=float((cells.lon_center-vt.a/2).min());east=float((cells.lon_center+vt.a/2).max())
    south=float((cells.lat_center+vt.e/2).min());north=float((cells.lat_center-vt.e/2).max())
    ref=PROJECT/'raw/ghsl_3ss/GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0_R8_C29.tif'
    with rasterio.open(ref) as src: tr=src.transform
    c0=int(np.floor((west-tr.c)/tr.a));c1=int(np.ceil((east-tr.c)/tr.a))
    r0=int(np.floor((north-tr.f)/tr.e));r1=int(np.ceil((south-tr.f)/tr.e))
    ft=tr*Affine.translation(c0,r0);shape=(r1-r0,c1-c0)
    tx,ty=edges(ft,shape)
    masks=np.zeros(vv.shape,dtype='uint8');masks[cells.row,cells.col]=1
    strong=np.zeros(vv.shape,dtype='uint8');strong[cells.row,cells.col]=(cells.vlm_mm_yr<=-5).astype('uint8')
    x=(tx[:-1]+tx[1:])/2;y=(ty[:-1]+ty[1:])/2
    domain=sample_array(masks,vt,x,y).astype(bool);screen=sample_array(strong,vt,x,y).astype(bool)
    exposures={}
    manifest=[dict(path=str(vlm),sha256=sha256(vlm))]
    for product,key in [('GHS_POP','population'),('GHS_BUILT_S','builtup_m2')]:
        paths=[PROJECT/f'raw/ghsl_3ss/{product}_E2020_GLOBE_R2023A_4326_3ss_V1_0_R8_C{c}.tif' for c in [28,29]]
        ds=[rasterio.open(p) for p in paths]
        a,mt=merge(ds,bounds=(tx[0],ty[0],tx[-1],ty[-1]),res=(tr.a,-tr.e))
        for d in ds:d.close()
        a=a[0,::-1].astype('float64')
        if a.shape!=shape or not np.allclose(list(mt)[:6],list(ft)[:6],rtol=0,atol=1e-10):raise ValueError('Crop alignment')
        valid=np.isfinite(a)&(a>=0)
        if key=='builtup_m2':valid&=a<65535
        if not valid[domain].all():raise ValueError('Exposure no-data inside analysis domain')
        exposures[key]=np.where(valid,a,0.)
        manifest += [dict(path=str(p),sha256=sha256(p)) for p in paths]
    area=pixel_areas(tx,ty)
    np.savez_compressed(FINE/'grid.npz',x_edges=tx,y_edges=ty,domain=domain,screen=screen,area_m2=area,**exposures)
    dump(FINE/'grid.json',dict(shape=list(shape),crs='EPSG:4326',transform=list(ft)[:6],row_order='south-to-north',domain_definition='GHSL 3ss pixels whose centers fall in the frozen 18077-cell VLM footprint; does not condition on velocity scale',primary_cells=int(domain.sum()),conditional_screen_cells=int(screen.sum()),inputs=manifest))
    return tx,ty,domain,screen,exposures,area

def main():
    tx,ty,domain,screen,exposures,area=make_grid()
    cc_paths=sorted((PROJECT/'raw/licsar').glob('*.geo.cc.tif'))
    if len(cc_paths)!=28:raise ValueError('Expected 28 legacy pairs')
    totals={k:np.zeros(domain.shape,dtype='float64') for k in ['unw','q20','q30','q40']}
    manifest=[];pair_rows=[]
    for n,p in enumerate(cc_paths):
        pair=p.name.split('.')[0];u=R/f'inputs_needed/licsar_chao_28_unw/{pair}.geo.unw.tif'
        dest=PAIR/(pair+'.npz');started=time.time()
        if dest.exists():
            with np.load(dest) as data: support={k:data[k] for k in totals}
        else:
            with rasterio.open(p) as cd,rasterio.open(u) as ud:
                if str(cd.crs)!='EPSG:4326' or str(ud.crs)!='EPSG:4326':raise ValueError('CRS')
                cc=cd.read(1);unw=ud.read(1)
                uv=(ud.read_masks(1)>0)&np.isfinite(unw)&(unw!=0)
                cv=(cd.read_masks(1)>0)&np.isfinite(cc)
                support=exact_joint_support(cc,cv,cd.transform,uv,ud.transform,tx,ty,CFG['coherence_thresholds'])
            np.savez_compressed(dest,**support)
        if any(np.any(support[a]+1e-6<support[b]) for a,b in [('unw','q20'),('q20','q30'),('q30','q40')]):raise ValueError('Nonmonotone masks')
        for k,f in support.items():
            totals[k]+=f
            for dn,d in [('geographic',domain),('inherited_strong',screen)]:
                for en,e in {**exposures,'area_m2':area}.items():
                    weight=e*d
                    pair_rows.append(dict(pair=pair,support=k,domain=dn,exposure=en,total=float(weight.sum()),supported=float((weight*f).sum()),unsupported=float((weight*(1-f)).sum())))
        manifest.append(dict(pair=pair,cc=str(p),unw=str(u),cc_sha256=sha256(p),unw_sha256=sha256(u),fine_file=dest.name,fine_sha256=sha256(dest)))
        print(f'{n+1}/28 {pair} {time.time()-started:.1f}s',flush=True)
    fields={k:(v/len(cc_paths)).astype('float32') for k,v in totals.items()}
    np.savez_compressed(FINE/'mean_support.npz',**fields)
    pd.DataFrame(pair_rows).to_csv(FINE/'pair_exposure.csv',index=False)
    dump(FINE/'pair_manifest.json',manifest)
    rows=[];curve=[]
    for k,f in fields.items():
        for dn,d in [('geographic',domain),('inherited_strong',screen)]:
            for en,e in {**exposures,'area_m2':area}.items():
                w=e*d;total=float(w.sum());sup=float((w*f).sum());uns=float((w*(1-f.astype('float64'))).sum())
                if not np.isclose(sup+uns,total,rtol=1e-10):raise ValueError('Conservation failure')
                rows.append(dict(support=k,domain=dn,exposure=en,total=total,supported=sup,unsupported=uns,weighted_support=sup/total))
                for t in np.linspace(0,1,101):curve.append(dict(support=k,domain=dn,exposure=en,threshold=t,below_threshold=float(w[f<t].sum()),share_below=float(w[f<t].sum()/total)))
    pd.DataFrame(rows).to_csv(FINE/'allocation_summary.csv',index=False)
    pd.DataFrame(curve).to_csv(FINE/'support_curves.csv',index=False)
    dump(FINE/'method.json',dict(method='Exact spherical rectangle intersections on union of cc and unw source boundaries, aggregated to GHSL3ss pixels; equal-weight mean over the 28 pairs',unwrapped_rule='GDAL mask AND finite AND nonzero, product-specific',support_interpretation='Mean area-pair quality support, not final-product validity or velocity accuracy',epoch=2020))
    print('BASELINE COMPLETE',flush=True)

if __name__=='__main__':main()
