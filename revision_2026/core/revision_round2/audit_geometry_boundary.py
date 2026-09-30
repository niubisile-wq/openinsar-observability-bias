"""Independent native-pixel overlap and metric allocation/boundary checks.

Does not import production common/experiment_allocation functions. Source raster
rectangles are intersected directly for deterministic sampled pixels. Metric
moments are aggregated with independently constructed source-fraction weights.
"""
from pathlib import Path
import json, math
import numpy as np,pandas as pd,rasterio,tifffile
from rasterio.features import shapes
from rasterio.warp import transform,reproject,Resampling
from affine import Affine
from scipy.sparse import csr_matrix
from shapely import box,covers,segmentize,prepare
from shapely.geometry import shape
from shapely.ops import unary_union,transform as transform_geom
from pyproj import Transformer

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;PROJECT=ROOT.parent;MASTER=PROJECT.parent.parent
DEST=HERE/'evidence';DEST.mkdir(exist_ok=True)
BACKUP=next(p.parent for p in (MASTER/'02_历史备份').glob('*/nature/revision_scientific_reports_v1')
            if (p/'inputs_needed/licsar_chao_28_unw').is_dir())

def area(l,b,r,t):
    return max(0,math.radians(r-l))*max(0,math.sin(math.radians(t))-math.sin(math.radians(b)))

def rects(ds,a,valid,l,b,r,t):
    tr=ds.transform
    c0=max(0,math.floor((l-tr.c)/tr.a));c1=min(ds.width,math.ceil((r-tr.c)/tr.a))
    r0=max(0,math.floor((t-tr.f)/tr.e));r1=min(ds.height,math.ceil((b-tr.f)/tr.e))
    result=[]
    for iy in range(r0,r1):
        for ix in range(c0,c1):
            left=max(l,tr.c+ix*tr.a);right=min(r,tr.c+(ix+1)*tr.a)
            bottom=max(b,tr.f+(iy+1)*tr.e);top=min(t,tr.f+iy*tr.e)
            if right>left and top>bottom:
                result.append((left,bottom,right,top,bool(valid[iy,ix]),float(a[iy,ix])))
    return result

def native_check(g):
    rng=np.random.default_rng(27192026)
    ids=np.flatnonzero(g['domain']);chosen=rng.choice(ids,64,replace=False)
    # Include corners of the finite footprint and low/high support examples.
    mean=np.load(ROOT/'results/fine/mean_support.npz')['q30']
    sorted_ids=ids[np.argsort(mean.ravel()[ids])]
    chosen=np.unique(np.r_[chosen,sorted_ids[:8],sorted_ids[-8:]])
    yy,xx=np.unravel_index(chosen,g['domain'].shape);rows=[];geo=[]
    x,y=g['x_edges'],g['y_edges']
    for cp in sorted((PROJECT/'raw/licsar').glob('*.geo.cc.tif')):
        pair=cp.name.split('.')[0];up=BACKUP/f'revision_scientific_reports_v1/inputs_needed/licsar_chao_28_unw/{pair}.geo.unw.tif'
        with np.load(ROOT/f'results/fine/pairs/{pair}.npz') as z:
            cache={k:z[k] for k in z.files}
        with rasterio.open(cp) as cd,rasterio.open(up) as ud:
            c=cd.read(1);u=ud.read(1);cv=cd.read_masks(1)>0;uv=(ud.read_masks(1)>0)&np.isfinite(u)&(u!=0)
            for path,ds in [(cp,cd),(up,ud)]:
                with tifffile.TiffFile(path) as tf:
                    tags=tf.pages[0].tags;tie=tags.get('ModelTiepointTag');scale=tags.get('ModelPixelScaleTag')
                    if tie and scale:
                        ti=tie.value;sc=scale.value;point=ds.tags().get('AREA_OR_POINT')=='Point'
                        left=ti[3]-ti[0]*sc[0]-(sc[0]/2 if point else 0)
                        top=ti[4]+ti[1]*sc[1]+(sc[1]/2 if point else 0)
                        diff=max(abs(left-ds.transform.c),abs(top-ds.transform.f))
                        assert diff<1e-9
                        geo.append(dict(pair=pair,kind='cc' if path==cp else 'unw',pixel_is_point=point,affine_max_difference_deg=diff))
            for iy,ix in zip(yy,xx):
                l,b,r,t=float(x[ix]),float(y[iy]),float(x[ix+1]),float(y[iy+1]);den=area(l,b,r,t)
                cr=rects(cd,c,cv,l,b,r,t);ur=rects(ud,u,uv,l,b,r,t)
                values={'unw':sum(area(*q[:4]) for q in ur if q[4])/den}
                for threshold in [.2,.3,.4]:
                    num=0.
                    for ca in cr:
                        if not ca[4] or ca[5]/255 < threshold:continue
                        for ua in ur:
                            if not ua[4]:continue
                            num+=area(max(ca[0],ua[0]),max(ca[1],ua[1]),min(ca[2],ua[2]),min(ca[3],ua[3]))
                    values[f'q{int(threshold*100)}']=num/den
                for k,v in values.items():
                    diff=v-float(cache[k][iy,ix]);assert abs(diff)<1e-6,(pair,k,iy,ix,diff)
                    rows.append(dict(pair=pair,row=int(iy),col=int(ix),support=k,direct_fraction=v,retained_fraction=float(cache[k][iy,ix]),difference=diff))
        print('native overlap checked',pair,flush=True)
    pd.DataFrame(rows).to_csv(DEST/'native_overlap_independent.csv',index=False)
    pd.DataFrame(geo).to_csv(DEST/'geotiff_registration.csv',index=False)
    return dict(pairs=28,sampled_pixels=len(chosen),comparisons=len(rows),max_abs_fraction_difference=max(abs(r['difference']) for r in rows),
        method='Direct source-rectangle nested intersections; deterministic random plus support-extreme pixels; all 28 pairs; independent of production remapper',
        tolerance='1e-6 absolute fraction, allowing retained float32 rounding; affine tolerance 1e-9 degrees')

def weight_matrix(src,dst):
    ri=[];ci=[];vv=[]
    for j,(l,r) in enumerate(zip(src[:-1],src[1:])):
        first=max(0,int(np.searchsorted(dst,l,side='right'))-1)
        for i in range(first,len(dst)-1):
            if dst[i]>=r:break
            length=min(r,dst[i+1])-max(l,dst[i])
            if length>0:ri.append(i);ci.append(j);vv.append(length/(r-l))
    return csr_matrix((vv,(ri,ci)),shape=(len(dst)-1,len(src)-1))

def allocation_check(g):
    path=ROOT/'results/allocation';m=np.load(path/'metric_geometry.npz');mx,my=m['x_edges'],m['y_edges']
    f=np.load(ROOT/'results/fine/mean_support.npz')['q30'].astype(float)
    fields={k:np.load(path/'metric50'/v,mmap_mode='r') for k,v in {
        'A':'geographic_area.npy','S':'geographic_q30_supported_area.npy',
        'P':'geographic_population.npy','H':'geographic_q30_population_cross.npy'}.items()}
    x,y=g['x_edges'],g['y_edges'];st=Affine(x[1]-x[0],0,x[0],0,-(y[1]-y[0]),y[-1])
    dt=Affine(50,0,mx[0],0,-50,my[-1]);fresh=[]
    a=g['area_m2']*g['domain'];p=g['population']*g['domain']
    for k,v in [('A',a),('S',a*f),('P',p),('H',p*f)]:
        dst=np.zeros(fields[k].shape,dtype=float)
        reproject(v[::-1].copy(),dst,src_transform=st,src_crs='EPSG:4326',dst_transform=dt,dst_crs='EPSG:32647',resampling=Resampling.sum,src_nodata=None,dst_nodata=0)
        diff=float(np.max(np.abs(dst[::-1]-fields[k])));assert diff<1e-6
        fresh.append(dict(moment=k,max_pixel_difference=diff,source_total=float(v.sum()),metric_total=float(dst.sum())))
    # Pixel-defined geographic boundary, densified to <= one native pixel before projection.
    polys=[shape(geom) for geom,val in shapes(g['domain'][::-1].astype('uint8'),transform=st) if val==1]
    geographic=segmentize(unary_union(polys),min(x[1]-x[0],y[1]-y[0]))
    proj=Transformer.from_crs(4326,32647,always_xy=True)
    domain=transform_geom(proj.transform,geographic)
    prepare(domain)
    ref_table=pd.read_csv(path/'scale_origin_comparison.csv');rows=[];checks=[]
    for size in [250,500,1000,2000]:
        for sx,sy in [(0,0),(size/2,0),(0,size/2),(size/2,size/2)]:
            tx=np.arange(math.floor((mx[0]-sx)/size)*size+sx,math.ceil((mx[-1]-sx)/size)*size+sx+size*.1,size)
            ty=np.arange(math.floor((my[0]-sy)/size)*size+sy,math.ceil((my[-1]-sy)/size)*size+sy+size*.1,size)
            wx,wy=weight_matrix(mx,tx),weight_matrix(my,ty)
            z={k:(wy@v)@wx.T for k,v in fields.items()};A,S,P,H=(z[k] for k in ['A','S','P','H'])
            keep=A>1e-6;sc=np.divide(S,A,out=np.zeros_like(S),where=keep);reference=P-H;uniform=P*(1-sc)
            xx,yy=np.meshgrid((tx[:-1]+tx[1:])/2,(ty[:-1]+ty[1:])/2)
            lon,lat=transform('EPSG:32647','EPSG:4326',xx.ravel(),yy.ravel())
            ix=np.searchsorted(x,lon,side='right')-1;iy=np.searchsorted(y,lat,side='right')-1
            valid=(ix>=0)&(ix<f.shape[1])&(iy>=0)&(iy<f.shape[0]);ids=np.flatnonzero(valid);ids=ids[g['domain'][iy[ids],ix[ids]]]
            center_support=np.zeros(ix.size);center_support[ids]=f[iy[ids],ix[ids]]
            center=P*(1-center_support.reshape(P.shape))
            rectangles=box(xx-size/2,yy-size/2,xx+size/2,yy+size/2)
            interior=keep&covers(domain,rectangles);boundary=keep&~interior
            center_out=keep&~np.isin(np.arange(ix.size),ids).reshape(P.shape)
            for method,estimate in [('area_uniform',uniform),('cell_center',center)]:
                old=ref_table[(ref_table.domain=='geographic')&(ref_table.support=='q30')&(ref_table.exposure=='population')&
                    (ref_table.cell_size_m==size)&(ref_table.shift_x_m==sx)&(ref_table.shift_y_m==sy)&(ref_table.method==method)].iloc[0]
                diff=float(estimate.sum()-old.estimated_unsupported);assert abs(diff)<1e-4
                checks.append(dict(size=size,sx=sx,sy=sy,method=method,recalculated=float(estimate.sum()),retained=float(old.estimated_unsupported),difference=diff))
                for name,sel in [('all',keep),('interior_complete',interior),('boundary_partial',boundary),('center_outside',center_out)]:
                    total=float(P[sel].sum());ref=float(reference[sel].sum());u=float(estimate[sel].sum());err=estimate[sel]-reference[sel]
                    rows.append(dict(size_m=size,shift_x=sx,shift_y=sy,method=method,subset=name,cells=int(sel.sum()),
                        population=total,reference_deficit=ref,estimate=u,signed_error=u-ref,
                        error_pct_reference=100*(u-ref)/ref if ref else None,error_pct_population=100*(u-ref)/total if total else None,
                        absolute_cell_error_sum=float(np.abs(err).sum()),positive_cells=int((err>1e-6).sum()),negative_cells=int((err< -1e-6).sum())))
            print('allocation and boundary',size,sx,sy,flush=True)
    pd.DataFrame(rows).to_csv(DEST/'boundary_allocation.csv',index=False)
    pd.DataFrame(checks).to_csv(DEST/'allocation_independent.csv',index=False)
    pd.DataFrame(fresh).to_csv(DEST/'fresh_metric_projection.csv',index=False)
    return dict(comparisons=len(checks),maximum_total_difference=max(abs(c['difference']) for c in checks),
        boundary_definition='Coarse UTM rectangles fully covered by the union of native GHSL domain pixels; geographic boundaries densified at native pixel spacing before projection',
        fresh_projection=fresh)

def main():
    with np.load(ROOT/'results/fine/grid.npz') as z:g={k:z[k] for k in z.files}
    import sys
    if '--boundary-only' in sys.argv:
        previous=pd.read_csv(DEST/'native_overlap_independent.csv')
        native=dict(pairs=int(previous.pair.nunique()),sampled_pixels=int(previous[['row','col']].drop_duplicates().shape[0]),
            comparisons=len(previous),max_abs_fraction_difference=float(previous.difference.abs().max()),
            method='Retained completed direct-overlap check from this round; continuing boundary calculation after performance optimization')
    else:
        native=native_check(g)
    report={'native':native,'allocation':allocation_check(g)}
    (DEST/'geometry_boundary_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
