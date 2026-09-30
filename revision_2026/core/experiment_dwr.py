"""E6 published-support polygons: exact union, conservation, scale comparison."""
import json,zipfile,shutil
from collections import defaultdict
from pathlib import Path
import numpy as np,pandas as pd,rasterio
from affine import Affine
from rasterio.warp import transform_bounds,reproject,transform,Resampling
from shapely.geometry import box
from shapely import union_all
from common import ROOT,OUT,dump,edges,pixel_areas,overlaps,sha256
from experiment_population import crop
from fetch_external import fetch
from experiment_allocation import partitions

DEST=OUT/'dwr';DEST.mkdir(exist_ok=True)
SOURCE=ROOT/'external/dwr'
ROI=[-119.9,36.6,-119.7,36.85]
CRS='EPSG:32611'

def grid():
    name='GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0_R6_C7'
    target=SOURCE/(name+'.tif')
    if not target.exists():
        fetch('https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A/GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss/V1-0/tiles/'+name+'.zip',SOURCE/(name+'.zip'))
        with zipfile.ZipFile(SOURCE/(name+'.zip')) as archive:
            member=next(n for n in archive.namelist() if n.endswith('.tif'))
            with archive.open(member) as src,target.open('wb') as out:shutil.copyfileobj(src,out)
    with rasterio.open(SOURCE/'GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0_R6_C7.tif') as ds:
        win=crop(ds,np.array([ROI[0],ROI[2]]),np.array([ROI[1],ROI[3]]));p=ds.read(1,window=win,masked=True)[::-1];tr=ds.window_transform(win)
        if np.ma.getmaskarray(p).any():raise ValueError('Population nodata in second AOI')
        p=p.data.astype('float64');x,y=edges(tr,p.shape)
    with rasterio.open(target) as ds:
        bw=crop(ds,np.array([ROI[0],ROI[2]]),np.array([ROI[1],ROI[3]]));b=ds.read(1,window=bw,masked=True)[::-1]
        if b.shape!=p.shape or np.ma.getmaskarray(b).any() or np.max(np.abs(np.array(ds.window_transform(bw))-np.array(tr)))>1e-9:raise ValueError('Built geometry/nodata')
        b=b.data.astype('float64')
    fullarea=pixel_areas(x,y)
    dx=np.maximum(0,np.minimum(x[1:],ROI[2])-np.maximum(x[:-1],ROI[0]))
    dy=np.maximum(0,np.sin(np.deg2rad(np.minimum(y[1:],ROI[3])))-np.sin(np.deg2rad(np.maximum(y[:-1],ROI[1]))))
    area=(6371008.8**2)*dy[:,None]*np.deg2rad(dx)[None,:]
    fractions=area/fullarea
    return x,y,p*fractions,b*fractions,area,tr

def polygon_support(x,y,area):
    features=json.loads((SOURCE/'roi_features.json').read_text(encoding='utf-8'))['features']
    pieces=defaultdict(list);rawareas=0.
    for f in features:
        rings=f['geometry']['rings']
        if len(rings)!=1:raise ValueError('Nonrectangular observation polygon')
        ring=np.asarray(rings[0]);xs=np.unique(ring[:,0]);ys=np.unique(ring[:,1])
        if len(xs)!=2 or len(ys)!=2:raise ValueError('Non-axis-aligned support: use general polygon intersection')
        left=max(xs[0],ROI[0]);right=min(xs[-1],ROI[2]);bottom=max(ys[0],ROI[1]);top=min(ys[-1],ROI[3])
        if left>=right or bottom>=top:continue
        rawareas+=(np.deg2rad(right-left)*(np.sin(np.deg2rad(top))-np.sin(np.deg2rad(bottom))))*6371008.8**2
        c0=max(0,np.searchsorted(x,left,side='right')-1);c1=min(len(x)-2,np.searchsorted(x,right,side='left'))
        r0=max(0,np.searchsorted(y,bottom,side='right')-1);r1=min(len(y)-2,np.searchsorted(y,top,side='left'))
        for r in range(r0,r1+1):
            lo=max(bottom,y[r]);hi=min(top,y[r+1])
            if lo>=hi:continue
            for c in range(c0,c1+1):
                l=max(left,x[c]);h=min(right,x[c+1])
                if l<h:pieces[(r,c)].append(box(np.deg2rad(l),np.sin(np.deg2rad(lo)),np.deg2rad(h),np.sin(np.deg2rad(hi))))
    covered=np.zeros(area.shape)
    for (r,c),polygons in pieces.items():covered[r,c]=union_all(polygons).area*6371008.8**2
    support=np.divide(covered,area,out=np.zeros_like(area),where=area>0)
    if support.max()>1+1e-6:raise ValueError('Support exceeds domain')
    dump(DEST/'polygon_audit.json',dict(features=len(features),unique_codes=len({f['attributes']['CODE'] for f in features}),raw_polygon_area_m2=rawareas,union_area_m2=float(covered.sum()),overlap_removed_m2=float(rawareas-covered.sum()),domain_area_m2=float(area.sum()),union_intersection='Each source rectangle clipped to AOI and GHSL pixels; union in longitude-radians/sine-latitude coordinates gives exact spherical rectangular area'))
    return np.clip(support,0,1)

def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--frozen-grid',action='store_true',help='Recompute E6 allocation from the released exact footprint integration.')
    if parser.parse_args().frozen_grid:
        with np.load(DEST/'grid.npz') as z:x,y,p,b,a,s=[z[k].copy() for k in ['x_edges','y_edges','population','builtup_m2','area_m2','support']]
        st=Affine(x[1]-x[0],0,x[0],0,-(y[1]-y[0]),y[-1])
    else:
        x,y,p,b,a,st=grid();s=polygon_support(x,y,a)
        np.savez_compressed(DEST/'grid.npz',x_edges=x,y_edges=y,population=p,builtup_m2=b,area_m2=a,support=s)
    left,bottom,right,top=transform_bounds('EPSG:4326',CRS,x[0],y[0],x[-1],y[-1],densify_pts=100)
    left,bottom=np.floor(np.array([left,bottom])/50)*50;right,top=np.ceil(np.array([right,top])/50)*50
    mx=np.arange(left,right+1,50);my=np.arange(bottom,top+1,50);dt=Affine(50,0,left,0,-50,top)
    def project(v):
        out=np.zeros((len(my)-1,len(mx)-1),dtype='float64')
        reproject(v[::-1].copy(),out,src_transform=st,src_crs='EPSG:4326',dst_transform=dt,dst_crs=CRS,resampling=Resampling.sum,src_nodata=None,dst_nodata=0)
        if not np.isclose(out.sum(),v.sum(),rtol=1e-8):raise ValueError('Projection conservation')
        return out[::-1]
    ma,ms=project(a),project(a*s);rows=[];overall=[]
    for en,e in [('population',p),('builtup_m2',b)]:
        ep,eh=project(e),project(e*s);reference=float((e*(1-s)).sum());total=float(e.sum())
        overall.append(dict(exposure=en,total=total,supported=float((e*s).sum()),unsupported=reference,support_share=float((e*s).sum()/total)))
        for size in [250,500,1000,2000]:
            for sx,sy in [(0,0),(size/2,0),(0,size/2),(size/2,size/2)]:
                tx,ty,agg=partitions(mx,my,size,sx,sy);A,S,P,H=[agg(v) for v in [ma,ms,ep,eh]]
                keep=A>1e-6;sc=np.divide(S,A,out=np.zeros_like(A),where=keep);ref=P-H
                xx,yy=np.meshgrid((tx[:-1]+tx[1:])/2,(ty[:-1]+ty[1:])/2);lon,lat=transform(CRS,'EPSG:4326',xx.ravel(),yy.ravel())
                c=np.searchsorted(x,lon,side='right')-1;r=np.searchsorted(y,lat,side='right')-1
                good=(c>=0)&(c<s.shape[1])&(r>=0)&(r<s.shape[0])&(np.array(lon)>=ROI[0])&(np.array(lon)<=ROI[2])&(np.array(lat)>=ROI[1])&(np.array(lat)<=ROI[3]);center=np.zeros(c.shape);center[good]=s[r[good],c[good]];center=center.reshape(xx.shape)
                for method,u in [('area_uniform',P*(1-sc)),('cell_center',P*(1-center))]:
                    err=(u-ref)[keep]
                    rows.append(dict(exposure=en,cell_size_m=size,shift_x_m=sx,shift_y_m=sy,method=method,total=total,reference_unsupported=reference,estimated_unsupported=float(u.sum()),signed_error=float(err.sum()),error_fraction_reference=float(err.sum()/reference),cell_absolute_error_sum=float(np.abs(err).sum()),covariance_integral=float((H-P*sc).sum())))
    pd.DataFrame(rows).to_csv(DEST/'scale_origin_comparison.csv',index=False);pd.DataFrame(overall).to_csv(DEST/'summary.csv',index=False)
    dump(DEST/'method.json',dict(roi=ROI,selection='Fresno fixed geographic rectangle, selected before outcome calculation',support='DWR 2026Q1 published actual observation-cell polygons; nonduplicated union',epoch='GHSL 2020 compared with published 2015-2026 observation availability; explicit epoch mismatch',crs=CRS,metric_integration_m=50,meaning='Published observation coverage, not a coherence proxy or interpolated velocity raster',no_unobserved_velocity_inference=True,limitation='No per-cell uncertainty or velocity was obtained from this polygon layer. It validates transferable support geometry/accounting, not independent physical accuracy.'))
    print('E6 DWR GEOMETRY COMPLETE',overall,flush=True)

if __name__=='__main__':main()
