"""Prepare two outcome-independent regions using frozen provider rasters."""
from pathlib import Path
import json
import hashlib
import numpy as np
import rasterio
from rasterio.windows import from_bounds,Window
from rasterio.warp import reproject,Resampling
from affine import Affine
import shapely
from pyproj import Transformer
ROOT=Path(__file__).resolve().parent.parent
STUDY=ROOT.parent/'01_当前返修工作区/InSAR_revision_20260926/strengthening'
SRC=ROOT/'inputs/aria_california'
OUT=ROOT/'inputs/new_regions'
OUT.mkdir(parents=True,exist_ok=True)
REGIONS=[('Los_Angeles','A064',[-118.4,33.85,-118.1,34.1],'20160214','20220513','P597',
          STUDY/'external/dwr/GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0_R6_C7.tif'),
         ('Davis_Sacramento','A137',[-121.9,38.4,-121.6,38.7],'20170303','20230829','P266',
          ROOT/'inputs/population/GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0_R6_C6.tif')]

def crop(path,bbox):
    with rasterio.open(path) as ds:
        requested=from_bounds(*bbox,transform=ds.transform)
        c0=max(0,int(np.floor(requested.col_off))-2);r0=max(0,int(np.floor(requested.row_off))-2)
        c1=min(ds.width,int(np.ceil(requested.col_off+requested.width))+2);r1=min(ds.height,int(np.ceil(requested.row_off+requested.height))+2)
        win=Window(c0,r0,c1-c0,r1-r0)
        arr=ds.read(1,window=win,masked=True).astype(float)
        return np.where(np.ma.getmaskarray(arr),np.nan,arr.data),ds.window_transform(win)

def overlaps(target_edges,source_edges):
    # Geographic rectangle areas use longitude span and sine(latitude) span.
    return np.maximum(0,np.minimum(target_edges[1:,None],source_edges[None,1:])-np.maximum(target_edges[:-1,None],source_edges[None,:-1]))

def main():
    projection=Transformer.from_crs('EPSG:4326','EPSG:3310',always_xy=True)
    for region,track,bbox,start,end,ref,pop_path in REGIONS:
        p,pt=crop(pop_path,bbox);p=np.nan_to_num(p,nan=0);p[p<0]=0
        vel,vt=crop(SRC/f'{track}_velocity_ref{ref}.tif',bbox)
        coh,ct=crop(SRC/f'{track}_avgSpatialCoh.tif',bbox)
        assert vel.shape==coh.shape and vt==ct
        assert abs(vt.a-0.000833334)<1e-10
        h,w=p.shape;sh,sw=vel.shape
        px=pt.c+np.arange(w+1)*pt.a;py=pt.f+np.arange(h+1)*pt.e
        vx=vt.c+np.arange(sw+1)*vt.a;vy=vt.f+np.arange(sh+1)*vt.e
        clipx0=np.maximum(px[:-1],bbox[0]);clipx1=np.minimum(px[1:],bbox[2])
        clipy0=np.maximum(py[1:],bbox[1]);clipy1=np.minimum(py[:-1],bbox[3])
        xwidth=np.maximum(0,clipx1-clipx0)
        ywidth=np.maximum(0,np.sin(np.deg2rad(clipy1))-np.sin(np.deg2rad(clipy0)))
        xwhole=np.diff(px);ywhole=np.sin(np.deg2rad(py[:-1]))-np.sin(np.deg2rad(py[1:]))
        fraction=np.outer(ywidth/ywhole,xwidth/xwhole)
        population=p*fraction
        dxover=np.maximum(0,np.minimum(clipx1[:,None],vx[None,1:])-np.maximum(clipx0[:,None],vx[None,:-1]))
        dyover=np.maximum(0,np.sin(np.deg2rad(np.minimum(clipy1[:,None],vy[None,:-1])))-np.sin(np.deg2rad(np.maximum(clipy0[:,None],vy[None,1:]))))
        supports={}
        finite=np.isfinite(vel)&(vel!=0)&np.isfinite(coh)
        for q in [.2,.3,.4]:
            observed=finite&(coh>=q)
            integral=(dyover@observed.astype(float))@dxover.T
            supports['support_q'+str(int(q*100))]=np.divide(integral,np.outer(ywidth,xwidth),out=np.zeros(p.shape),where=np.outer(ywidth,xwidth)>0)
            assert supports['support_q'+str(int(q*100))].max()<=1+1e-9
        yy,xx=np.indices(p.shape);keep=fraction.ravel()>0
        x0=clipx0[xx.ravel()[keep]];x1=clipx1[xx.ravel()[keep]]
        y0=clipy0[yy.ravel()[keep]];y1=clipy1[yy.ravel()[keep]]
        cells=shapely.box(x0,y0,x1,y1)
        cells=shapely.transform(cells,lambda xy:np.column_stack(projection.transform(xy[:,0],xy[:,1])))
        area=shapely.area(cells)
        arrays={'population':population.ravel()[keep],'area_m2':area,'geometry_wkb':shapely.to_wkb(cells,hex=True).astype(str),
                'x_center_m':shapely.get_x(shapely.centroid(cells)), 'y_center_m':shapely.get_y(shapely.centroid(cells))}
        arrays.update({k:v.ravel()[keep] for k,v in supports.items()})
        np.savez_compressed(OUT/(region+'.npz'),**arrays)
        protocol={'name':region,'track':track,'bbox_wgs84':bbox,'observed_start':start,'observed_end':end,
                  'date_provenance':'Archived completed notebook velocity figure date range, verified against Sangha et al. 2026 Table 1; config run_notes can be stale for A137.',
                  'LOS_sign':'Positive toward satellite','velocity_unit':'m/year in archived source; no physical velocity classes are inferred here',
                  'reference_station':ref,'population_epoch':2020,'population_release':'GHS-POP R2023A 3 arc-second cell counts',
                  'population_units':float(population.sum()),'area_m2':float(area.sum()),'native_cells':len(area),
                  'mean_coherence_units':'unitless, archived mean spatial coherence',
                  'support_rule':'Finite, nonzero archived relative LOS velocity and finite mean spatial coherence >= threshold; zero treated as excluded according to provider notebook.',
                  'geometric_method':'Exact geographic rectangle overlap with longitude/sine-latitude area weights, clipped to fixed ROI; GHSL density homogeneous within each source cell. Project clipped source cells to EPSG:3310 for reporting-unit overlap.',
                  'checks':{'region_chosen_before_outcomes':True,'velocity_coherence_grids_identical':True,'support_range_0_1':True,'no_center_only_support_assignment':True},
                  'limits':['Not a newly processed SAR inversion or new independent population truth.',
                            'Mean spatial coherence and finite rate masks give a product-support criterion, not external physical accuracy.',
                            'Regional population-weighted comparisons are conditional on GHS-POP allocations and their 2020 epoch.',
                            'Fractional support is integrated at GHSL scale then held constant inside each clipped GHSL cell when forming larger reporting units.'],
                  'source_sha256':{str(path.relative_to(ROOT.parent)).replace('\\','/'):hashlib.sha256(path.read_bytes()).hexdigest() for path in
                                   [pop_path,SRC/f'{track}_velocity_ref{ref}.tif',SRC/f'{track}_avgSpatialCoh.tif']}}
        (OUT/(region+'_protocol.json')).write_bytes((json.dumps(protocol,indent=2)+'\n').encode())
        print(region,len(area),'population',population.sum(),'support ranges',[(k,float(v.min()),float(v.max())) for k,v in supports.items()],flush=True)

if __name__=='__main__':main()
