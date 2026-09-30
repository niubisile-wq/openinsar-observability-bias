import os
from pathlib import Path
import os,json,hashlib
import numpy as np,pandas as pd
import rasterio
from rasterio.features import rasterize
from affine import Affine
BASE=Path(os.environ['BASE']); PROJ=Path(os.environ['INSAR_STRENGTHENING_ROOT']); BLOCK_DIR=BASE/'census_blocks'/'fresno_roi'; OUT=BASE/'census_block_analysis';OUT.mkdir(exist_ok=True)
ROI=(-119.9,36.6,-119.7,36.85);R=6371008.8
def load_blocks():
 fs=[]
 for f in sorted(BLOCK_DIR.glob('*.geojson')):fs.extend(json.loads(f.read_text(encoding='utf8'))['features'])
 if len(fs)!=6999 or len({x['properties']['GEOID'] for x in fs})!=6999:raise RuntimeError('Census pages incomplete or duplicated')
 complete=[]
 for f in fs:
  g=f['geometry'];polys=[g['coordinates']] if g['type']=='Polygon' else g['coordinates'];pts=[p for poly in polys for ring in poly for p in ring]
  xs=[p[0] for p in pts];ys=[p[1] for p in pts]
  if min(xs)>=ROI[0]-1e-10 and max(xs)<=ROI[2]+1e-10 and min(ys)>=ROI[1]-1e-10 and max(ys)<=ROI[3]+1e-10:complete.append(f)
 return fs,complete
def main():
 z=np.load(PROJ/'results'/'dwr'/'grid.npz');x=z['x_edges'];y=z['y_edges'];pop=z['population'];area=z['area_m2'];support=z['support'];nr,nc=pop.shape;dx=float(np.median(np.diff(x)));dy=float(np.median(np.diff(y)))
 all_blocks,blocks=load_blocks();nb=len(blocks);dwr=json.loads((PROJ/'external'/'dwr'/'roi_features.json').read_text(encoding='utf8'))['features'];dwr_shapes=[]
 for f in dwr:
  rings=f['geometry']['rings']
  if len(rings)!=1:raise RuntimeError('unexpected DWR polygon ring')
  dwr_shapes.append(({'type':'Polygon','coordinates':[rings[0]]},1))
 rows=[]
 for factor in [5,10]:
  shp=(nr*factor,nc*factor);tr=Affine(dx/factor,0,x[0],0,-dy/factor,y[-1]);bshapes=[(f['geometry'],i+1) for i,f in enumerate(blocks)]
  labels=rasterize(bshapes,out_shape=shp,transform=tr,fill=0,dtype='int32',all_touched=False)
  obs=rasterize(dwr_shapes,out_shape=shp,transform=tr,fill=0,dtype='uint8',all_touched=False).astype(bool)
  lon=x[0]+(np.arange(nc*factor)+.5)*dx/factor;lat=y[-1]-(np.arange(nr*factor)+.5)*dy/factor
  inroi=((lon>=ROI[0])&(lon<=ROI[2]))[None,:]&((lat>=ROI[1])&(lat<=ROI[3]))[:,None]
  weights=R**2*np.deg2rad(dx/factor)*np.deg2rad(dy/factor)*np.cos(np.deg2rad(lat))[:,None]*inroi
  ba=np.bincount(labels.ravel(),weights=weights.ravel(),minlength=nb+1)[1:];oa=np.bincount(labels.ravel(),weights=(weights*obs).ravel(),minlength=nb+1)[1:]
  cov=np.divide(oa,ba,out=np.zeros(nb),where=ba>0);cp=np.array([f['properties']['P0010001'] or 0 for f in blocks],dtype='float64');aland=np.array([f['properties']['ALAND'] or 0 for f in blocks],dtype='float64');awater=np.array([f['properties']['AWATER'] or 0 for f in blocks],dtype='float64');tot=cp.sum()
  cshare=float(np.dot(cp,cov)/tot);medshare=float(cp[cov>=.5].sum()/tot);nineshare=float(cp[cov>=.9].sum()/tot)
  lab4=np.flipud(labels).reshape(nr,factor,nc,factor);obs4=np.flipud(obs).reshape(nr,factor,nc,factor);w4=np.flipud(weights).reshape(nr,factor,nc,factor)
  den=np.sum(w4,axis=(1,3));bw=np.sum(w4*(lab4>0),axis=(1,3));bwo=np.sum(w4*((lab4>0)&obs4),axis=(1,3));fr=np.divide(bw,den,out=np.zeros_like(bw),where=den>0);fro=np.divide(bwo,den,out=np.zeros_like(bwo),where=den>0)
  gp=float(np.sum(pop*fr));gpo=float(np.sum(pop*fro));gs=gpo/gp if gp else float('nan')
  rows.append(dict(supersampling=factor,all_Census_blocks_intersecting_ROI=len(all_blocks),complete_blocks=nb,complete_block_Census_population=int(tot),complete_block_land_area_km2=aland.sum()/1e6,rasterized_complete_block_area_km2=ba.sum()/1e6,median_rasterized_to_Census_area_ratio=float(np.median(np.divide(ba,aland+awater,out=np.zeros(nb),where=(aland+awater)>0))),blocks_with_rasterized_area_pct=100*np.mean(ba>0),DWR_observation_polygons=len(dwr),DWR_area_in_complete_blocks_km2=oa.sum()/1e6,Census_population_weighted_uniform_support_share_pct=100*cshare,Census_population_in_blocks_with_median_observed_area_ge_50pct_pct=100*medshare,Census_population_in_blocks_with_observed_area_ge_90pct_pct=100*nineshare,native_GHSL_population_in_complete_blocks=gp,native_GHSL_population_weighted_DWR_support_share_pct=100*gs,difference_GHSL_vs_Census_uniform_support_share_pp=100*(gs-cshare),median_block_DWR_area_coverage_pct=100*np.median(cov),blocks_with_any_DWR_area=int(np.sum(cov>0)),blocks_with_ge50pct_DWR_area=int(np.sum(cov>=.5)),blocks_with_ge90pct_DWR_area=int(np.sum(cov>=.9))))
  if factor==10:
   pd.DataFrame(dict(GEOID=[f['properties']['GEOID'] for f in blocks],census_2020_P0010001=cp,ALAND_m2=aland,AWATER_m2=awater,rasterized_area_m2=ba,DWR_observed_area_m2=oa,DWR_observed_area_share=cov,Census_uniform_supported_population=cp*cov,median_observable_block=(cov>=.5),Census_population_in_median_observable_blocks=cp*(cov>=.5))).to_csv(OUT/'fresno_block_support_2020.csv',index=False)
  print('done',factor,flush=True)
 pd.DataFrame(rows).to_csv(OUT/'fresno_block_support_summary.csv',index=False)
 hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(BLOCK_DIR.glob('*.geojson'))}
 meta=dict(source_service='Esri ArcGIS Living Atlas USA_Census_2020_DHC_Blocks, Block layer; Census 2020 TIGER geography; population P0010001',query='GEOID prefix 06019, exact project ROI, WGS84 intersect; 6,999 unique block GEOIDs; complete block subset 6,788',privacy='Census DHC counts use Census Bureau differential privacy per service metadata; GHSL may share census-derived inputs.',epochs='Census/GHSL 2020; DWR published observation polygons 2015-2026; this assesses support only, no velocity field.',geometry_method='Census block and DWR rectangle polygons rasterized at 5x and 10x the 3-arcsec GHSL grid; spherical row-area weights; complete blocks retained to avoid including population outside ROI.',estimands='Census population weighted mean DWR observed-area fraction under uniform within-block allocation; population in blocks with >=50% or >=90% observed area; native GHSL population weighted support on same complete-block domain.',limitations=['Not an Ohenhen-style median velocity exposure classification: DWR polygons contain no per-cell LOS velocities.','Subpixel rasterization approximates polygon intersections; factor sensitivity is reported.','Census counts are privacy-protected; population-source independence is unverified.','Complete-block restriction excludes boundary-crossing Census blocks.'],ROI=list(ROI),GHSL_grid_sha256=hashlib.sha256((PROJ/'results'/'dwr'/'grid.npz').read_bytes()).hexdigest(),Census_response_sha256=hashes,DWR_feature_count=len(dwr))
 (OUT/'fresno_block_support_method.json').write_text(json.dumps(meta,indent=2),encoding='utf8');print(pd.DataFrame(rows).to_string(index=False),flush=True)
if __name__=='__main__':main()
