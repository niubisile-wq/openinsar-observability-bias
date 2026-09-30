import os
"""Independent exact spherical-area overlap audit for 1-km UTM grids.

Each source 0.001-degree cell is represented in lon-radians/sin(latitude), an
equal-area spherical coordinate system. Densified inverse-UTM cell edges are
intersected directly with those source cells, independent of raster reprojection.
"""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from shapely.geometry import box,Polygon
from shapely import STRtree
from pyproj import Transformer
from affine import Affine
import rasterio
from rasterio.warp import Resampling,reproject

HERE = Path(os.environ.get('INSAR_PAIRED_ROOT', Path(__file__).resolve().parent))
STUDY=Path(os.environ['INSAR_STRENGTHENING_ROOT'])
SRC=STUDY/'results'/'timeseries'/'matched_quality.npz'; GEOM=STUDY/'results'/'timeseries'/'geometry.json'
PILOT=HERE/'paired_support_pilot'; OUT=HERE/'direct_geometry_audit'
CRS='EPSG:32647'; R=6371008.8

def main():
 OUT.mkdir(exist_ok=True); geom=json.loads(GEOM.read_text(encoding='utf-8')); tr=Affine(*geom['transform'])
 with np.load(SRC) as z:
  xe=z['x_edges']; ye=z['y_edges']; p=z['population'].astype(float); a=z['area_m2'].astype(float)
  fs={'matched437_mean_support':z['matched_support'].astype(float),'final_quality_mask':z['full_valid'].astype(float)}
 # Source cells in equal-area coordinates, north-to-south ordering.
 lon=np.deg2rad(xe); lat=np.sin(np.deg2rad(ye)); polys=[]; index=[]
 for r in range(p.shape[0]):
  j=p.shape[0]-1-r
  for c in range(p.shape[1]):
   polys.append(box(lon[c],lat[j],lon[c+1],lat[j+1])); index.append((r,c))
 tree=STRtree(polys); index=np.asarray(index); togeo=Transformer.from_crs(CRS,'EPSG:4326',always_xy=True)
 # Make four grid origins, match the 1-km primary scale.
 x0,y0,x1,y1=rasterio.warp.transform_bounds('EPSG:4326',CRS,xe[0],ye[0],xe[-1],ye[-1],densify_pts=100)
 rows=[]; accum={n:dict(P=0.,A=0.,AP=0.,PP=0.) for n in fs}
 for sx,sy in [(0.,0.),(500.,0.),(0.,500.),(500.,500.)]:
  xs=np.arange(np.floor((x0-sx)/1000)*1000+sx,np.ceil((x1-sx)/1000)*1000+sx+1,1000)
  ys=np.arange(np.floor((y0-sy)/1000)*1000+sy,np.ceil((y1-sy)/1000)*1000+sy+1,1000)
  for yy in ys[:-1]:
   for xx in xs[:-1]:
    # Densify all four projected straight edges before inverse projection.
    n=24; e=[]
    for t in np.linspace(0,1,n,endpoint=False):e.append((xx+1000*t,yy))
    for t in np.linspace(0,1,n,endpoint=False):e.append((xx+1000,yy+1000*t))
    for t in np.linspace(0,1,n,endpoint=False):e.append((xx+1000*(1-t),yy+1000))
    for t in np.linspace(0,1,n,endpoint=False):e.append((xx,yy+1000*(1-t)))
    ex,ey=togeo.transform(*np.array(e).T); cell=Polygon(np.column_stack([np.deg2rad(ex),np.sin(np.deg2rad(ey))]))
    if not cell.is_valid: cell=cell.buffer(0)
    ids=tree.query(cell,predicate='intersects'); P=A=0.; numer={n:0. for n in fs}; anumer={n:0. for n in fs}
    for q in ids:
     r,c=index[q]; ov=cell.intersection(polys[q]).area
     if ov<=0:continue
     fraction=ov/polys[q].area; P+=p[r,c]*fraction; A+=a[r,c]*fraction
     for name,f in fs.items():
      numer[name]+=p[r,c]*fraction*f[r,c];anumer[name]+=a[r,c]*fraction*f[r,c]
    if P<=0:continue
    for name in fs:
     rows.append(dict(origin_x_m=sx,origin_y_m=sy,cell_x=xx,cell_y=yy,population=P,area_m2=A,
       support=name,population_supported=numer[name],area_supported=anumer[name]))
  print('origin complete',sx,sy,flush=True)
 # Summarize direct integration, then compare with the independent rasterio 25m
 # aggregate source totals in the pilot ledger. The direct calculation is also
 # the geometric reference for uniform-area allocation effects.
 out=[]
 for name,f in fs.items():
  native=float((p*(1-f)).sum()); total=float(p.sum())
  for (sx,sy),d in pd.DataFrame([r for r in rows if r['support']==name]).groupby(['origin_x_m','origin_y_m']):
   valid_area=d.area_m2>0; pop=float(d.population.sum()); supp=float(d.population_supported.sum())
   A=float(d.area_m2.sum()); As=float(d.area_supported.sum()); uniform=float((d.population*(1-np.divide(d.area_supported,d.area_m2,out=np.zeros(len(d)),where=d.area_m2>0))).sum())
   out.append(dict(support=name,origin_x_m=sx,origin_y_m=sy,source_population=total,geometric_population=pop,
     population_mass_drift=pop-total,source_native_deficit=native,geometric_supported_population=supp,
     geometric_native_deficit=pop-supp,deficit_pct_total=100*(pop-supp)/total,
     uniform_area_deficit=uniform,uniform_difference_pp_total=100*(uniform-(pop-supp))/total,
     geometric_area_m2=A,source_area_m2=float(a.sum()),area_mass_drift=A-float(a.sum()),
     area_support=float(As/A)))
 pd.DataFrame(rows).to_csv(OUT/'cell_intersections.csv',index=False)
 df=pd.DataFrame(out);df.to_csv(OUT/'summary.csv',index=False)
 meta=dict(source=str(SRC),source_sha256=hashlib.sha256(SRC.read_bytes()).hexdigest(),geometry='Exact Shapely intersections in spherical equal-area lon-radians/sine-latitude coordinates; 24-segment densification per inverse-projected UTM cell edge; source-cell uniform density/support assumption',crs=CRS,cell_size_m=1000,origins=[[0,0],[500,0],[0,500],[500,500]],population_total=float(p.sum()),area_total_m2=float(a.sum()),cell_count=len(rows),interpretation='Independent source-to-coarse-grid geometry cross-check; uniformity within each 0.001-degree input cell is assumed. Not independent household or motion truth.')
 (OUT/'method.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False),encoding='utf-8')
 print(df.to_string(index=False))
if __name__=='__main__':main()
