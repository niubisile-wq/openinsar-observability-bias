"""Exact polygon intersections for conditional DWR masking on official Census blocks."""
from pathlib import Path
import os
import sys
import json
import hashlib
import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.warp import reproject,Resampling
from affine import Affine
import shapely
from shapely.strtree import STRtree
from pyproj import Transformer

ROOT=Path(__file__).resolve().parent.parent
REV=ROOT.parent/'github_original_source/revision_2026'
STUDY=ROOT.parent/'01_当前返修工作区/InSAR_revision_20260926/strengthening'
OUT=ROOT/'inputs/fresno_real_blocks'
os.environ['INSAR_STRENGTHENING_ROOT']=str(STUDY)
os.environ['INSAR_FOLLOWUP_ROOT']=str(ROOT/'intermediates')
sys.path.insert(0,str(REV/'followup'))
from fresno_capop_support_sensitivity import census_counts,complete_blocks

def main():
    blocks,records=complete_blocks(census_counts());b=np.asarray(blocks,dtype=object)
    p=np.array([v for _,v in records],dtype=float)
    with rasterio.open(OUT/'DWR_20201001_20211001.tif') as ds:
        field=ds.read(1);tr=ds.transform;h,w=field.shape
        valid=np.isfinite(field)&(field!=ds.nodata)
    duration=json.loads((OUT/'DWR_rate_provenance.json').read_text(encoding='utf-8'))['duration_years']
    rates=field.astype(float)*304.8/duration
    footprints=json.loads((STUDY/'external/dwr/roi_features.json').read_text(encoding='utf-8'))['features']
    footprint=rasterize([(dict(type='Polygon',coordinates=f['geometry']['rings']),1) for f in footprints],
                        out_shape=(h,w),transform=tr,fill=0,dtype='uint8',all_touched=False).astype(bool)
    assert int(footprint.sum())==38738
    valid&=footprint
    src=STUDY/'results/dwr/grid.npz'
    with np.load(src) as z: pop=z['population'];xe=z['x_edges'];ye=z['y_edges']
    src_tr=Affine(xe[1]-xe[0],0,xe[0],0,-(ye[1]-ye[0]),ye[-1])
    ghs=np.zeros((h,w),dtype=float)
    reproject(pop[::-1].copy(),ghs,src_transform=src_tr,src_crs='EPSG:4326',
              dst_transform=tr,dst_crs='EPSG:4326',resampling=Resampling.sum,src_nodata=None,dst_nodata=0)
    massdiff=abs(float(ghs.sum())-float(pop.sum()))/float(pop.sum())
    if massdiff>.005: raise RuntimeError('GHSL reprojection mass check failed: '+str(massdiff))
    rows,cols=np.indices((h,w));lon=tr.c+(cols.ravel()+.5)*tr.a;lat=tr.f+(rows.ravel()+.5)*tr.e
    cells=shapely.box(lon-tr.a/2,lat+tr.e/2,lon+tr.a/2,lat-tr.e/2)
    projection=Transformer.from_crs('EPSG:4326','EPSG:3310',always_xy=True)
    cells=shapely.transform(cells,lambda xy:np.column_stack(projection.transform(xy[:,0],xy[:,1])))
    points=shapely.points(*projection.transform(lon,lat))
    tree=STRtree(b)
    ci,bi=tree.query(cells,predicate='intersects')
    area=shapely.area(shapely.intersection(cells[ci],b[bi]))
    keep=area>1e-6;ci=ci[keep];bi=bi[keep];area=area[keep]
    cell_area=shapely.area(cells)
    block_area=np.bincount(bi,weights=area,minlength=len(b))
    block_geom_area=shapely.area(b)
    relative=np.abs(block_area-block_geom_area)/block_geom_area
    assert float(relative.max())<1e-5, relative.max()
    graw=ghs.ravel()[ci]*area/cell_area[ci]
    sums=np.bincount(bi,weights=graw,minlength=len(b))
    eligible=(sums>0)&(p>0)
    allocation=np.zeros(len(ci));ok=eligible[bi]
    allocation[ok]=p[bi[ok]]*graw[ok]/sums[bi[ok]]
    total=float(allocation.sum());assert abs(total-p[eligible].sum())<1e-6
    available=ok&valid.ravel()[ci]
    cindex,bindex=tree.query(points,predicate='within')
    assert len(np.unique(cindex))==len(cindex)
    cp=eligible[bindex]&valid.ravel()[cindex]
    cindex=cindex[cp];bindex=bindex[cp]
    np.savez_compressed(OUT/'real_block_masking_inputs.npz',
                        pixel_shape=np.array([h,w]),rates_mm_per_year=rates,
                        original_valid=valid,intersection_pixel=ci[available],intersection_block=bi[available],
                        intersection_area_m2=area[available],allocation_weight=allocation[available],
                        center_pixel=cindex,center_block=bindex,
                        block_census_population=p,block_eligible=eligible,
                        block_geoid=np.array([v[0] for v in records]),
                        pixel_x_m=shapely.get_x(points).reshape(h,w),pixel_y_m=shapely.get_y(points).reshape(h,w))
    metadata={'Census_complete_blocks':len(records),'Census_complete_population':float(p.sum()),
              'Census_population_with_positive_GHSL_mass':float(p[eligible].sum()),
              'population_allocation_in_known_rate_footprint':float(allocation[available].sum()),
              'population_allocation_outside_known_rate_footprint':float(total-allocation[available].sum()),
              'known_DWR_rate_pixels':int(valid.sum()),'exact_intersections_in_known_rate_footprint':int(available.sum()),
              'GHSL_reprojection_relative_mass_error':massdiff,'max_polygon_relative_area_closure_error':float(relative.max()),
              'source_sha256':{str(path.relative_to(ROOT.parent)).replace('\\','/'):hashlib.sha256(path.read_bytes()).hexdigest() for path in
                               [src,OUT/'DWR_20201001_20211001.tif',STUDY/'external/round3/ca2020.pl.zip',STUDY/'external/round3/census_blocks/tl_2020_06019_tabblock20.shp']},
              'estimand':'GHSL-normalized Census allocations restricted to known DWR-rate cell intersections; original unobserved locations excluded from truth.',
              'baseline_rule':'Normalize GHSL over complete blocks before restricting to known-rate footprint; do not force all Census population into observed cells.',
              'median_rule':'Unweighted published-rate pixel-center median within each irregular block; absence of a center remains undefined.',
              'limits':['DWR interpolated vertical products and footprint have shared lineage and different vintages; this is a conditional method transfer, not independent physical validation.',
                        'GHSL remains modeled population, not household ground truth.',
                        'Rates are constant inside their published cells; GHSL is resampled conservatively onto that grid.',
                        'Original missing rates are never assigned known reference classes.']}
    (OUT/'real_block_geometry_protocol.json').write_bytes((json.dumps(metadata,indent=2)+'\n').encode())
    print(json.dumps({k:v for k,v in metadata.items() if k not in ['source_sha256','limits']}),flush=True)

if __name__=='__main__': main()
