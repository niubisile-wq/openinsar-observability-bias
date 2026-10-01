import os
from pathlib import Path
import json, csv
import numpy as np
import rasterio
from rasterio.features import rasterize
from dwr_units import feet_to_feet_per_year,duration_years,WINDOWS,DAYS_PER_YEAR

BASE = Path(os.environ.get('INSAR_FOLLOWUP_ROOT', Path(__file__).resolve().parent))
PROJ = Path(os.environ['INSAR_STRENGTHENING_ROOT'])
BLOCK_DIR = BASE / 'census_blocks' / 'fresno_roi'
OUT = BASE / 'census_block_velocity_analysis'
OUT.mkdir(exist_ok=True)
ROI = (-119.9, 36.6, -119.7, 36.85)
PERIODS = {
    '2015-10 to 2016-10 annual': ('rate', BASE / 'dwr_velocity' / 'rate_2015_2016.tif'),
    '2020-10 to 2021-10 annual': ('rate', BASE / 'dwr_velocity' / 'rate_2020_2021.tif'),
    '2023-10 to 2024-10 annual': ('rate', BASE / 'dwr_velocity' / 'rate_2023_2024.tif'),
    '2025-04 to 2026-04 annual': ('rate', BASE / 'dwr_velocity' / 'rate_2025_2026.tif'),
    '2015-10 to 2021-10 temporal median of annual rates': ('rate', BASE / 'dwr_velocity' / 'rate_median_2015_2021.tif'),
    '2015-06 to 2021-01 cumulative average rate': ('cumulative', BASE / 'dwr_velocity' / 'total_since_20150613_20210101.tif'),
}
YEARS_2015_2021 = duration_years('total_since_20150613_20210101.tif')
THRESHOLDS = [-10.0, -5.0, -3.0, 0.0]  # mm yr-1; classes follow Ohenhen et al. 2025

def read_blocks():
    feats = []
    for p in sorted(BLOCK_DIR.glob('*.geojson')):
        feats.extend(json.loads(p.read_text(encoding='utf-8'))['features'])
    if len(feats) != 6999 or len({f['properties']['GEOID'] for f in feats}) != 6999:
        raise RuntimeError('Census GeoJSON pages are incomplete or duplicated')
    keep = []
    for f in feats:
        coords = f['geometry']['coordinates']
        polys = [coords] if f['geometry']['type'] == 'Polygon' else coords
        pts = [pt for poly in polys for ring in poly for pt in ring]
        if (min(p[0] for p in pts) >= ROI[0]-1e-10 and max(p[0] for p in pts) <= ROI[2]+1e-10
            and min(p[1] for p in pts) >= ROI[1]-1e-10 and max(p[1] for p in pts) <= ROI[3]+1e-10):
            keep.append(f)
    if len(keep) != 6788:
        raise RuntimeError(f'Expected 6,788 complete Census blocks, found {len(keep)}')
    return keep

def rings_geojson(f):
    rings = f['geometry']['rings']
    return {'type': 'Polygon', 'coordinates': rings}

def classify(v):
    if v >= 0: return '>=0'
    if v >= -3: return '[-3,0)'
    if v >= -5: return '[-5,-3)'
    if v >= -10: return '[-10,-5)'
    return '<-10'

def main():
    blocks = read_blocks()
    dwr = json.loads((PROJ/'external'/'dwr'/'roi_features.json').read_text(encoding='utf-8'))['features']
    with rasterio.open(next(iter(PERIODS.values()))[1]) as ref:
        profile = ref.profile.copy(); transform = ref.transform; shape = (ref.height, ref.width)
        if (ref.width, ref.height) != (177, 279): raise RuntimeError('Unexpected DWR raster dimensions')
        if abs(transform.a-0.001133584162363299) > 1e-10: raise RuntimeError('Raster grid not native-aligned')
    block_shapes = []
    block_props = []
    for i, f in enumerate(blocks, 1):
        block_shapes.append((f['geometry'], i))
        block_props.append((f['properties']['GEOID'], int(f['properties']['P0010001']), int(f['properties']['ALAND']), int(f['properties']['AWATER'])))
    block_id = rasterize(block_shapes, out_shape=shape, transform=transform, fill=0, dtype='int32', all_touched=False)
    support_shapes = [(rings_geojson(f), 1) for f in dwr]
    support = rasterize(support_shapes, out_shape=shape, transform=transform, fill=0, dtype='uint8', all_touched=False).astype(bool)
    if int(support.sum()) != len(dwr):
        raise RuntimeError(f'DWR cell/pixel alignment check failed: {support.sum()} cells for {len(dwr)} polygons')
    # Also preserve the points as a CSV so the resulting block statistics can be independently audited.
    summaries=[]; block_rows=[]; pixel_rows=[]
    geoid = np.array([p[0] for p in block_props]); pop = np.array([p[1] for p in block_props], dtype=np.int64)
    for label, (kind, path) in PERIODS.items():
        with rasterio.open(path) as ds:
            arr = ds.read(1).astype('float64')
            nodata = ds.nodata if ds.nodata is not None else -9999
            valid_arr = np.isfinite(arr) & (arr != nodata)
            rate_ft_yr = feet_to_feet_per_year(arr,path.name)
            valid = support & valid_arr & (block_id > 0)
            r,c = np.where(valid); bidx = block_id[r,c]-1
            mm = rate_ft_yr[r,c] * 304.8  # DWR raw values are feet; converted to mm/year.
            vals = [[] for _ in blocks]
            for bi,v in zip(bidx,mm): vals[int(bi)].append(float(v))
            med = np.array([np.median(v) if v else np.nan for v in vals])
            npx = np.array([len(v) for v in vals])
            cls = np.array([classify(v) if np.isfinite(v) else 'no valid pixels' for v in med])
            included = npx > 0
            total_pop = int(pop.sum()); covered_pop = int(pop[included].sum())
            counts={k:int(pop[included & (cls==k)].sum()) for k in ['>=0','[-3,0)','[-5,-3)','[-10,-5)','<-10']}
            pixel_class=np.where(mm>=0,0,np.where(mm>=-3,1,np.where(mm>=-5,2,np.where(mm>=-10,3,4))))
            # Assign each observed pixel's share of a block population uniformly within that block.
            pixel_pop=np.divide(pop[bidx],npx[bidx],where=npx[bidx]>0)
            pixel_mass={k:float(pixel_pop[pixel_class==j].sum()) for j,k in enumerate(['>=0','[-3,0)','[-5,-3)','[-10,-5)','<-10'])}
            summaries.append(dict(period=label,source_pixels_on_DWR_footprint=int(valid.sum()),blocks_with_valid_pixels=int(included.sum()),complete_block_population=int(total_pop),population_in_blocks_with_valid_pixels=covered_pop,share_of_complete_Census_population_in_valid_blocks_pct=covered_pop/total_pop*100,median_pixel_rate_mm_per_year=float(np.median(mm)),population_assigned_by_block_median_ge_0=counts['>=0'],population_assigned_by_block_median_minus3_to_0=counts['[-3,0)'],population_assigned_by_block_median_minus5_to_minus3=counts['[-5,-3)'],population_assigned_by_block_median_minus10_to_minus5=counts['[-10,-5)'],population_assigned_by_block_median_lt_minus10=counts['<-10'],uniform_within_block_population_ge_0=pixel_mass['>=0'],uniform_within_block_population_minus3_to_0=pixel_mass['[-3,0)'],uniform_within_block_population_minus5_to_minus3=pixel_mass['[-5,-3)'],uniform_within_block_population_minus10_to_minus5=pixel_mass['[-10,-5)'],uniform_within_block_population_lt_minus10=pixel_mass['<-10']))
            for i in range(len(blocks)):
                block_rows.append(dict(period=label,GEOID=geoid[i],census_population_2020=pop[i],valid_DWR_rate_pixels=npx[i],median_vertical_rate_mm_yr=med[i],rate_class=cls[i]))
            for rr,cc,v,bi in zip(r,c,mm,bidx): pixel_rows.append(dict(period=label,row=int(rr),col=int(cc),GEOID=geoid[bi],vertical_rate_mm_yr=float(v)))
    for name, rows in [('fresno_census_block_velocity_summary.csv',summaries),('fresno_census_block_velocity_by_block.csv',block_rows),('fresno_dwr_rate_pixels_by_block.csv',pixel_rows)]:
        with (OUT/name).open('w',newline='',encoding='utf-8-sig') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    metadata={'DWR_service':'DWR TRE ALTAMIRA annual rate and cumulative vertical displacement ImageServer, raw GeoTIFF export','raster_grid':'native grid alignment: 177 x 279, EPSG:4326, 0.0011335841623633 degree x 0.0008998200359928 degree; nearest-neighbor export snapped to service pixel boundaries with adjustAspectRatio=false','DWR_support_geometry':'38,738 DWR 100 m observation-cell polygons from project external/dwr/roi_features.json, rasterized with pixel-centre inclusion; exact validation gives 38,738 polygons to 38,738 raster cells','population':'2020 Census DHC P1 total (P0010001), Census block GEOID prefix 06019; 6,788 complete blocks','estimand':'median DWR vertical displacement rate across valid published DWR observation-cell centers assigned to each Census block; whole block population is assigned to this median-rate class','classification_mm_yr':{'VLM >= 0':'not subsiding','-3 <= VLM < 0':'0 to -3','-5 <= VLM < -3':'-3 to -5','-10 <= VLM < -5':'-5 to -10','VLM < -10':'below -10'},'threshold_source':'Ohenhen et al., Nature Cities 2, 543-554 (2025), https://doi.org/10.1038/s44284-025-00240-y','limitations':['DWR rasters are interpolated products from DWR/TRE ALTAMIRA Sentinel-1 observations, not the manuscript original LOS field; this is a transfer benchmark, not independent InSAR validation.','DWR values are vertical rates (feet/year converted to mm/year), not LOS rates; they share provider lineage with the DWR observation footprint.','Raster values are sampled only at published DWR 100m observation-cell footprints, then median-aggregated within Census blocks.','Observation-cell locations and annual-rate rasters have slightly different product vintages; this is documented as a cross-product sensitivity, not a pixelwise accuracy assessment.','The 2015-2021 temporal median summarizes six overlapping annual-rate rasters; the cumulative end-point rate is a different temporal estimand from a fitted linear trend.','Census DHC population is differentially private and population source independence is unverified.','A complete Census block with no DWR observation cells is excluded from rate-class population assignment and reported as uncovered.'],'periods':list(PERIODS.keys()),'raster_pixel_count':int(shape[0]*shape[1]),'DWR_cell_feature_count':len(dwr),'DWR_rasterized_support_cell_count':int(support.sum()),'complete_block_count':len(blocks),'complete_block_population':int(pop.sum()),'cumulative_rate_duration_years':YEARS_2015_2021}
    metadata['unit_convention']='Raw feet for each documented interval, divided by days/365.2425 before classification; the six-window composite is already annualized.'
    metadata['input_windows']=WINDOWS
    metadata['days_per_year']=DAYS_PER_YEAR
    metadata['limitations'][1]='DWR interval displacements are annualized vertical quantities, not LOS; provider lineage is shared with the footprint.'
    metadata['limitations'][4]='The temporal median uses six consecutive nonoverlapping October-to-October windows; the cumulative endpoint rate is a different estimand.'
    (OUT/'method.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
    print(json.dumps(summaries,indent=2))

if __name__=='__main__': main()
