import os
"""Compare 1-km population grids with 2020 Census block totals in Fresno.

This is a coarse-unit benchmark. Population counts from raster cells are
allocated to block polygons by intersection area, assuming uniform density
within each raster cell. It is not independent household-level truth.
"""
from pathlib import Path
import glob,json,hashlib,zipfile
import numpy as np,pandas as pd,rasterio
from rasterio.io import MemoryFile
from rasterio.windows import from_bounds
from pyproj import Transformer
from shapely.geometry import shape,box,Polygon
from shapely.ops import transform as shp_transform
from shapely.strtree import STRtree

ROOT=Path(os.environ.get('INSAR_FOLLOWUP_ROOT', Path(__file__).resolve().parent))
BLOCK_DIR=ROOT/'census_blocks'/'fresno_roi'
INPUT=ROOT/'inputs'/'fresno_population'
OUT=ROOT/'fresno_population_benchmark'
OUT.mkdir(exist_ok=True)
to_54009=Transformer.from_crs('EPSG:4326','ESRI:54009',always_xy=True)

def load_blocks():
    geoms=[]; rows=[]
    roi=(-119.9,36.6,-119.7,36.85)
    for f in sorted(BLOCK_DIR.glob('blocks_*.geojson')):
        doc=json.loads(f.read_text(encoding='utf-8'))
        for feat in doc['features']:
            pr=feat['properties']; g=shape(feat['geometry'])
            # Match the published Fresno complete-block domain exactly.
            x0,y0,x1,y1=g.bounds
            if x0<roi[0]-1e-10 or y0<roi[1]-1e-10 or x1>roi[2]+1e-10 or y1>roi[3]+1e-10:
                continue
            if g.is_empty or not g.is_valid: g=g.buffer(0)
            gp=shp_transform(to_54009.transform,g)
            geoms.append(gp)
            rows.append(dict(GEOID=str(pr['GEOID']),Census_population=int(pr['P0010001'] or 0),
                             land_area_m2=float(pr['ALAND'] or 0),water_area_m2=float(pr['AWATER'] or 0)))
    df=pd.DataFrame(rows)
    if df.GEOID.duplicated().any(): raise ValueError('Duplicate Census block GEOID')
    return df,geoms

def overlay_grid(values,transform,crs,geoms,tree,df,kind,nodata):
    """Area-weight each population cell into intersecting Census blocks."""
    out=np.zeros(len(df),dtype='float64')
    h,w=values.shape
    if str(crs)=='ESRI:54009':
        pixw=float(transform.a); pixh=abs(float(transform.e))
        if not np.isclose(pixw,pixh): raise ValueError('Expected square GHSL pixels')
        for r in range(h):
            row=values[r]
            for c in np.flatnonzero(np.isfinite(row)&(row!=nodata)&(row>0)):
                x0=transform.c+c*transform.a; y1=transform.f+r*transform.e
                cell=box(x0,y1-pixh,x0+pixw,y1); cellarea=cell.area
                ix=tree.query(cell,predicate='intersects')
                if len(ix)==0: continue
                shares=[]
                for i in ix:
                    frac=cell.intersection(geoms[int(i)]).area/cellarea
                    if frac>1e-12: out[int(i)]+=row[c]*frac; shares.append(frac)
                if sum(shares)>1.00001: raise ValueError(('Overlapping Census geometry',kind,r,int(c),sum(shares)))
    elif str(crs)=='EPSG:4326':
        for r in range(h):
            row=values[r]
            for c in np.flatnonzero(np.isfinite(row)&(row!=nodata)&(row>0)):
                x0=transform.c+c*transform.a; x1=x0+transform.a
                y1=transform.f+r*transform.e; y0=y1+transform.e
                # Densify the geographic cell perimeter before projection; straight
                # edges in longitude/latitude are curved in the equal-area CRS.
                xs=np.linspace(min(x0,x1),max(x0,x1),9)
                ys=np.linspace(min(y0,y1),max(y0,y1),9)
                ring=[(x,min(y0,y1)) for x in xs]+[(max(x0,x1),y) for y in ys[1:]]+[(x,max(y0,y1)) for x in xs[-2::-1]]+[(min(x0,x1),y) for y in ys[-2:0:-1]]
                cell_geo=Polygon(ring)
                cell=shp_transform(to_54009.transform,cell_geo); cellarea=cell.area
                ix=tree.query(cell,predicate='intersects')
                if len(ix)==0: continue
                shares=[]
                for i in ix:
                    frac=cell.intersection(geoms[int(i)]).area/cellarea
                    if frac>1e-12: out[int(i)]+=row[c]*frac; shares.append(frac)
                if sum(shares)>1.00001: raise ValueError(('Overlapping Census geometry',kind,r,int(c),sum(shares)))
    else: raise ValueError(f'Unexpected CRS {crs}')
    return out

def metrics(obs,pred):
    err=pred-obs; total=float(obs.sum())
    corr=float(np.corrcoef(obs,pred)[0,1]) if np.std(obs)>0 and np.std(pred)>0 else float('nan')
    rank=float(pd.Series(obs).corr(pd.Series(pred),method='spearman'))
    return dict(predicted_total=float(pred.sum()),total_bias=float(pred.sum()-total),total_bias_pct=100*float(pred.sum()-total)/total,
        MAE_per_block=float(np.mean(np.abs(err))),RMSE_per_block=float(np.sqrt(np.mean(err**2))),
        WAPE_pct=100*float(np.abs(err).sum())/total,mean_signed_bias_per_block=float(err.mean()),
        Pearson_r=corr,Spearman_r=rank,blocks_with_prediction_zero=int(np.sum((pred==0)&(obs>0))))

def main():
    df,geoms=load_blocks()
    if len(df)!=6788 or int(df.Census_population.sum())!=523953:
        raise ValueError(f'Unexpected complete-block domain: {len(df)} blocks, {df.Census_population.sum()} Census residents')
    obs=df.Census_population.to_numpy(float)
    total=float(obs.sum()); print('blocks',len(df),'Census total',total,'land km2',df.land_area_m2.sum()/1e6)
    tree=STRtree(geoms)
    # GHSL: one 1-km Mollweide tile containing the complete study footprint.
    z=zipfile.ZipFile(INPUT/'ghsl_2020_1km_fresno_tile.zip')
    tif=next(x for x in z.namelist() if x.lower().endswith('.tif'))
    mem=MemoryFile(z.read(tif)); src=mem.open(); ghsrc=src
    xmin,ymin,xmax,ymax=zip(*(g.bounds for g in geoms))
    b=(min(xmin),min(ymin),max(xmax),max(ymax))
    win=from_bounds(*b,transform=src.transform).round_offsets().round_lengths()
    gh=src.read(1,window=win); ght=src.window_transform(win)
    if src.nodata is not None: gh=np.where(gh==src.nodata,0,gh)
    gh_counts=overlay_grid(gh,ght,src.crs,geoms,tree,df,'GHSL_2020_1km',src.nodata)
    src.close();mem.close()
    # WorldPop: read only the small Fresno window from its public 1-km US GeoTIFF.
    wp_path=INPUT/'worldpop_2020_1km_us_unadj.tif'
    with rasterio.open(wp_path) as src:
        # Inverse-transform study bounds from the equal-area block envelope.
        from_54009=Transformer.from_crs('ESRI:54009','EPSG:4326',always_xy=True)
        wgs=from_54009.transform_bounds(*b,densify_pts=64)
        win=from_bounds(*wgs,transform=src.transform).round_offsets().round_lengths()
        win=win.intersection(rasterio.windows.Window(0,0,src.width,src.height))
        wp=src.read(1,window=win); wpt=src.window_transform(win)
        wp_counts=overlay_grid(wp,wpt,src.crs,geoms,tree,df,'WorldPop_2020_1km_UNadj',src.nodata)
    df['GHSL_E2020_1km_block_count']=gh_counts
    df['WorldPop_2020_1km_UNadj_block_count']=wp_counts
    df['uniform_area_block_count']=total*df.land_area_m2.to_numpy()/df.land_area_m2.sum()
    models={'GHSL raw':gh_counts,'WorldPop raw':wp_counts,
            'GHSL normalized to Census total':gh_counts*(total/gh_counts.sum()),
            'WorldPop normalized to Census total':wp_counts*(total/wp_counts.sum()),
            'Uniform over Census land area':df.uniform_area_block_count.to_numpy()}
    rows=[]
    for name,pred in models.items(): rows.append(dict(model=name,normalization='none' if 'raw' in name else ('Census total' if 'normalized' in name else 'Census total; land-area weights'),**metrics(obs,pred)))
    summary=pd.DataFrame(rows)
    df.to_csv(OUT/'fresno_census_block_population_model_counts.csv',index=False)
    summary.to_csv(OUT/'fresno_census_block_population_model_metrics.csv',index=False)
    census_sha={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in BLOCK_DIR.glob('blocks_*.geojson')}
    meta=dict(blocks=len(df),Census_total=total,Census_source_field='2020 DHC P0010001; differential privacy applies',
      geometry_transfer='Raster-cell counts allocated to Census blocks in ESRI:54009 by exact polygon intersection area fraction; assumes uniform within source raster cell. GHSL 1-km cells are equal-area; WorldPop 30-arcsecond cells are transformed from EPSG:4326 to the equal-area projection before intersection.',
      WorldPop_license='CC BY 4.0',WorldPop_doi='10.5258/SOTON/WP00671',
      GHSL_product='GHS-POP E2020 GLOBE R2023A, 1 km, ESRI:54009, V1.0',
      source_sha256={'WorldPop':hashlib.sha256(wp_path.read_bytes()).hexdigest(),'GHSL_tile_zip':hashlib.sha256((INPUT/'ghsl_2020_1km_fresno_tile.zip').read_bytes()).hexdigest(),'Census_block_geojson':census_sha},
      interpretation='Coarse block-total benchmark only. Census, GHSL and WorldPop may share census inputs. Does not independently validate within-block household locations; normalized models isolate relative allocation shape but are calibrated to the Census total.',
      ghsl_tile=tif,worldpop_window_shape=list(wp.shape),models=list(models))
    (OUT/'method.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    print(summary.to_string(index=False))
    print('model totals',gh_counts.sum(),wp_counts.sum(),'pixels in windows',gh.shape,wp.shape)

if __name__=='__main__':main()
