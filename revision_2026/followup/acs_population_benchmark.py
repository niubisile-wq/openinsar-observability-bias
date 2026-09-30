"""Compare Fresno 1-km population grids with independent ACS block-group estimates.

ACS 2017-2021 estimates are survey estimates (90% MOE), not ground truth.
The analysis selects only 2020 block groups wholly contained in the fixed Fresno
rectangle, aggregates 2020 Census block counts and allocates raster-cell counts
by intersection area. It evaluates totals and relative allocation separately.
"""
from pathlib import Path
from zipfile import ZipFile
from collections import defaultdict
import csv, hashlib, json, os
import numpy as np
import rasterio
from openpyxl import load_workbook
from pyproj import Transformer
from rasterio.windows import from_bounds
from shapely.geometry import box, Polygon, shape
from shapely.ops import transform as shp_transform
from shapely.strtree import STRtree
import fiona
from scipy.stats import spearmanr

ROOT=Path(os.environ.get('INSAR_PAIRED_ROOT',Path(__file__).resolve().parent))
INPUT=ROOT/'inputs'/'acs_2021'
POP=ROOT/'inputs'/'fresno_population'
OUT=ROOT/'acs_population_benchmark'
ROI=(-119.9,36.6,-119.7,36.85)
TO_EQUAL=Transformer.from_crs('EPSG:4326','ESRI:54009',always_xy=True)


def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()


def read_acs():
 compact=ROOT/'inputs'/'acs_compact'/'fresno_acs_values.json'
 if compact.is_file() and not (INPUT/'ca_seq0002.zip').exists():
  frozen=json.loads(compact.read_text(encoding='utf-8'))
  values=frozen['acs']
  for value in values.values():value['replicates']=np.asarray(value['replicates'],dtype=np.float64)
  return values,frozen['metadata']
 lookup=load_workbook(INPUT/'lookup.xlsx',read_only=True,data_only=True).active
 match=None
 for row in lookup.iter_rows(min_row=2,values_only=True):
  if row[1]=='B01003' and int(row[2])==2 and int(row[4])==130:
   match=row; break
 lookup.parent.close()
 if match is None: raise ValueError('B01003 sequence and field position do not match the retained lookup')
 # ACS sequence fields: six record identifiers followed by 124 cells; B01003 total is last.
 index=int(match[4])-1
 geo_wb=load_workbook(INPUT/'mini_geo.xlsx',read_only=True,data_only=True)
 geows=geo_wb['ca']; logrec_to_geoid={}
 for state,logrec,geoid,name in geows.iter_rows(min_row=2,values_only=True):
  if geoid and str(geoid).startswith('15000US06019'):
   logrec_to_geoid[str(logrec).zfill(7)]=str(geoid)[7:]
 geo_wb.close()
 values={}
 for prefix in ['e','m']:
  archive=INPUT/'ca_seq0002.zip'
  with ZipFile(archive) as z:
   member=f'{prefix}20215ca0002000.txt'
   with z.open(member) as raw:
    for row in csv.reader((line.decode('utf-8') for line in raw)):
     if len(row)<=index: raise ValueError(f'ACS sequence row too short: {len(row)}')
     logrec=row[5].zfill(7)
     geoid=logrec_to_geoid.get(logrec)
     if geoid:
      values.setdefault(geoid,{})[prefix]=int(row[index]) if row[index] else 0
 # Census variance-replicate file provides estimate, published MOE and 80 replicate estimates.
 vre={}
 with ZipFile(INPUT/'replicate_bg_B01003_06.csv.zip') as z:
  member=z.namelist()[0]
  with z.open(member) as raw:
   reader=csv.DictReader(line.decode('utf-8') for line in raw)
   for row in reader:
    full=str(row.get('GEOID') or '')
    geoid=full[-12:] if full else ''
    if geoid.startswith('06019') and row.get('ESTIMATE'):
     reps=np.array([float(row[f'Var_Rep{i}']) for i in range(1,81)],dtype=np.float64)
     vre[geoid]={'estimate':int(row['ESTIMATE']),'moe':float(row['MOE']),'replicates':reps}
 for geoid,data in values.items():
  if geoid in vre and data.get('e')!=vre[geoid]['estimate']:
   raise ValueError(f'Sequence/VRE estimate mismatch for {geoid}')
  if geoid in vre and data.get('m')!=int(round(vre[geoid]['moe'])):
   raise ValueError(f'Sequence/VRE published MOE mismatch for {geoid}')
 return {g:{**v,**vre[g]} for g,v in values.items() if 'e' in v and 'm' in v and g in vre}, {'sequence':2,'start_position':130,'cell_index_zero_based':index,'table':'B01003','variable':'Total population','moe_confidence':'90 percent','variance_replicates':80,'replicate_factor':4/80,'moe_factor':1.645}


def load_groups():
 acs,acs_meta=read_acs(); shapes={}; props={}
 # TIGER/Line 2020 California block groups; GEOID vintage matches 2020 boundaries.
 tiger_dir=OUT/'tiger_bg_source'; tiger_dir.mkdir(exist_ok=True)
 if (INPUT/'tl_2020_06_bg.zip').is_file():
  with ZipFile(INPUT/'tl_2020_06_bg.zip') as z:
   shpname=next(n for n in z.namelist() if n.lower().endswith('.shp'))
   stem=Path(shpname).stem
   for suffix in ['.shp','.dbf','.shx','.prj']:
    source=f'{stem}{suffix}'
    (tiger_dir/source).write_bytes(z.read(source))
  geography=tiger_dir/f'{stem}.shp'
 else:
  geography=ROOT/'inputs'/'acs_compact'/'fresno_block_groups.geojson'
 with fiona.open(geography) as src:
  for feat in src:
   p=feat['properties']; geoid=str(p.get('GEOID') or '')
   if not geoid.startswith('06019'): continue
   g=shape(feat['geometry']); x0,y0,x1,y1=g.bounds
   if x0<ROI[0]-1e-10 or y0<ROI[1]-1e-10 or x1>ROI[2]+1e-10 or y1>ROI[3]+1e-10: continue
   if geoid not in acs: continue
   shapes[geoid]=shp_transform(TO_EQUAL.transform,g); props[geoid]=p
 if not shapes: raise ValueError('No complete Fresno block groups matched ACS and fixed ROI')
 # Census 2020 block counts summed to matching block groups (existing verified Census P1 block input).
 block_pop=defaultdict(int); block_counts=defaultdict(int)
 # The retained S20 audit table contains 2020 PL94-171 P1 totals for all complete ROI blocks.
 import pandas as pd
 blocks=pd.read_csv(ROOT/'per_block.csv',dtype={'GEOID20':str})
 for _,row in blocks.iterrows():
  geoid=str(row['GEOID20']); bg=geoid[:12]
  if bg in shapes:
   block_pop[bg]+=int(row['census_population_2020']); block_counts[bg]+=1
 missing=set(shapes)-set(block_pop)
 if missing: raise ValueError(f'Census block counts missing for {len(missing)} selected ACS block groups')
 return shapes,props,acs,acs_meta,block_pop,block_counts


def raster_counts(path,geoms,kind,nodata):
 tree=STRtree(list(geoms.values())); geoid_list=list(geoms)
 out={g:0.0 for g in geoid_list}
 with rasterio.open(path) as src:
  bounds=box(*TO_EQUAL.transform_bounds(*ROI,densify_pts=64))
  if str(src.crs)=='EPSG:4326':
   w=from_bounds(*ROI,transform=src.transform).round_offsets().round_lengths()
  else:
   w=from_bounds(*bounds.bounds,transform=src.transform).round_offsets().round_lengths()
  w=w.intersection(rasterio.windows.Window(0,0,src.width,src.height))
  a=src.read(1,window=w); t=src.window_transform(w); crs=str(src.crs); nd=src.nodata if src.nodata is not None else nodata
  pix_area=None
  if crs=='ESRI:54009': pix_area=abs(t.a*t.e)
  for r in range(a.shape[0]):
   for c in np.flatnonzero(np.isfinite(a[r])&(a[r]!=nd)&(a[r]>0)):
    if crs=='ESRI:54009':
     x0=t.c+c*t.a; y1=t.f+r*t.e; cell=box(x0,y1+t.e,x0+t.a,y1)
    elif crs=='EPSG:4326':
     x0=t.c+c*t.a; x1=x0+t.a; y1=t.f+r*t.e; y0=y1+t.e
     xs=np.linspace(min(x0,x1),max(x0,x1),9); ys=np.linspace(min(y0,y1),max(y0,y1),9)
     ring=[(x,min(y0,y1)) for x in xs]+[(max(x0,x1),y) for y in ys[1:]]+[(x,max(y0,y1)) for x in xs[-2::-1]]+[(min(x0,x1),y) for y in ys[-2:0:-1]]
     cell=shp_transform(TO_EQUAL.transform,Polygon(ring)); pix_area=cell.area
    else: raise ValueError(f'Unexpected CRS {crs} in {kind}')
    ix=tree.query(cell,predicate='intersects')
    for i in ix:
     inter=cell.intersection(geoms[geoid_list[int(i)]]).area
     if inter>0: out[geoid_list[int(i)]]+=float(a[r,c])*inter/pix_area
 return out


def metrics(obs,pred,moe,total_moe):
 obs=np.asarray(obs,float); pred=np.asarray(pred,float); moe=np.asarray(moe,float)
 err=pred-obs; total=float(obs.sum()); pt=float(pred.sum())
 corr=float(spearmanr(obs,pred).statistic) if len(obs)>1 else float('nan')
 # Census 90% MOEs for nonoverlapping block groups aggregate by root-sum-square.
 scale=total/pt if pt>0 else float('nan'); norm=pred*scale
 return {'groups':len(obs),'acs_total_estimate':total,'acs_total_90pct_moe_exact_sdr':total_moe,'predicted_total':pt,
         'signed_total_error':pt-total,'signed_total_error_pct_of_acs':100*(pt-total)/total if total else None,
         'total_error_within_exact_acs_90pct_moe':bool(abs(pt-total)<=total_moe),
         'group_spearman':corr,'raw_group_mae':float(np.mean(np.abs(err)),),
         'raw_group_mean_abs_pct_error':float(np.mean(np.abs(err)/np.maximum(obs,1))*100),
         'raw_groups_within_acs_90pct_moe':int(np.sum(np.abs(err)<=moe)),
         'raw_groups_within_acs_moe_fraction':float(np.mean(np.abs(err)<=moe)),
         'acs_total_normalized_scale_factor':scale,
         'normalized_group_mae':float(np.mean(np.abs(norm-obs))),
         'normalized_group_mean_abs_pct_error':float(np.mean(np.abs(norm-obs)/np.maximum(obs,1))*100),
         'normalized_groups_within_acs_90pct_moe':int(np.sum(np.abs(norm-obs)<=moe))}


def main():
 OUT.mkdir(parents=True,exist_ok=True)
 geoms,props,acs,acs_meta,census,nblocks=load_groups()
 ghsl_zip=POP/'ghsl_2020_1km_fresno_tile.zip'
 with ZipFile(ghsl_zip) as z:
  tif=next(n for n in z.namelist() if n.lower().endswith('.tif'))
  ghsl_path=OUT/Path(tif).name
  if not ghsl_path.exists(): ghsl_path.write_bytes(z.read(tif))
 wp_path=POP/'worldpop_2020_1km_us_unadj.tif'
 ghsl=raster_counts(ghsl_path,geoms,'GHSL',-9999)
 world=raster_counts(wp_path,geoms,'WorldPop',-9999)
 ids=sorted(geoms)
 acs_est=np.array([acs[g]['e'] for g in ids],float); acs_moe=np.array([acs[g]['m'] for g in ids],float)
 replicate_totals=np.sum(np.stack([acs[g]['replicates'] for g in ids]),axis=0)
 exact_variance=(4.0/80.0)*float(np.square(replicate_totals-acs_est.sum()).sum())
 exact_total_moe90=1.645*np.sqrt(exact_variance)
 results=[]
 for name,pred in [('2020 Census P1 block total',np.array([census[g] for g in ids],float)),('GHSL 2020 1-km raw',np.array([ghsl[g] for g in ids])),('WorldPop 2020 1-km UN-adjusted raw',np.array([world[g] for g in ids]))]:
  results.append({'model':name,**metrics(acs_est,pred,acs_moe,exact_total_moe90)})
 rows=[]
 for i,g in enumerate(ids):
  rows.append({'GEOID':g,'name':str(props[g].get('NAMELSAD') or ''),'census_block_count':nblocks[g],'census_2020_P1':census[g],
    'ACS_2017_2021_B01003_estimate':int(acs[g]['e']),'ACS_2017_2021_B01003_MOE_90pct':int(acs[g]['m']),
    'GHSL_2020_1km_count_allocated_by_area':ghsl[g],'WorldPop_2020_1km_UNadjusted_count_allocated_by_area':world[g]})
 import pandas as pd
 pd.DataFrame(rows).to_csv(OUT/'blockgroup_counts.csv',index=False)
 pd.DataFrame(results).to_csv(OUT/'model_summary.csv',index=False)
 input_files=[INPUT/'lookup.xlsx',INPUT/'mini_geo.xlsx',INPUT/'ca_seq0002.zip',INPUT/'tl_2020_06_bg.zip',ghsl_zip,wp_path]
 method={'study_roi_lonlat':ROI,'selected_complete_block_groups':len(ids),'definition':'2020 Census block-group polygons whose full bounding boxes lie inside the predeclared Fresno rectangle; exact geometries overlaid in ESRI:54009.',
  'acs_source':'2017-2021 ACS 5-year Summary File, California sequence 0002, Table B01003 total population estimate and 90% margin of error; geography LOGRECNO joined through California 5-year Mini Geo.',
  'acs_mapping':acs_meta,'moe_aggregation':'Exact SDR variance replicate aggregation: sum each of 80 B01003 block-group replicate estimates over the selected disjoint block groups; variance=(4/80)*sum((replicate_sum-point_estimate_sum)^2); 90% MOE=1.645*sqrt(variance).',
  'raster_overlay':'Source raster population counts allocated to block-group polygons by exact intersection area fraction, assuming uniform population within each 1-km source pixel. GHSL ESRI:54009; WorldPop EPSG:4326 cell edges densified before equal-area projection.',
  'temporal_note':'ACS 2017-2021 five-year estimates overlap 2020 but are a period estimate; GHSL/WorldPop are 2020 surfaces and Census P1 is the 2020 decennial count.',
  'limitations':['Block groups are much coarser than the raster and cannot validate within-block-group placement.','ACS is a survey estimate with MOE, not a gold standard; 2020 Census and ACS estimates have different error structures.','Source lineage independence between GHSL/WorldPop and ACS is not established.','The selected complete block groups form a subset of the Fresno rectangle, not the full city/county.','Within-pixel area allocation is an explicit modeling assumption.'],
  'inputs':{p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in [*input_files,INPUT/'replicate_bg_B01003_06.csv.zip',INPUT/'replicate_table_list.csv',*(ROOT/'inputs'/'acs_compact').glob('*')] if p.is_file()}}
 (OUT/'method.json').write_text(json.dumps(method,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
 print(pd.DataFrame(results).to_string(index=False)); print('selected BGs',len(ids),'ACS exact aggregate 90% MOE',exact_total_moe90)

if __name__=='__main__': main()

