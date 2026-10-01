import os
from pathlib import Path
import json, hashlib
import numpy as np, pandas as pd, rasterio
from rasterio.features import rasterize
from rasterio.transform import Affine
from dwr_units import feet_to_feet_per_year,WINDOWS,DAYS_PER_YEAR
base=Path(os.environ.get('INSAR_FOLLOWUP_ROOT', Path(__file__).resolve().parent))
proj=Path(os.environ['INSAR_STRENGTHENING_ROOT'])
files=list((base/'census_blocks'/'fresno_roi').glob('*.geojson'))
features=[]
for p in files: features.extend(json.loads(p.read_text(encoding='utf-8'))['features'])
blocks=[]
roi=(-119.9,36.6,-119.7,36.85)
for f in features:
 pts=[pt for poly in ([f['geometry']['coordinates']] if f['geometry']['type']=='Polygon' else f['geometry']['coordinates']) for ring in poly for pt in ring]
 xs=[q[0] for q in pts]; ys=[q[1] for q in pts]
 if min(xs)>=roi[0]-1e-10 and max(xs)<=roi[2]+1e-10 and min(ys)>=roi[1]-1e-10 and max(ys)<=roi[3]+1e-10: blocks.append(f)
assert len(blocks)==6788
geom=json.loads((proj/'external'/'dwr'/'roi_features.json').read_text(encoding='utf-8'))['features']
periods=[('2015-10_2016-10','rate_2015_2016.tif'),('2016-10_2017-10','annual_oct_2016_2017.tif'),('2017-10_2018-10','annual_oct_2017_2018.tif'),('2018-10_2019-10','annual_oct_2018_2019.tif'),('2019-10_2020-10','annual_oct_2019_2020.tif'),('2020-10_2021-10','rate_2020_2021.tif')]
ref=base/'dwr_velocity'/periods[0][1]
with rasterio.open(ref) as ds: trans=ds.transform; shape=(ds.height,ds.width); crs=ds.crs
bids=rasterize([(f['geometry'],i+1) for i,f in enumerate(blocks)],out_shape=shape,transform=trans,fill=0,dtype='int32')
support=rasterize([({'type':'Polygon','coordinates':f['geometry']['rings']},1) for f in geom],out_shape=shape,transform=trans,fill=0,dtype='uint8').astype(bool)
assert int(support.sum())==len(geom)==38738
rows=[];allclasses=[]
for label,name in periods:
 with rasterio.open(base/'dwr_velocity'/name) as ds:
  assert ds.transform==trans and (ds.height,ds.width)==shape and ds.crs==crs
  a=ds.read(1).astype(float); valid=support & np.isfinite(a)&(a!=-9999)&(bids>0)
  idx=np.where(valid); bidx=bids[idx]-1; v=feet_to_feet_per_year(a[idx],name)*304.8
  grouped=[[] for _ in blocks]
  for i,x in zip(bidx,v):grouped[int(i)].append(float(x))
  med=np.array([np.median(x) if x else np.nan for x in grouped]); nobs=np.array([len(x) for x in grouped])
  rows.append(pd.DataFrame({'GEOID':[x['properties']['GEOID'] for x in blocks],'census_population_2020':[int(x['properties']['P0010001'] or 0) for x in blocks],'period':label,'valid_cells':nobs,'median_rate_mm_yr':med}))
  allclasses.append(med)
df=pd.concat(rows,ignore_index=True)
wide=df.pivot(index='GEOID',columns='period',values='median_rate_mm_yr')
pop=df.drop_duplicates('GEOID').set_index('GEOID').census_population_2020
common=wide.dropna().join(pop.rename('population'))
class_order=['>=0','[-3,0)','[-5,-3)','[-10,-5)','<-10']
def cl(v):return np.select([v>=0,v>=-3,v>=-5,v>=-10],class_order[:4],default='<-10')
classes=pd.DataFrame(index=common.index)
for p,_ in periods:classes[p]=cl(common[p].values)
classes['population']=common.population
persist=np.sum(np.column_stack([common[p].values < -5 for p,_ in periods]),axis=1)
classes['years_below_minus5']=persist
classes['same_class_all_six']=np.all(classes[[p for p,_ in periods]].values==classes[[periods[0][0]]].values,axis=1)
classes['population_mode_class']=classes[[p for p,_ in periods]].mode(axis=1)[0]
N=classes.population.sum()
rows2=[]
for k in range(7):
 m=classes.years_below_minus5.eq(k);rows2.append({'annual_windows_below_minus5':k,'blocks':int(m.sum()),'population_2020':int(classes.loc[m,'population'].sum()),'population_share_common_pct':100*classes.loc[m,'population'].sum()/N})
pd.DataFrame(rows2).to_csv(base/'census_block_velocity_analysis'/'fresno_six_year_below5_persistence.csv',index=False)
# Full per-block classifications and intermediate annual block medians.
classes.to_csv(base/'census_block_velocity_analysis'/'fresno_six_annual_block_classifications.csv')
df.to_csv(base/'census_block_velocity_analysis'/'fresno_six_annual_block_medians.csv',index=False)
summary={'blocks_total':len(blocks),'blocks_valid_in_all_six_annual_windows':len(common),'population_2020_in_all_six_valid_blocks':int(N),'population_share_of_all_complete_blocks_pct':float(N/sum(int(x['properties']['P0010001'] or 0) for x in blocks)*100),'population_share_with_all_six_classes_identical_pct':float(classes.loc[classes.same_class_all_six,'population'].sum()/N*100),'population_share_all_six_annual_medians_below_minus5_pct':float(classes.loc[classes.years_below_minus5.eq(6),'population'].sum()/N*100),'population_share_never_below_minus5_pct':float(classes.loc[classes.years_below_minus5.eq(0),'population'].sum()/N*100),'population_share_crossing_minus5_at_least_once_pct':float(classes.loc[classes.years_below_minus5.between(1,5),'population'].sum()/N*100),'annual_windows':[x[0] for x in periods],'rate_class_thresholds_mm_yr':class_order,'support_cells':len(geom),'DWR_geojson_sha256':hashlib.sha256((proj/'external'/'dwr'/'roi_features.json').read_bytes()).hexdigest(),'DWR_vertical_rates_sha256':{name:hashlib.sha256((base/'dwr_velocity'/name).read_bytes()).hexdigest() for _,name in periods},'limitations':['Only blocks with valid DWR cells in all six annual windows are included.','Annual-window block medians are classified using fixed 2020 Census counts.','This describes temporal variability in a product, not independent validation or persistent hazard.','DWR rate cells are interpolated vertical products, not manuscript LOS.']}
summary['input_windows']={name:WINDOWS[name] for _,name in periods}
summary['annualization_days_per_year']=DAYS_PER_YEAR
(base/'census_block_velocity_analysis'/'fresno_six_year_stability_method.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2));print(pd.DataFrame(rows2).to_string(index=False));print('\nmedian class by annual window population shares')
for p,_ in periods:
 clp=cl(common[p].values)
 shares={key:float(classes.loc[clp==key,'population'].sum()/N*100) for key in class_order}
 print(p,shares)
