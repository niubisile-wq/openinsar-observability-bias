import os
from pathlib import Path
import pandas as pd, json, hashlib
base=Path(os.environ.get('INSAR_FOLLOWUP_ROOT', Path(__file__).resolve().parent))
analysis=base/'census_block_velocity_analysis'
src=analysis/'fresno_census_block_velocity_by_block.csv'
df=pd.read_csv(src)
a='2015-10 to 2021-10 temporal median of annual rates'; b='2025-04 to 2026-04 annual'
x=df[df.period.eq(a)].set_index('GEOID'); y=df[df.period.eq(b)].set_index('GEOID')
cols=['census_population_2020','valid_DWR_rate_pixels','median_vertical_rate_mm_yr','rate_class']
c=x[cols].join(y[cols],how='inner',lsuffix='_2015_21',rsuffix='_2025_26')
c=c[(c.valid_DWR_rate_pixels_2015_21>0)&(c.valid_DWR_rate_pixels_2025_26>0)].copy()
def classify(v):
 return pd.cut(v,[-float('inf'),-10,-5,-3,0,float('inf')],labels=['<-10','[-10,-5)','[-5,-3)','[-3,0)','>=0'],right=False).astype(str)
c['class_2015_21']=classify(c.median_vertical_rate_mm_yr_2015_21)
c['class_2025_26']=classify(c.median_vertical_rate_mm_yr_2025_26)
pop=c.census_population_2020_2015_21.astype(float); n=float(pop.sum())
severe_old=c.class_2015_21.isin(['<-10','[-10,-5)']); severe_new=c.class_2025_26.isin(['<-10','[-10,-5)'])
metrics=[('Common blocks with valid rates in both periods',len(c)),('2020 Census population in common blocks',int(n)),('Population share below -5 mm/yr, 2015-2021 (%)',100*pop[severe_old].sum()/n),('Population share below -5 mm/yr, 2025-2026 (%)',100*pop[severe_new].sum()/n),('Population share changing class (%)',100*pop[c.class_2015_21.ne(c.class_2025_26)].sum()/n),('Population share below -5 in both periods (%)',100*pop[severe_old&severe_new].sum()/n),('Population share switching from >=-5 to <-5 (%)',100*pop[(~severe_old)&severe_new].sum()/n),('Population share switching from <-5 to >=-5 (%)',100*pop[severe_old&(~severe_new)].sum()/n)]
pd.DataFrame(metrics,columns=['metric','value']).to_csv(base/'census_block_velocity_analysis'/'fresno_population_temporal_stability_summary.csv',index=False)
order=['>=0','[-3,0)','[-5,-3)','[-10,-5)','<-10']
mat=pd.crosstab(c.class_2015_21,c.class_2025_26,values=pop,aggfunc='sum').fillna(0).reindex(index=order,columns=order,fill_value=0)
mat.index.name='2015-2021 block median class';mat.columns.name='2025-2026 block median class'
mat.to_csv(base/'census_block_velocity_analysis'/'fresno_population_rate_class_transition_2015_21_to_2025_26.csv')
c.reset_index().to_csv(base/'census_block_velocity_analysis'/'fresno_block_temporal_class_transitions.csv',index=False)
meta={'estimand':'population-weighted block rate-class transition using the same complete Census blocks with valid DWR cells in both selected periods','early_period':a,'late_period':b,'population_vintage':'2020 Census P0010001 held fixed in both periods to isolate rate-class changes; late-period values are not a contemporaneous exposure estimate','thresholds_mm_per_year':['>=0','[-3,0)','[-5,-3)','[-10,-5)','<-10'],'common_blocks':int(len(c)),'common_population_2020':int(n),'input_csv':str(src),'input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'limitations':['Different temporal estimands (six-year pixelwise median vs single annual raster).','The 2025-2026 rate is a single annual product and must not be interpreted as long-term trend or persistent hazard.','DWR annual rates are interpolated vertical products, not the manuscript LOS field.','Census 2020 population is intentionally held fixed; it is not a population estimate for 2025-2026.','This is an exploratory sensitivity that demonstrates temporal class instability, not independent validation.']}
(base/'census_block_velocity_analysis'/'fresno_temporal_stability_method.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
print(pd.DataFrame(metrics,columns=['metric','value']).to_string(index=False));print('\n'+mat.to_string())
