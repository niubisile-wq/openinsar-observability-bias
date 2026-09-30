import os
"""Spatial decomposition of population/support allocation effects.

Uses fixed Bangkok LiCSBAS matched-stack outputs. This is a finite-domain
diagnostic, not cross-validation, a confidence interval, or deformation truth.
"""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from pyproj import Transformer

ROOT=Path(os.environ.get('INSAR_FOLLOWUP_ROOT', Path(__file__).resolve().parent))
SRC=ROOT/'inputs'/'bangkok_matched'/'matched_quality.npz'
OUT=ROOT/'bangkok_spatial_decomposition'
OUT.mkdir(exist_ok=True)
T=Transformer.from_crs('EPSG:4326','EPSG:32647',always_xy=True)

def grouped(label, x, n):
    return np.bincount(label.ravel(),weights=np.asarray(x).ravel(),minlength=n)

def run_grid(xy, size):
    ix=np.floor(xy[0]/size).astype('int64')
    iy=np.floor(xy[1]/size).astype('int64')
    pairs=np.column_stack((ix.ravel(),iy.ravel()))
    unique, inv=np.unique(pairs,axis=0,return_inverse=True)
    label=inv.reshape(ix.shape)
    return label, unique, len(unique)

def main():
    with np.load(SRC) as z:
        x=z['x_edges']; y=z['y_edges']
        pop=z['population'].astype('float64'); area=z['area_m2'].astype('float64')
        fields={'437_pair_mean_input':z['matched_support'].astype('float64'),
                'final_quality_mask':z['full_valid'].astype('float64')}
    lon,lat=np.meshgrid((x[:-1]+x[1:])/2,(y[:-1]+y[1:])/2)
    east,north=T.transform(lon,lat)
    xy=(np.asarray(east),np.asarray(north))
    total=float(pop.sum())
    summaries=[]; by5=[]; by1=[]
    for name,support in fields.items():
        for scale in [1000,5000]:
            label,cells,n=run_grid(xy,scale)
            P=grouped(label,pop,n); A=grouped(label,area,n)
            H=grouped(label,pop*support,n); S=grouped(label,area*support,n)
            area_fraction=np.divide(S,A,out=np.zeros(n),where=A>0)
            pop_fraction=np.divide(H,P,out=np.zeros(n),where=P>0)
            delta=H-P*area_fraction  # uniform deficit minus native deficit
            for i,(ix,iy) in enumerate(cells):
                if P[i] <= 0: continue
                rec=dict(support=name,scale_m=scale,grid_easting_index=int(ix),grid_northing_index=int(iy),
                    origin_easting_m=int(ix*scale),origin_northing_m=int(iy*scale),
                    population_units=P[i],native_support_fraction=pop_fraction[i],
                    area_support_fraction=area_fraction[i],native_missing_support_population=P[i]-H[i],
                    uniform_missing_support_population=P[i]*(1-area_fraction[i]),
                    uniform_minus_native_missing_population=delta[i],
                    contribution_pp_of_domain_population=100*delta[i]/total,
                    native_low_support=bool(pop_fraction[i]<.5),area_low_support=bool(area_fraction[i]<.5))
                if scale==1000: by1.append(rec)
                else: by5.append(rec)
            if scale in (1000,5000):
                pos=float(delta[delta>0].sum()); neg=float(delta[delta<0].sum())
                native_missing=float(P.sum()-H.sum())
                uniform_missing=float(P.sum()-(P*area_fraction).sum())
                summaries.append(dict(support=name,population_units=total,
                    reporting_unit_m=scale,
                    native_missing_population=native_missing,
                    uniform_missing_population=uniform_missing,
                    native_missing_pct=100*native_missing/total,
                    uniform_missing_pct=100*uniform_missing/total,
                    uniform_minus_native_pp=100*float(delta.sum())/total,
                    positive_tile_effect_units=pos,negative_tile_effect_units=neg,
                    cancellation_fraction=1-abs(pos+neg)/(pos+abs(neg)) if pos+abs(neg)>0 else 0,
                    spatial_tiles=n,low_support_threshold=.5,
                    population_flagged_low_by_native=float(P[pop_fraction<.5].sum()),
                    population_flagged_low_by_area=float(P[area_fraction<.5].sum()),
                    population_in_discordant_low_support_units=float(P[(pop_fraction<.5)!=(area_fraction<.5)].sum()),
                    discordant_1km_unit_count=int(np.sum((pop_fraction<.5)!=(area_fraction<.5)))))
    pd.DataFrame(by1).to_csv(OUT/'bangkok_1km_support_local_effects.csv',index=False)
    pd.DataFrame(by5).to_csv(OUT/'bangkok_5km_support_local_effects.csv',index=False)
    pd.DataFrame(summaries).to_csv(OUT/'bangkok_spatial_effect_summary.csv',index=False)
    meta=dict(source='matched_quality.npz',source_sha256=hashlib.sha256(SRC.read_bytes()).hexdigest(),
      method='Fixed projected UTM Zone 47N grids aligned to integer kilometre coordinates; exact sums of supplied per-pixel area, GHSL allocation, and support. At each reporting unit, compare population-weighted support with area-weighted support; 5-km signed contributions decompose the AOI aggregate.',
      support_fields=list(fields),reporting_scales_m=[1000,5000],support_cutoff=.5,
      meanings=dict(native='input GHSL allocation units multiplied by pixel support',uniform='same unit population multiplied by within-unit area-average support',
        low_support='unit support fraction below 50%; exploratory support-category threshold, not a hazard threshold'),
      caveats=['This is a finite-domain decomposition, not spatial cross-validation or an uncertainty interval.','All pixels of the matched LiCSBAS processing grid are retained; it does not impute or establish motion in invalid locations.','Population raster units are allocation values, not household truth.','UTM grid origins are fixed at integer kilometres; other grid origins are represented in the prior factorial sensitivity, not this regional decomposition.'])
    (OUT/'method.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    print(pd.DataFrame(summaries).to_string(index=False))

if __name__=='__main__': main()
