"""Descriptive queen-contiguity Moran's I for Fresno per-block support deficits."""
from pathlib import Path
import os
import json
import numpy as np
import pandas as pd
import shapely
from shapely.strtree import STRtree
try:
    from capop_support_sensitivity import census_counts, complete_blocks
except ModuleNotFoundError:
    from fresno_capop_support_sensitivity import census_counts, complete_blocks

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
RESULTS=Path(os.environ.get('INSAR_FOLLOWUP_ROOT', HERE))/'fresno_population_allocation'
OUT=RESULTS/'spatial_autocorrelation.json'

def moran(x,i,j,n):
    x=np.asarray(x,dtype='float64')
    z=x-x.mean()
    s0=2*len(i)
    numerator=2*float(np.sum(z[i]*z[j]))
    denominator=float(np.sum(z*z))
    if s0<=0 or denominator<=0: raise RuntimeError('Moran I is undefined for this graph/field.')
    return n/s0*numerator/denominator

def main():
    RESULTS.mkdir(parents=True,exist_ok=True)
    counts=census_counts(); geoms,records=complete_blocks(counts)
    frame=pd.read_csv(RESULTS/'per_block.csv',dtype={'GEOID20':str}).set_index('GEOID20')
    indices=[k for k,(gid,pop) in enumerate(records) if bool(frame.loc[gid,'included_common_population_domain'])]
    gids=[records[k][0] for k in indices]
    polygons=[geoms[k] for k in indices]
    population=np.array([records[k][1] for k in indices],dtype='float64')
    tree=STRtree(polygons)
    pairs=tree.query(polygons,predicate='touches')
    i,j=pairs; keep=i<j; i=i[keep]; j=j[keep]
    n=len(polygons); degree=np.bincount(np.concatenate([i,j]),minlength=n)
    ghsl=np.array([frame.loc[gid,'support_deficit_GHSL_normalized_within_block_units'] for gid in gids])
    capop=np.array([frame.loc[gid,'support_deficit_CA_POP_normalized_within_block_units'] for gid in gids])
    values={
      'GHSL_unsupported_fraction':moran(ghsl/population,i,j,n),
      'CA_POP_unsupported_fraction':moran(capop/population,i,j,n),
      'CA_POP_minus_GHSL_unsupported_fraction':moran((capop-ghsl)/population,i,j,n)}
    result={
      'method':'Global Moran I with binary queen-contiguity weights; queen adjacency includes shared edges or vertices; each unordered pair counted once.',
      'domain':'Fresno fixed ROI; common population-allocation domain only; complete 2020 Census blocks.',
      'n_blocks':n,'undirected_neighbor_pairs':int(len(i)),'isolates':int((degree==0).sum()),
      'degree_quantiles':dict(zip(['min','q25','median','q75','q95','max'],[float(v) for v in np.quantile(degree,[0,.25,.5,.75,.95,1])])),
      'fields':values,
      'interpretation':'Descriptive spatial dependence for this fixed domain. No permutation p-value, spatial-sampling confidence interval, or extrapolation to other regions is reported.'}
    if n!=5713 or len(i)!=17083 or (degree==0).sum()!=12: raise RuntimeError('Unexpected queen-neighbor graph; inspect geometries.')
    OUT.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
if __name__=='__main__': main()
