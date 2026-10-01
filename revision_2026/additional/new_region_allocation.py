from pathlib import Path
import argparse
import json
import hashlib
import numpy as np
import pandas as pd
import shapely
from shapely.strtree import STRtree

def analyse(inputs,out):
    out.mkdir(parents=True,exist_ok=True);summaries=[];local=[];baselines=[]
    for region in ['Los_Angeles','Davis_Sacramento']:
        z=np.load(inputs/(region+'.npz'),allow_pickle=False)
        p=z['population'];a=z['area_m2'];cells=shapely.from_wkb(z['geometry_wkb'])
        density=p/a;total=float(p.sum())
        bounds=shapely.total_bounds(cells);tree=STRtree(cells)
        for q in [20,30,40]:
            s=z[f'support_q{q}'];native=float(np.sum(p*(1-s)))
            baselines.append({'region':region,'coherence_threshold':q/100,'population_units':total,
                              'area_support_share_pct':100*np.sum(a*s)/a.sum(),
                              'population_support_share_pct':100*(total-native)/total,
                              'native_deficit_allocation_units':native})
        for size in [250,500,1000,2000]:
            for ox,oy in [(0,0),(size/2,0),(0,size/2),(size/2,size/2)]:
                gx=np.arange(int(np.floor((bounds[0]+ox)/size)),int(np.ceil((bounds[2]+ox)/size)))
                gy=np.arange(int(np.floor((bounds[1]+oy)/size)),int(np.ceil((bounds[3]+oy)/size)))
                xx,yy=np.meshgrid(gx,gy);xx=xx.ravel();yy=yy.ravel()
                units=shapely.box(xx*size-ox,yy*size-oy,(xx+1)*size-ox,(yy+1)*size-oy)
                ti,ci=tree.query(units,predicate='intersects')
                area=shapely.area(shapely.intersection(units[ti],cells[ci]));ok=area>1e-6;ti=ti[ok];ci=ci[ok];area=area[ok]
                n=len(units);aa=np.bincount(ti,weights=area,minlength=n)
                pp=np.bincount(ti,weights=area*density[ci],minlength=n);positive=aa>0
                assert abs(pp.sum()-total)<1e-5
                assert abs(aa.sum()-a.sum())/a.sum()<1e-9
                for q in [20,30,40]:
                    s=z[f'support_q{q}'];sa=np.bincount(ti,weights=area*s[ci],minlength=n)
                    native=np.bincount(ti,weights=area*density[ci]*(1-s[ci]),minlength=n)
                    uniform=pp*(1-np.divide(sa,aa,out=np.zeros(n),where=aa>0));delta=uniform-native
                    mr=np.divide(pp,aa,out=np.zeros(n),where=aa>0);ms=np.divide(sa,aa,out=np.zeros(n),where=aa>0)
                    second_r=np.divide(np.bincount(ti,weights=area*density[ci]**2,minlength=n),aa,out=np.zeros(n),where=aa>0)
                    second_s=np.divide(np.bincount(ti,weights=area*s[ci]**2,minlength=n),aa,out=np.zeros(n),where=aa>0)
                    bound=aa*np.sqrt(np.maximum(0,second_r-mr**2)*np.maximum(0,second_s-ms**2))
                    assert np.all(np.abs(delta)<=bound+1e-6)
                    n_top=max(1,int(np.ceil(.1*positive.sum())));active=np.flatnonzero(positive)
                    topnative=set(active[np.argsort(native[active],kind='stable')[-n_top:]])
                    topuniform=set(active[np.argsort(uniform[active],kind='stable')[-n_top:]])
                    summary={'region':region,'coherence_threshold':q/100,'scale_m':size,'origin_x_m':ox,'origin_y_m':oy,
                             'reporting_cells':int(positive.sum()),'population_units':total,'native_deficit_units':float(native.sum()),
                             'uniform_deficit_units':float(uniform.sum()),'uniform_minus_native_pp':100*delta.sum()/total,
                             'local_L1_error_pp':100*np.abs(delta).sum()/total,'sum_local_Cauchy_bound_pp':100*bound.sum()/total,
                             'top_decile_overlap':len(topnative&topuniform)/n_top,
                             'positive_error_cells':int(np.count_nonzero(delta>1e-8)),'negative_error_cells':int(np.count_nonzero(delta< -1e-8))}
                    summaries.append(summary)
                    for i in active:
                        local.append({'region':region,'coherence_threshold':q/100,'scale_m':size,'origin_x_m':ox,'origin_y_m':oy,
                                      'tile_x':int(xx[i]),'tile_y':int(yy[i]),'population_units':pp[i],
                                      'native_deficit_units':native[i],'uniform_deficit_units':uniform[i],
                                      'signed_error_units':delta[i],'Cauchy_bound_units':bound[i]})
        print('New region analysed',region,flush=True)
    frame=pd.DataFrame(summaries);assert len(frame)==96
    frame.to_csv(out/'new_region_allocation_summary.csv',index=False)
    pd.DataFrame(local).to_csv(out/'new_region_allocation_by_tile.csv',index=False)
    pd.DataFrame(baselines).to_csv(out/'new_region_support_baselines.csv',index=False)
    protocol={'regions':2,'summary_designs':96,'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(inputs.glob('*.npz'))},
              'checks':{'area_and_population_conserved_per_partition':True,'per_tile_Cauchy_bound_verified':True},
              'limits':['Transfer of an allocation/support identity to two published InSAR product regions; no new velocity truth is supplied.',
                        'Fixed region and threshold/scale/origin designs retain positive, negative and small effects.',
                        'Top-decile overlaps describe the fixed allocation/support deficit rankings, not physical risk priorities.']}
    (out/'new_region_allocation_protocol.json').write_bytes((json.dumps(protocol,indent=2)+'\n').encode())
    print(frame.groupby(['region','coherence_threshold']).uniform_minus_native_pp.agg(['min','max']).to_string(),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path);p.add_argument('--output',type=Path);args=p.parse_args();root=Path(__file__).resolve().parent.parent
    analyse(args.inputs or root/'inputs/new_regions',args.output or root/'results/new_regions')
