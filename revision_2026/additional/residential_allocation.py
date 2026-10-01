from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd

def analyse(inputs,out):
    out.mkdir(parents=True,exist_ok=True)
    z=np.load(inputs);p=z['block_population'];nb=len(p);bi=z['intersection_block'];area=z['intersection_area_m2']
    gh=z['ghsl_mass'];ca=z['capop_mass'];s=z['dwr_support_fraction'];classified=z['classified_fraction']
    masks={'strict':z['strict_residential_fraction'],'broader_mixed_use':z['broad_residential_fraction']}
    bsum=lambda a:np.bincount(bi,weights=a,minlength=nb)
    totals={'Area_uniform':bsum(area),'GHSL_100m':bsum(gh),'CA_POP_100m':bsum(ca)}
    sums=[];byblock=[];eligibility=[]
    for definition,res in masks.items():
        resarea=area*res;resden=bsum(resarea)
        common=(p>0)&(totals['GHSL_100m']>0)&(totals['CA_POP_100m']>0)&(resden>0)
        full_p=float(p.sum());common_p=float(p[common].sum())
        coverage_gate_met=common_p>=.9*full_p
        weights={'Area_uniform':area,'GHSL_100m':gh,'CA_POP_100m':ca,
                 'Municipal_residential_uniform':resarea,
                 'GHSL_restricted_to_municipal_residential':gh*res,
                 'CA_POP_restricted_to_municipal_residential':ca*res}
        # If a restricted model has zero mass, report it and shrink the fair joint domain.
        for raw in weights.values():common&=bsum(raw)>0
        common_p=float(p[common].sum())
        for name,raw in weights.items():
            den=bsum(raw);pair_ok=common[bi]
            alloc=np.zeros(len(bi));alloc[pair_ok]=p[bi[pair_ok]]*raw[pair_ok]/den[bi[pair_ok]]
            assert abs(alloc.sum()-common_p)<1e-6
            deficit=bsum(alloc*(1-s))
            sums.append({'residential_definition':definition,'allocation_rule':name,'complete_Census_population':full_p,
                         'common_Census_population':common_p,'excluded_Census_population':full_p-common_p,
                         'common_population_share_pct':100*common_p/full_p,'common_blocks':int(common.sum()),
                         'initial_90pct_coverage_gate_met':coverage_gate_met,
                         'reallocated_Census_units':float(alloc.sum()),'support_deficit_units':float(deficit.sum()),
                         'support_deficit_share_pct':100*deficit.sum()/common_p})
            for i in np.flatnonzero(common):
                byblock.append({'residential_definition':definition,'allocation_rule':name,'GEOID20':str(z['block_geoid'][i]),
                                'Census_population_units':p[i],'support_deficit_units':deficit[i]})
        # Measure land-use eligibility on the mapped part; unmapped area is not classified as nonresidential.
        for name,raw,den in [('Area_uniform',area,totals['Area_uniform']),('GHSL_100m',gh,totals['GHSL_100m']),('CA_POP_100m',ca,totals['CA_POP_100m'])]:
            eligible=(p>0)&(den>0);ok=eligible[bi]
            alloc=np.zeros(len(bi));alloc[ok]=p[bi[ok]]*raw[ok]/den[bi[ok]]
            total=float(alloc.sum());inside=float(np.sum(alloc*res));nonres=float(np.sum(alloc*np.maximum(0,classified-res)))
            unmapped=float(np.sum(alloc*(1-classified)))
            assert abs(inside+nonres+unmapped-total)<1e-5
            eligibility.append({'residential_definition':definition,'population_allocation_model':name,
                                'Census_units_denominator':total,'inside_residential_allocation_units':inside,
                                'mapped_nonresidential_allocation_units':nonres,'unmapped_classification_allocation_units':unmapped,
                                'residential_share_of_all_units_pct':100*inside/total,
                                'mapped_nonresidential_share_of_all_units_pct':100*nonres/total,
                                'unmapped_share_of_all_units_pct':100*unmapped/total,
                                'residential_share_within_mapped_units_pct':100*inside/(inside+nonres)})
    pd.DataFrame(sums).to_csv(out/'residential_allocation_summary.csv',index=False)
    pd.DataFrame(byblock).to_csv(out/'residential_allocation_by_block.csv',index=False)
    pd.DataFrame(eligibility).to_csv(out/'residential_landuse_eligibility.csv',index=False)
    protocol={'input_sha256':hashlib.sha256(inputs.read_bytes()).hexdigest(),'definitions':2,'allocation_rules':6,
              'common_domain_rule':'Positive Census/GHSL/CA-POP and residential reference mass for every compared allocation; report excluded units.',
              'checks':{'Census_totals_conserved_per_common_block':True,'unmapped_classification_separate_from_nonresidential':True},
              'limits':['Independent evidence is land-use classification, not independent resident counts or household placement.',
                        'Eligibility disagreement is not a false-population count: mixed use, group quarters, grid mixing and vintage mismatch can contribute.',
                        'No residential rule is presented as a universal replacement for satellite population models.',
                        'Spatially averaged DWR support inside 100m cells cannot resolve house-level coverage.']}
    (out/'residential_allocation_protocol.json').write_bytes((json.dumps(protocol,indent=2)+'\n').encode())
    print(pd.DataFrame(sums).to_string(index=False),flush=True)
    print(pd.DataFrame(eligibility).to_string(index=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path);p.add_argument('--output',type=Path);args=p.parse_args();root=Path(__file__).resolve().parent.parent
    analyse(args.inputs or root/'inputs/residential_reference/residential_allocation_inputs.npz',args.output or root/'results/residential_allocation')
