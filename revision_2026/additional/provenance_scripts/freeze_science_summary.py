"""Freeze the six executed analyses, with limitations and actual output hashes."""
from pathlib import Path
import hashlib
import json
import pandas as pd
ROOT=Path(__file__).resolve().parent.parent
TABLES={
 'E7':['real_block_masking/real_block_unmasked_baseline.csv','real_block_masking/real_block_masking_summary.csv'],
 'E8':['spatial_acquisition/spatial_calibration_baselines.csv','spatial_acquisition/spatial_calibration_summary.csv','spatial_acquisition/acquisition_date_influence_summary.csv'],
 'E9':['training_quality/quality_coverage_error.csv'],
 'E10':['new_regions/new_region_support_baselines.csv','new_regions/new_region_allocation_summary.csv'],
 'E11':['residential_allocation/residential_allocation_summary.csv','residential_allocation/residential_landuse_eligibility.csv'],
 'E12':['gnss_transfer/gnss_transfer_summary.csv'],
}
LIMITS={
 'E7':'Artificial masks on published known rates, not truth in genuine no-data cells. Undefined blocks and common-computable denominators retained.',
 'E8':'Held-out recalibration of released population grids, not original population-model training CV. Acquisition-date deletion changes support summaries, not re-inverted velocities; no iid confidence intervals.',
 'E9':'Training-only post-filters on fixed inversion; 44 excluded-edge residuals remain internally dependent, not absolute velocity accuracy.',
 'E10':'Two regions and all 96 scale/origin/quality cases retained; small effects do not establish a universal magnitude. Published averaged coherence is not pair-specific support.',
 'E11':'Municipal land-use eligibility, not household truth; non-2020 vintage and unmapped areas retained. Common domain covers only 89.13-89.22% of complete Census population, failing the initial 90% diagnostic coverage gate.',
 'E12':'Archived independent GNSS comparison on external published California LOS fields, not fresh daily fits or Bangkok validation. Same period and spatial frame, different temporal models; all large negative results retained.',
}
def main():
    ledger=json.loads((ROOT/'completion_ledger.json').read_text(encoding='utf-8'))
    records={}
    for req in ledger['requirements']:
        ident=req['id']
        if ident not in TABLES:continue
        req['status']='experiment_executed_with_limits'
        req['limits']=LIMITS[ident]
        req['evidence']=[]
        records[ident]={}
        for name in TABLES[ident]:
            path=ROOT/'results'/name
            frame=pd.read_csv(path)
            assert len(frame)>0
            req['evidence'].append({'path':'results/'+name,'rows':len(frame),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
            records[ident][name]=json.loads(frame.to_json(orient='records'))
    ledger['completion_proven']=False
    ledger['pending']='Manuscript/SI/responses/figures/locators and new public numerical reproduction; raw source provenance packaging.'
    (ROOT/'completion_ledger.json').write_bytes((json.dumps(ledger,ensure_ascii=False,indent=2)+'\n').encode('utf-8'))
    (ROOT/'scientific_number_ledger.json').write_bytes((json.dumps({'table_records':records,'limits':LIMITS},ensure_ascii=False,indent=2)+'\n').encode('utf-8'))
    print({k:sum(len(v) for v in vals.values()) for k,vals in records.items()})
if __name__=='__main__':main()
