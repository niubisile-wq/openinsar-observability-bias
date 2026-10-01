"""Independent municipal land-use overlap on fixed official Census geometry."""
from pathlib import Path
import os
import sys
import hashlib
import json
import numpy as np
import shapely
from shapely.strtree import STRtree
from pyproj import Transformer
ROOT=Path(__file__).resolve().parent.parent
REV=ROOT.parent/'github_original_source/revision_2026'
STUDY=ROOT.parent/'01_当前返修工作区/InSAR_revision_20260926/strengthening'
OUT=ROOT/'inputs/residential_reference'
os.environ['INSAR_STRENGTHENING_ROOT']=str(STUDY)
os.environ['INSAR_FOLLOWUP_ROOT']=str(ROOT/'intermediates')
sys.path.insert(0,str(REV/'followup'))
from fresno_capop_support_sensitivity import census_counts,complete_blocks,target_grid

# Semantic category definitions are frozen before reading any coverage or allocation outcomes.
STRICT={'rl','rml','rm','rmh','rh','rr','rmht','run'}
BROAD=STRICT|{'raur','rmx','cmx','cce'}

def main():
    protocol={'strict_residential_codes':sorted(STRICT),'broader_eligibility_codes':sorted(BROAD),
              'code_choice':'Existing residential categories plus mobile-home parks; broader sensitivity adds explicitly mixed-use/residential-agricultural categories.',
              'target':'Residential land-use eligibility, not resident counts, households or demographic locations.'}
    (OUT/'category_protocol.json').write_bytes((json.dumps(protocol,indent=2)+'\n').encode())
    features=[]
    for path in sorted(OUT.glob('landuse_*.geojson')):features.extend(json.loads(path.read_text(encoding='utf-8'))['features'])
    assert len(features)==123493
    geoms=shapely.from_geojson(np.array([json.dumps(f['geometry']) for f in features]))
    projection=Transformer.from_crs('EPSG:4326','EPSG:3310',always_xy=True)
    geoms=shapely.transform(geoms,lambda xy:np.column_stack(projection.transform(xy[:,0],xy[:,1])))
    invalid=~shapely.is_valid(geoms);invalid_n=int(invalid.sum())
    geoms[invalid]=shapely.make_valid(geoms[invalid])
    codes=np.array([str(f['properties'].get('ELU') or '').lower() for f in features])
    classified=codes!='';strict=np.isin(codes,list(STRICT));broad=np.isin(codes,list(BROAD))
    # FeatureRefreshDate is a layer refresh time, not the actual vintage of every classification.
    refresh=sorted(set(str(f['properties'].get('LayerRefreshDate')) for f in features))
    blocks,records=complete_blocks(census_counts());b=np.asarray(blocks,dtype=object)
    cap,ghs,are,obsa,tr,warp=target_grid(blocks)
    h,w=cap.shape;yy,xx=np.indices(cap.shape)
    cx=tr.c+(xx.ravel()+.5)*tr.a;cy=tr.f+(yy.ravel()+.5)*tr.e
    cells=shapely.box(cx-tr.a/2,cy+tr.e/2,cx+tr.a/2,cy-tr.e/2)
    ci,bi=STRtree(b).query(cells,predicate='intersects')
    clips=shapely.intersection(cells[ci],b[bi]);area=shapely.area(clips)
    keep=area>1e-6;ci=ci[keep];bi=bi[keep];clips=clips[keep];area=area[keep]
    fractions=np.zeros((len(clips),3));tree=STRtree(geoms)
    for offset in range(0,len(clips),5000):
        cc,pi=tree.query(clips[offset:offset+5000],predicate='intersects')
        fragments=shapely.intersection(clips[offset:offset+5000][cc],geoms[pi])
        good=shapely.area(fragments)>1e-7;cc=cc[good];pi=pi[good];fragments=fragments[good]
        order=np.argsort(cc,kind='stable');cc=cc[order];pi=pi[order];fragments=fragments[order]
        unique,starts,counts=np.unique(cc,return_index=True,return_counts=True)
        for index,start,count in zip(unique,starts,counts):
            ix=slice(start,start+count);pieces=fragments[ix];src_ids=pi[ix]
            for column,mask in enumerate([classified,strict,broad]):
                eligible=mask[src_ids]
                if eligible.any():fractions[offset+index,column]=shapely.area(shapely.union_all(pieces[eligible]))/area[offset+index]
        print('Residential exact overlaps',min(offset+5000,len(clips)),'/',len(clips),flush=True)
    assert fractions.max()<=1+1e-6
    fractions=np.clip(fractions,0,1)
    assert np.all(fractions[:,1]<=fractions[:,2]+1e-8) and np.all(fractions[:,2]<=fractions[:,0]+1e-8)
    cell_area=tr.a*(-tr.e);p=np.array([v for _,v in records],dtype=float)
    support=np.clip(np.divide(obsa,are,out=np.zeros_like(obsa),where=are>0),0,1).ravel()[ci]
    np.savez_compressed(OUT/'residential_allocation_inputs.npz',
                        block_geoid=np.array([v[0] for v in records]),block_population=p,
                        intersection_block=bi,intersection_area_m2=area,
                        ghsl_mass=ghs.ravel()[ci]*area/cell_area,capop_mass=cap.ravel()[ci]*area/cell_area,
                        classified_fraction=fractions[:,0],strict_residential_fraction=fractions[:,1],
                        broad_residential_fraction=fractions[:,2],dwr_support_fraction=support)
    provenance={'source':'https://services2.arcgis.com/WkBUojyNPhsWOk1W/arcgis/rest/services/Existing_Land_Use/FeatureServer/19',
                'source_owner':'FresnoGIS','source_item':'bc57c53874e14b1e850d2324b7a568ca',
                'raw_features':len(features),'invalid_geometries_repaired':invalid_n,'layer_refresh_values':refresh,
                'union_method':'Union classified/residential parcel fragments separately inside each exact block/native-cell intersection, avoiding overlap double counting.',
                'independent_dimension':'Municipal residential-use classification distinct from satellite-built-up eligibility; no reference population counts or household occupancy is supplied.',
                'warp_checks':warp,'exact_intersections':len(clips),
                'source_sha256':{path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(OUT.glob('landuse_*.geojson'))},
                'limits':['Unknown actual land-use vintage for each parcel; layer refresh and 2024 item modification are not a 2020 ground-truth date.',
                          'Shared Census totals remain a constraint; population-model source independence is not asserted.',
                          'Commercial institutional/group-quarter residence and mixed use complicate residential eligibility.',
                          'DWR support is held constant inside each 100m projected population cell after conservative area reprojection; subcell household support is not identified.']}
    (OUT/'residential_geometry_provenance.json').write_bytes((json.dumps(provenance,indent=2)+'\n').encode())
    print('Prepared independent residential classification',len(features),invalid_n,'repairs',flush=True)

if __name__=='__main__':main()
