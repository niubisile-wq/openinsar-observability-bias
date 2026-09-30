"""E0 factorial audit separating encoding, registration and validity rules."""
import itertools,json,hashlib,requests
import numpy as np,pandas as pd,rasterio,tifffile
from common import PROJECT,BACKUP,ROOT,OUT,dump,sha256,edges,separable_resample

DEST=OUT/'lineage';DEST.mkdir(exist_ok=True)

def main():
    source=PROJECT/'archive/nature/11_submission_ready_v1/source_data/chao_phraya_area_weighted_exposure_cells.csv'
    cells=pd.read_csv(source);x=cells.lon_center.values;y=cells.lat_center.values
    combinations=list(itertools.product(['raw_byte','normalized255'],['legacy_tiepoint','gdal_pixelispoint'],[False,True]))
    counts={c:np.zeros(len(cells),dtype=int) for c in combinations};upath=BACKUP/'revision_scientific_reports_v1/inputs_needed/licsar_chao_28_unw'
    for path in sorted((PROJECT/'raw/licsar').glob('*.geo.cc.tif')):
        with rasterio.open(path) as cd,rasterio.open(upath/path.name.replace('.cc.','.unw.')) as ud:
            ca=cd.read(1);ua=ud.read(1);uv=(ud.read_masks(1)>0)&np.isfinite(ua)&(ua!=0)
            ur,uc=rasterio.transform.rowcol(ud.transform,x,y);ur=np.array(ur);uc=np.array(uc)
            good=(ur>=0)&(ur<ud.height)&(uc>=0)&(uc<ud.width);unw=np.zeros(len(cells),bool);unw[good]=uv[ur[good],uc[good]]
            with tifffile.TiffFile(path) as tf:
                tags=tf.pages[0].geotiff_tags;sx,sy,_=tags['ModelPixelScale'];_,_,_,x0,y0,_=tags['ModelTiepoint']
            for registration in ['legacy_tiepoint','gdal_pixelispoint']:
                if registration=='legacy_tiepoint':
                    cc=np.floor((x-x0)/sx).astype(int);rr=np.floor((y0-y)/sy).astype(int)
                else:rr,cc=map(np.asarray,rasterio.transform.rowcol(cd.transform,x,y))
                inside=(rr>=0)&(rr<cd.height)&(cc>=0)&(cc<cd.width);v=np.zeros(len(cells),dtype='uint8');v[inside]=ca[rr[inside],cc[inside]]
                for encoding in ['raw_byte','normalized255']:
                    passed=inside&(v>=(.3 if encoding=='raw_byte' else 77))
                    for joint in [False,True]:counts[(encoding,registration,joint)]+=passed&(unw if joint else True)
    rows=[]
    for combination,n in counts.items():
        encoding,registration,joint=combination
        for domain,sel in [('geographic',np.ones(len(cells),bool)),('inherited_strong',cells.vlm_mm_yr.le(-5).values)]:
            for exposure,column in [('population','population_weighted'),('builtup_m2','builtup_m2_weighted')]:
                e=cells[column].values*sel
                for statistic,f in [('center_majority_binary',(n>=14).astype(float)),('center_mean_frequency',n/28)]:rows.append(dict(encoding=encoding,registration=registration,joint_unwrapped_validity=joint,domain=domain,exposure=exposure,statistic=statistic,total=float(e.sum()),unsupported=float((e*(1-f)).sum()),low_center_cells=int((sel&(n<14)).sum())))
        cells['_'.join(map(str,combination))]=n
    # The submission-ready table is reproducible under normalized encoding and
    # legacy tiepoint placement; another archive branch used raw-byte values.
    if not np.array_equal(counts[('normalized255','legacy_tiepoint',False)],cells.observable_count.values):raise ValueError('Legacy count reconstruction mismatch')
    pd.DataFrame(rows).to_csv(DEST/'encoding_registration_validity_factorial.csv',index=False)
    cells.to_csv(DEST/'center_counts.csv',index=False)
    archived=[]
    for relative in ['03_exposure_closure/chao_phraya_area_weighted_exposure_censoring','11_submission_ready_v1/source_data','11_submission_ready_v1_stage/source_data']:
        p=PROJECT/'archive/nature'/relative/'chao_phraya_area_weighted_exposure_cells.csv';d=pd.read_csv(p);sel=(d.vlm_mm_yr<=-5)&(d.observable_count<14)
        archived.append(dict(path=str(p),sha256=sha256(p),selected_cells=int(sel.sum()),population=float(d.loc[sel,'population_weighted'].sum()),builtup_km2=float(d.loc[sel,'builtup_m2_weighted'].sum()/1e6)))
    pd.DataFrame(archived).to_csv(DEST/'archive_branches.csv',index=False)
    vlm=PROJECT/'archive/nature/03_exposure_closure/chao_phraya/gridVLM/chaoPhraya_vlm.tif'
    with rasterio.open(vlm) as ds:
        a=ds.read(1);v=a[cells.row,cells.col]
        units=dict(path=str(vlm),sha256=sha256(vlm),band_units=ds.units,scales=ds.scales,offsets=ds.offsets,tags=ds.tags(),band_tags=ds.tags(1),cells=int(len(v)),raw_below_minus5=int((v<=-5).sum()),inherited_times10_below_minus5=int((v*10<=-5).sum()),raw_quantiles=np.quantile(v,[0,.05,.5,.95,1]).tolist(),status='The archived x10 scale is not independently established by the public GeoTIFF unit metadata; retained only as a conditional legacy screen. Main analyses use geographic footprint without a velocity threshold.',source_record='https://zenodo.org/records/18092985',source_processing='Ohenhen et al. WabInSAR, distinct from LiCSAR support sample')
    dump(DEST/'velocity_unit_investigation.json',units)
    iranmeta=DEST/'iran_10815578_metadata.json'
    if not iranmeta.exists():
        response=requests.get('https://zenodo.org/api/records/10815578',timeout=45);response.raise_for_status();dump(iranmeta,response.json())
    official=json.loads(iranmeta.read_text(encoding='utf-8'));lookup={f['key']:f for f in official['files']};iran=[]
    for path in sorted((PROJECT/'archive/nature/05_iran_insar_probe').glob('*.tif')):
        h=hashlib.md5()
        with path.open('rb') as stream:
            for block in iter(lambda:stream.read(2**20),b''):h.update(block)
        actual='md5:'+h.hexdigest();expected=lookup.get(path.name,{}).get('checksum')
        iran.append(dict(filename=path.name,checksum=actual,official_checksum=expected,matches=actual==expected,used_in_new_exposure_experiments=False))
    dump(DEST/'iran_actual_file_match.json',iran)
    dump(DEST/'interpretation.json',dict(original_submission_population=archived[1]['population'],conclusion='Original headline is reproducible from submission-ready tables. Conflicting archive branches exist; do not claim the submitted total is unreproducible.',factorial='Full 2x2x2 comparison holds archived exposure weights and geographic cell centers fixed. Encoding, registration and added unwrapped validity are changed independently.',estimands='A binary majority rule sums entire low-frequency cell populations. Mean frequency deficit averages unsupported population over pairs. Fine area-pair support also resolves area intersections. These are different estimands and must not be presented as an additive code-correction waterfall.',population='Exposure allocation units; neither confirmed affected people nor unique people never observed during the whole period.'))
    print('E0 LINEAGE COMPLETE',flush=True)

if __name__=='__main__':main()
