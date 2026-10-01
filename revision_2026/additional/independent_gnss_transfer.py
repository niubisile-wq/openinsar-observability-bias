"""Audit independent GNSS observations for the two fixed transfer regions.

GNSS LOS estimates are extracted from the provider's completed notebooks.
This reanalysis is independent of our population/support calculations, not a
fresh fit of every daily GNSS series or a validation of Bangkok missing rates.
"""
from pathlib import Path
import argparse
import json
import hashlib
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window
from pyproj import Geod

REGIONS=[('Los_Angeles','A064',[-118.4,33.85,-118.1,34.1],'P597','MHMS','20160214','20220513'),
         ('Davis_Sacramento','A137',[-121.9,38.4,-121.6,38.7],'P266','UCD1','20170303','20230829')]

def sample_stations(inputs,track,reference,selected,radius,geod,sample_source):
    """Read either complete source raster or exact archived 7x7 source windows."""
    estimates={};station_info={};excluded=[]
    path=inputs/'aria_california'/f'{track}_velocity_ref{reference}.tif'
    full=(sample_source=='full' or (sample_source=='auto' and path.is_file()))
    if full:
        with rasterio.open(path) as ds:
            for _,station in selected.iterrows():
                row,col=ds.index(station.longitude,station.latitude)
                if not (radius<=row<ds.height-radius and radius<=col<ds.width-radius):
                    excluded.append({'station':station.station,'pixel_radius':radius,'reason':'outside raster'});continue
                vals=ds.read(1,window=Window(col-radius,row-radius,2*radius+1,2*radius+1)).astype(float)
                xx,yy=ds.xy(row,col)
                station_info[station.station]={'sample_row':row,'sample_col':col,'pixel_center_longitude':xx,
                                              'pixel_center_latitude':yy,'values':vals}
    else:
        with np.load(inputs/'gnss_reference'/f'{track}_station_windows.npz',allow_pickle=False) as capsule:
            windows={str(s):w for s,w in zip(capsule['stations'],capsule['source_windows'])}
        meta=pd.read_csv(inputs/'gnss_reference'/f'{track}_station_window_metadata.csv').set_index('station')
        for _,station in selected.iterrows():
            name=station.station
            assert name in windows and name in meta.index
            vals=windows[name][3-radius:4+radius,3-radius:4+radius].astype(float)
            station_info[name]={**meta.loc[name].to_dict(),'values':vals}
    for _,station in selected.iterrows():
        name=station.station
        if name not in station_info:continue
        info=station_info[name];vals=info.pop('values')
        finite=np.isfinite(vals)
        if name!=reference:finite&=vals!=0
        if not finite.any():
            excluded.append({'station':name,'pixel_radius':radius,'reason':'no published finite nonzero LOS pixel'});continue
        estimates[name]=float(np.median(vals[finite])*1000)
        offset=float(geod.inv(station.longitude,station.latitude,info['pixel_center_longitude'],info['pixel_center_latitude'])[2])
        station_info[name]={'sample_row':int(info['sample_row']),'sample_col':int(info['sample_col']),
                            'station_to_pixel_center_m':offset,'finite_sample_pixels':int(finite.sum())}
    return estimates,station_info,excluded


def analyse(inputs,out,sample_source='auto'):
    out.mkdir(parents=True,exist_ok=True);rows=[];summaries=[];exclusions=[];clusters=[]
    geod=Geod(ellps='WGS84')
    for region,track,bbox,reference,local_cal,start,end in REGIONS:
        stations=pd.read_csv(inputs/'gnss_reference'/f'{track}_archived_rates_with_coordinates.csv')
        roi=stations.longitude.between(bbox[0],bbox[2])&stations.latitude.between(bbox[1],bbox[3])
        selected=stations[roi|stations.station.eq(reference)].copy()
        assert selected.station.eq(reference).sum()==1 and selected.station.eq(local_cal).sum()==1
        # Local frame-control site is nearest fixed box center among spatially/time-eligible archived sites.
        center=[(bbox[0]+bbox[2])/2,(bbox[1]+bbox[3])/2]
        regional=stations[roi].copy()
        distances=geod.inv(np.full(len(regional),center[0]),np.full(len(regional),center[1]),regional.longitude.to_numpy(),regional.latitude.to_numpy())[2]
        assert regional.iloc[np.argmin(distances)].station==local_cal
        groups={};next_group=0
        for idx,station in regional.sort_values('station').iterrows():
            if station.station in groups:continue
            distances=geod.inv(np.full(len(regional),station.longitude),np.full(len(regional),station.latitude),regional.longitude.to_numpy(),regional.latitude.to_numpy())[2]
            near=regional.loc[distances<=100,'station'].tolist()
            for name in near:groups[name]=next_group
            next_group+=1
        for station,group in groups.items():clusters.append({'region':region,'station':station,'collocation_cluster_100m':group})
        for radius in [0,1,3]:
                estimates,station_info,missing=sample_stations(inputs,track,reference,selected,radius,geod,sample_source)
                exclusions.extend({'region':region,**x} for x in missing)
                if reference not in estimates or local_cal not in estimates:raise RuntimeError('Predefined calibration/reference site is not evaluable')
                g=selected.set_index('station').provider_gnss_los_velocity_mm_yr.to_dict()
                for frame_control in [reference,local_cal]:
                    block=[]
                    for _,station in selected.iterrows():
                        name=station.station
                        # No reference/frame-control site counted in validation.
                        if name in [reference,local_cal] or name not in estimates:continue
                        los=estimates[name]-estimates[frame_control];gnss=g[name]-g[frame_control]
                        value={'region':region,'track':track,'station':name,'frame_control':frame_control,'pixel_radius':radius,
                               'start_date':start,'end_date':end,'collocation_cluster_100m':groups[name],
                               'published_insar_relative_los_mm_yr':los,'archived_gnss_relative_los_mm_yr':gnss,
                               'signed_error_mm_yr':los-gnss,**station_info[name]}
                        rows.append(value);block.append(value)
                    frame=pd.DataFrame(block)
                    if len(frame)<3:raise RuntimeError('Fewer than three predefined validation stations')
                    cluster=frame.groupby('collocation_cluster_100m').signed_error_mm_yr.median()
                    error=frame.signed_error_mm_yr.to_numpy()
                    summaries.append({'region':region,'track':track,'frame_control':frame_control,'pixel_radius':radius,
                                      'spatially_eligible_archived_stations':len(regional),
                                      'validation_stations':len(frame),'validation_collocation_clusters':len(cluster),
                                      'mean_signed_error_mm_yr':float(error.mean()),'MAE_mm_yr':float(np.abs(error).mean()),
                                      'RMSE_mm_yr':float(np.sqrt(np.mean(error**2))),
                                      'collocation_cluster_RMSE_mm_yr':float(np.sqrt(np.mean(cluster.to_numpy()**2))),
                                      'min_signed_error_mm_yr':float(error.min()),'max_signed_error_mm_yr':float(error.max()),
                                      'max_station_to_pixel_center_m':float(frame.station_to_pixel_center_m.max())})
        print('Independent GNSS transfer',region,'spatial candidates',len(regional),flush=True)
    pd.DataFrame(rows).to_csv(out/'gnss_transfer_by_station.csv',index=False)
    pd.DataFrame(summaries).to_csv(out/'gnss_transfer_summary.csv',index=False)
    pd.DataFrame(exclusions,columns=['region','station','pixel_radius','reason']).to_csv(out/'gnss_transfer_exclusions.csv',index=False)
    pd.DataFrame(clusters).to_csv(out/'gnss_transfer_collocation_clusters.csv',index=False)
    protocol={'reference_type':'Independent NGL GNSS LOS estimates archived by Sangha et al. 2026, corresponding to each track period and source geometry.',
              'reference_processing':'Provider notebook cell 38 omits the model argument. Archived MintPy time_func.py confirms default polynomial degree 1, without seasonal/step terms; InSAR trend includes annual/semiannual and earthquake-step terms. This temporal-model mismatch is a limitation.',
              'calibration_validation_disjoint':True,'local_control_selection':'Nearest box-center archived station, using coordinates alone: MHMS / UCD1; these and provider reference P597 / P266 excluded from every reported validation metric.',
              'frame_method':'Subtract the same named station from both LOS products; no regression/ramp fit using validation stations.',
              'sample_method':'Containing native raster cell using floor through rasterio.index; 1/3-pixel-radius median sensitivity; no rounded pixel-corner indexing.',
              'dependence':'Station errors, shared reference and nearby monuments are dependent; metrics are descriptive, with no iid confidence interval.',
              'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((inputs/'gnss_reference').glob('*archived_rates_with_coordinates.csv'))},
              'checks':{'no_validation_site_used_for_frame_control':True,'local_control_chosen_from_location_only':True,
                        'all_eligible_archived_stations_and_unavailable_samples_retained':True,'collocated_sites_not_counted_as_independent_replication':True},
              'limits':['Reanalysis of a published external LOS product, not independent validation of our Bangkok inversion or real missing cells.',
                        'Archived GNSS trend outputs preserve the provider station screening and processing; raw daily-series refitting has not been claimed.',
                        'Exact GNSS-projection geometry and time-function choices belong to the provider; archived notebook and implementation preserved for audit.',
                        'Monument/radar scattering targets, source screening, vintage and processing choices may differ.',
                        'Small residuals would not validate population placement or support masks, and large residuals are retained.']}
    (out/'gnss_transfer_protocol.json').write_bytes((json.dumps(protocol,indent=2)+'\n').encode())
    print(pd.DataFrame(summaries).to_string(index=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path);p.add_argument('--output',type=Path);p.add_argument('--sample-source',choices=['auto','full','windows'],default='auto');args=p.parse_args();root=Path(__file__).resolve().parent.parent
    analyse(args.inputs or root/'inputs',args.output or root/'results/gnss_transfer',args.sample_source)
