from pathlib import Path
import hashlib, json, math
import numpy as np
import pandas as pd
import h5py
import rasterio
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results' / 'gnss_feasibility'
OUT.mkdir(parents=True, exist_ok=True)
STATIONS = {
    'CUSV': {'lon': 100.5339209025, 'lat': 13.7359146691},
    'CUUT': {'lon': 100.533935, 'lat': 13.735991},
}
heading_deg = -168.882966516
window_days = 7

with open(ROOT / 'results' / 'timeseries' / 'geometry.json', encoding='utf-8') as f:
    grid = json.load(f)
a, _, x0, _, e, y0 = grid['transform']
with h5py.File(ROOT / 'results' / 'timeseries' / 'TS_full' / 'cum_filt.h5', 'r') as f:
    dates = pd.to_datetime([str(v) for v in f['imdates'][:]], format='%Y%m%d')
    cum = f['cum'][:]
    mask = f['mask'][:]
    velocity = f['vel'][:]

geometry_path = ROOT / 'external' / 'licsar_metadata' / '062D_07629_131313.geo.U.tif'
source_meta = {}
with rasterio.open(geometry_path) as src:
    for name, loc in STATIONS.items():
        loc['u'] = float(next(src.sample([(loc['lon'], loc['lat'])]))[0])
        loc['incidence_deg'] = float(np.rad2deg(np.arccos(loc['u'])))
        az = np.deg2rad(heading_deg) + np.pi / 2
        inc = np.arccos(loc['u'])
        loc['los_enu'] = [float(np.sin(inc) * np.sin(az)),
                          float(np.sin(inc) * np.cos(az)), loc['u']]

for station, loc in STATIONS.items():
    gnss_path = ROOT / 'external' / 'gnss' / f'{station}.tenv3'
    d = pd.read_csv(gnss_path, sep=r'\s+', comment='#')
    d['date'] = pd.to_datetime(d['YYMMMDD'], format='%y%b%d')
    for component, integer, fraction in [
        ('E', '_e0(m)', '__east(m)'),
        ('N', '____n0(m)', '_north(m)'),
        ('U', 'u0(m)', '____up(m)'),
    ]:
        d[component] = (d[integer] + d[fraction]) * 1000.0
    d['los_mm'] = d[list('ENU')].to_numpy() @ np.asarray(loc['los_enu'])
    d['los_formal_sigma_mm'] = np.sqrt(
        (d['sig_e(m)'].to_numpy() * 1000 * loc['los_enu'][0])**2
        + (d['sig_n(m)'].to_numpy() * 1000 * loc['los_enu'][1])**2
        + (d['sig_u(m)'].to_numpy() * 1000 * loc['los_enu'][2])**2
    )

    col = int(round((loc['lon'] - x0) / a))
    row = int(round((loc['lat'] - y0) / e))
    loc['row'], loc['col'] = row, col
    loc['grid_center_lon'] = float(x0 + (col + 0.5) * a)
    loc['grid_center_lat'] = float(y0 + (row + 0.5) * e)
    loc['pixel_offset_m'] = float(np.hypot(
        (loc['grid_center_lon'] - loc['lon']) * 108000,
        (loc['grid_center_lat'] - loc['lat']) * 111000))
    loc['valid_pixel'] = bool(mask[row, col] > 0)
    loc['pixel_velocity_mm_yr'] = float(velocity[row, col])
    cube = cum[:, row-1:row+2, col-1:col+2]
    valid = mask[row-1:row+2, col-1:col+2] > 0
    loc['valid_pixels_3x3'] = int(valid.sum())
    point = cum[:, row, col].astype(float)
    neighborhood = np.nanmedian(np.where(valid[None, :, :], cube, np.nan), axis=(1, 2))

    sampled = []
    for date in dates:
        w = d[(d.date >= date - pd.Timedelta(days=window_days))
              & (d.date <= date + pd.Timedelta(days=window_days))]
        if w.empty:
            sampled.append((np.nan, np.nan, 0))
        else:
            sampled.append((float(w.los_mm.median()),
                            float(w.los_formal_sigma_mm.median()), len(w)))
    g = np.array([v[0] for v in sampled], dtype=float)
    gs = np.array([v[1] for v in sampled], dtype=float)
    ng = np.array([v[2] for v in sampled], dtype=int)

    rows_out = []
    for i, date in enumerate(dates):
        rows_out.append({
            'station': station, 'date': date.date().isoformat(),
            'gnss_projected_los_mm': g[i], 'gnss_formal_sigma_mm': gs[i],
            'gnss_days_in_15d_window': int(ng[i]),
            'insar_point_cum_mm': float(point[i]),
            'insar_3x3_median_cum_mm': float(neighborhood[i]),
            'mask_valid_point': int(loc['valid_pixel']),
        })
    pd.DataFrame(rows_out).to_csv(OUT / f'{station}_acquisition_match.csv', index=False)

    good = np.isfinite(g) & np.isfinite(point)
    inds = np.where(good)[0]
    result = {
        'station': station,
        'coordinate': {'longitude': loc['lon'], 'latitude': loc['lat']},
        'gnss_daily_series': {
            'n_rows': int(len(d)), 'start': d.date.min().date().isoformat(),
            'end': d.date.max().date().isoformat(),
            'matched_acquisitions': int(good.sum()), 'total_acquisitions': int(len(dates)),
            'matching_rule': f'median of GNSS solutions within +/-{window_days} days',
        },
        'insar_match': {
            'nearest_grid_center': [loc['grid_center_lon'], loc['grid_center_lat']],
            'point_offset_m': loc['pixel_offset_m'],
            'valid_point': loc['valid_pixel'],
            'valid_pixels_3x3': loc['valid_pixels_3x3'],
            'point_velocity_mm_yr': loc['pixel_velocity_mm_yr'],
        },
        'los_unit_vector_ENU': loc['los_enu'],
        'local_incidence_deg_from_geo_U': loc['incidence_deg'],
        'metrics': {},
        'interpretation': 'Feasibility screen only; not a physical validation. InSAR is relative to an internal reference area while GNSS is in a global reference frame; a co-located GNSS reference or documented common-frame correction is unavailable. The station is also rooftop-mounted and has a nonzero pixel offset.'
    }
    for label, ins in [('point', point), ('3x3_median', neighborhood)]:
        ok = np.isfinite(g) & np.isfinite(ins)
        ids = np.where(ok)[0]
        base = ids[0]
        x = ins[ok] - ins[base]
        y = g[ok] - g[base]
        t = np.array([(dates[i] - dates[base]).days / 365.25 for i in ids])
        xres = x - np.polyval(np.polyfit(t, x, 1), t)
        yres = y - np.polyval(np.polyfit(t, y, 1), t)
        result['metrics'][label] = {
            'n': int(ok.sum()), 'first_date': dates[base].date().isoformat(),
            'last_date': dates[ids[-1]].date().isoformat(),
            'insar_linear_rate_mm_yr': float(np.polyfit(t, x, 1)[0]),
            'gnss_projected_linear_rate_mm_yr': float(np.polyfit(t, y, 1)[0]),
            'centered_pearson_r': float(pearsonr(x, y).statistic),
            'centered_rmse_mm': float(np.sqrt(np.mean((x-y)**2))),
            'detrended_pearson_r': float(pearsonr(xres, yres).statistic),
            'detrended_rmse_mm': float(np.sqrt(np.mean((xres-yres)**2))),
            'median_projected_gnss_formal_sigma_mm': float(np.nanmedian(gs[ok])),
        }
    (OUT / f'{station}_summary.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    source_meta[station] = {
        'series_url': f'https://geodesy.unr.edu/gps_timeseries/IGS20/tenv3/IGS20/{station}.tenv3',
        'station_page': f'https://geodesy.unr.edu/NGLStationPages/stations/{station}.sta',
        'series_sha256': hashlib.sha256(gnss_path.read_bytes()).hexdigest(),
        'matched_acquisitions': int(good.sum()),
        'metrics_point': result['metrics']['point'],
        'metrics_3x3': result['metrics']['3x3_median'],
        'local_incidence_deg': loc['incidence_deg'],
        'los_enu': loc['los_enu'],
        'point_offset_m': loc['pixel_offset_m'],
        'valid_pixels_3x3': loc['valid_pixels_3x3'],
    }

(OUT / 'source_manifest.json').write_text(json.dumps({
    'generated_utc_date': '2026-09-27',
    'n_gloss': 'NGL IGS20 tenv3 station time series',
    'stations': source_meta,
    'licsar_geometry': {
        'file': 'external/licsar_metadata/062D_07629_131313.geo.U.tif',
        'sha256': hashlib.sha256(geometry_path.read_bytes()).hexdigest(),
        'description': 'Up component of LiCSAR LOS unit vector; East/North components derived from metadata heading for feasibility screen.',
    },
    'timeseries': {
        'file': 'results/timeseries/TS_full/cum_filt.h5',
        'dates': [dates.min().date().isoformat(), dates.max().date().isoformat()],
        'n_acquisitions': int(len(dates)),
    },
}, indent=2), encoding='utf-8')
print(json.dumps(source_meta, indent=2))

