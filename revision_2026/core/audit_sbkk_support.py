from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
import h5py

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results' / 'gnss_feasibility'
OUT.mkdir(parents=True, exist_ok=True)
# RTSD/National CORS Data Center published SBKK coordinate (JICA report, Table 5)
lat = 13 + 47/60 + 34.59825/3600
lon = 100 + 35/60 + 47.37945/3600
with open(ROOT / 'results' / 'timeseries' / 'geometry.json', encoding='utf-8') as f:
    grid = json.load(f)
a, _, x0, _, e, y0 = grid['transform']
with h5py.File(ROOT / 'results' / 'timeseries' / 'TS_full' / 'cum_filt.h5', 'r') as f:
    dates = pd.to_datetime([str(v) for v in f['imdates'][:]], format='%Y%m%d')
    cum = f['cum'][:]
    mask = f['mask'][:]
    velocity = f['vel'][:]
rows, cols = np.indices(mask.shape)
xs = x0 + (cols + 0.5) * a
ys = y0 + (rows + 0.5) * e
distance_m = np.hypot((xs - lon) * 108000, (ys - lat) * 111000)
within_100m = distance_m <= 100
supported = within_100m & (mask > 0)
nearest_row, nearest_col = np.unravel_index(np.argmin(distance_m), distance_m.shape)
point_offset = float(distance_m[nearest_row, nearest_col])
series = np.nanmedian(np.where(supported[None, :, :], cum, np.nan), axis=(1, 2))
time_years = (dates - dates[0]).days / 365.25
valid_dates = np.isfinite(series)
rate = float(np.polyfit(time_years[valid_dates], series[valid_dates], 1)[0])
output = pd.DataFrame({
    'date': dates.date.astype(str),
    'median_cum_los_mm_within_100m': series,
    'valid_pixel_centers_within_100m': int(supported.sum()),
})
output.to_csv(OUT / 'SBKK_100m_LiCSAR.csv', index=False)
summary = {
    'station': 'SBKK (RTSD CORS)',
    'coordinate_dd': {'latitude': lat, 'longitude': lon},
    'coordinate_source': 'JICA National CORS Data Center report, Table 5; source https://openjicareport.jica.go.jp/pdf/12386280_01.pdf',
    'insar_stack': {
        'n_acquisitions': int(len(dates)),
        'first_date': dates.min().date().isoformat(),
        'last_date': dates.max().date().isoformat(),
        'nearest_pixel_center_offset_m': point_offset,
        'nearest_pixel_center_valid': bool(mask[nearest_row, nearest_col] > 0),
        'pixel_centers_within_100m': int(within_100m.sum()),
        'valid_pixel_centers_within_100m': int(supported.sum()),
        'median_velocity_los_mm_yr': float(np.nanmedian(velocity[supported])),
        'q25_velocity_los_mm_yr': float(np.nanpercentile(velocity[supported], 25)),
        'q75_velocity_los_mm_yr': float(np.nanpercentile(velocity[supported], 75)),
        'median_cumulative_los_rate_mm_yr': rate,
        'cumulative_series_range_mm': float(np.nanmax(series) - np.nanmin(series)),
    },
    'scope': 'Product-specific support audit only. No SBKK GNSS time series was recovered or compared. The LiCSAR LOS rate is not a vertical rate and is not directly comparable with the published SBKK ellipsoidal-height rate.',
}
(OUT / 'SBKK_support_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
print(json.dumps(summary, indent=2))
