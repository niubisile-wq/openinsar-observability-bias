"""One declared year-length convention for DWR interval displacement products."""
from datetime import date

DAYS_PER_YEAR=365.2425
WINDOWS={
 'rate_2015_2016.tif':('2015-10-01','2016-10-01'),
 'annual_oct_2016_2017.tif':('2016-10-01','2017-10-01'),
 'annual_oct_2017_2018.tif':('2017-10-01','2018-10-01'),
 'annual_oct_2018_2019.tif':('2018-10-01','2019-10-01'),
 'annual_oct_2019_2020.tif':('2019-10-01','2020-10-01'),
 'rate_2020_2021.tif':('2020-10-01','2021-10-01'),
 'rate_2023_2024.tif':('2023-10-01','2024-10-01'),
 'rate_2025_2026.tif':('2025-04-01','2026-04-01'),
 'total_since_20150613_20210101.tif':('2015-06-13','2021-01-01')}

def duration_years(filename):
    start,end=WINDOWS[filename]
    return (date.fromisoformat(end)-date.fromisoformat(start)).days/DAYS_PER_YEAR

def feet_to_feet_per_year(values,filename):
    # Derived composite has already been annualized separately in every window.
    if filename=='rate_median_2015_2021.tif':return values
    return values/duration_years(filename)
