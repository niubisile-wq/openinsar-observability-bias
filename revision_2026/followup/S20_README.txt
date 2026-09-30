Supplementary S20: Fresno fine-resolution population-allocation sensitivity

Code: followup/fresno_capop_support_sensitivity.py
Results: reference_outputs/followup/fresno_population_allocation/
Input URL and SHA-256 records: followup/fresno_population_source_manifest.json

The script expects the strengthening project layout at the package root:
  external/round3/ca2020.pl.zip
  external/round3/census_blocks/tl_2020_06019_tabblock20.shp (+ sidecars)
  external/round3/CAPOP_2020_100m_TOTAL.tif
  results/dwr/grid.npz
The U.S. Census PL94-171 P1 file and 2020 TIGER/Line tabulation-block polygons are official sources.
The CA-POP 100-m grid is Depsky et al. (2022), DOI 10.5281/zenodo.5874927, derived by dasymetric
allocation from the same 2020 Census block counts. Inputs are not redistributed in this code archive.
Use the manifest URLs and hashes, then place the county subset of the California TIGER/Line 2020
Tabulation Block file in the listed path. Run with the recorded follow-up Python environment.

Interpretation: the comparison conditions both 100-m grids on the same block totals. It is allocation-model
sensitivity, not independent population validation. Three blocks (45 people) lack GHSL mass and are excluded
from the common comparison domain. Geometry-source block counts differ by one from S16-S19, while complete-
block population totals agree at 523,953. No inferential p-values are reported for the fixed regional domain.
