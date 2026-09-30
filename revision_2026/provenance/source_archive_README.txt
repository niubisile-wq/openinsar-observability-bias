InSAR scientific revision code, version 0.4.0

Manuscript: Quantifying how spatial allocation and sampling alter population-weighted InSAR support summaries.
Authors: Zixuan Liu and Wei Xiong. Journal submission ID: 20e66410-c472-42dc-ab99-957c24a43f96.

This release consolidates the underlying code for the current scientific revision, including Supplementary S12-S20.
It is a source-code archive with input manifests and reference summaries, not a redistribution of third-party data.
Cite the version-specific DOI on the record hosting this release. The prior DOI 10.5281/zenodo.21444768 identifies older code only.

Directory map
  core/: primary allocation, artificial missingness, sampling, population models, cross-sensitivity,
         matched LiCSBAS processing, figure generation and numerical audits.
  paired_support/: paired input-support/final-mask transfer, direct geometry audit (S13), observed LOS masking (S14).
  followup/: covariance factorial checks (S15), Census support (S16), DWR rate transfer and temporal/persistence
             sensitivity (S17), Bangkok spatial decomposition (S18), Fresno population benchmark (S19), 100-m model sensitivity, and queen-contiguity spatial dependence diagnostic (S20).
  input_manifests/: hashes and logical paths for frozen radar inputs. Historical local prefixes are normalized.
  reference_outputs/: non-spatial aggregate results and method records used to check a reproduction.

Reproduction
Use Python 3.10. The core and LiCSBAS environments are recorded separately under core/. Do not merge the
two pinned NumPy versions. Install the environment appropriate to the analysis; the follow-up analysis uses
the recorded LiCSBAS environment plus pyproj. Full package versions are in environment_licsbas.json.

python reproduce.py --help
python reproduce.py check-code
python reproduce.py followup --strengthening-root PATH_TO_STRENGTHENING --workspace PATH_TO_FOLLOWUP_INPUTS
python reproduce.py population-benchmark --workspace PATH_TO_FOLLOWUP_INPUTS

Input layout
The strengthening root contains results/timeseries/matched_quality.npz and geometry.json, results/dwr/grid.npz,
external/dwr/roi_features.json and external/worldpop_thailand_2020_UNadj.tif. For full raw reprocessing, fetch
the radar catalogue and population inputs using core/fetch_* scripts and the frozen manifests, and install
LiCSBAS2 separately. The full matched processing inputs are not bundled.

For S20, run python followup/fresno_capop_support_sensitivity.py then python followup/spatial_dependence_audit.py after acquiring the S20 inputs; the latter reports descriptive queen-contiguity global Moran I without spatial sampling inference.

The follow-up workspace contains census_blocks/fresno_roi/blocks_*.geojson (6,999 unique queried GEOIDs),
inputs/bangkok_matched/matched_quality.npz, inputs/fresno_population/ (the exact GHSL/WorldPop grids),
dwr_velocity/ (documented annual and endpoint rasters), and external/round3/ (Census P1, TIGER block subset, CA-POP 100-m grid for S20). See INPUT_ACQUISITION.txt and method JSONs.
External inputs are required: this archive is not a self-contained raw-data reproduction bundle.

Path changes
Source operations, thresholds, random seeds and calculations are retained. Only follow-up and paired-support
path configuration is adapted to environment variables. source_provenance.json records both source hashes.
The dispatcher writes analysis outputs into the selected workspace. It does not deposit code or submit a manuscript.
Historical prose generators under core/ are retained as source history; do not run them over the current manuscript.

Validation
check-code verifies Python syntax and archive hashes without input data. Data-based rerun checks are documented
in the submission handoff verification. These checks establish computational consistency for the defined inputs;
they do not establish independent geodetic validation or true exposed-person counts.

Round 3 response-strengthening analyses (added 2026-09-30)
Set PATH_TO_STRENGTHENING to the retained workspace containing results/timeseries/matched_quality.npz.
python reproduce.py partial-identification --strengthening-root PATH_TO_STRENGTHENING --workspace OUTPUT
python reproduce.py published-method-benchmark --strengthening-root PATH_TO_STRENGTHENING --workspace OUTPUT
python reproduce.py reporting-tolerance --workspace OUTPUT
The first reports finite-domain LOS-class allocation bounds under the actual final mask. The second compares census-unit median-rate, within-unit observed-population prevalence and centre sampling under fixed artificial masks on originally valid LOS pixels. The third reads the retained Fresno scale_method_sensitivity.csv and evaluates user-declared illustrative tolerances. These are sensitivity analyses; none validates deformation or population in genuinely unsupported locations. See the local manuscript Supplementary S21-S23 and method.json files.

ACS survey benchmark: run reproduce.py acs-population-benchmark --workspace <analysis-workspace>. The workspace must contain inputs/acs_2021, inputs/fresno_population, and per_block.csv. This compares 2020 population grids with ACS 2017-2021 block-group estimates; it is a coarse survey benchmark, not household-location validation.
