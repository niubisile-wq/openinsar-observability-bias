# Reproduce the 2026 scientific revision

**Manuscript:** *Quantifying how spatial allocation and sampling alter population-weighted InSAR support summaries* (Scientific Reports revision).

**Version scope:** v0.5.1 integrates six additional experiments into the manuscript, Supplementary S25-S30, all reviewer responses, scientific figures and regenerated line/page locators. The fixed scientific release is [v0.5.0](https://github.com/niubisile-wq/openinsar-observability-bias/releases/tag/v0.5.0), commit `f70124f1cdc11e43848ce97c01f49fad6f7efb62`, including all portable inputs, 21 additional numerical comparisons and two selected-source assets. See the [additional experiment guide](additional/README.md) and run `python revision_2026/reproduce.py additional`. The preceding S1-S24 scientific baseline remains v0.4.1, commit `5bf98f19abbed981989f9c97ef9b118eea524c5b`. The [v0.5.1 document asset](https://github.com/niubisile-wq/openinsar-observability-bias/releases/download/v0.5.1/reviewer-documents-v0.5.1.zip) contains the six final PDFs and two independently recompiled, self-contained LaTeX archives. Document verification records identify the exact current snapshot.

This directory is the current revision. The older `14_esin_strengthened_v1/` manuscript and Zenodo DOI `10.5281/zenodo.21444768` describe an earlier version. That DOI does **not** archive the new experiments. The updated Zenodo deposit remains pending.

The checkout contains current manuscript/Supplementary Information PDFs and LaTeX, all retained experiment code, frozen analysis inputs, full reference tables including replicate-level results, input provenance, and numerical comparison commands. The supplied arrays are processed analysis inputs. Recomputing from these arrays reproduces the reported statistical analyses; it does not independently reproduce every raw-interferogram processing step or establish external geodetic accuracy.

## Start here

Use **Python 3.10**, a 64-bit environment, and the pinned dependency file. From the repository root:

Linux/macOS (with `python3.10` available):

```sh
python3.10 -m venv .venv
.venv/bin/python -m pip install -r revision_2026/requirements_reviewer_lock.txt
.venv/bin/python revision_2026/reproduce.py check-code
.venv/bin/python revision_2026/reproduce.py quick --workspace reproduction_output
```

Windows PowerShell (no activation-policy change is needed):

```powershell
py -3.10 -m venv .venv
& .\.venv\Scripts\python.exe -m pip install -r revision_2026/requirements_reviewer_lock.txt
& .\.venv\Scripts\python.exe revision_2026/reproduce.py check-code
& .\.venv\Scripts\python.exe revision_2026/reproduce.py quick --workspace reproduction_output
```

For the commands below, `python` means the executable of that installed environment.

The `quick` command performs actual data-based reruns of S21–S24 and compares six scientific tables with the frozen reference values. It uses the bundled processed LOS arrays, official Census survey-value/geography subsets, and the attributed Fresno population rasters. It requires no data-provider account or original Windows drive paths.

The end of a successful run prints `"status": "PASS"`. The machine-readable comparison, timing, row counts and maximum numerical differences are saved to `reproduction_output/reproduction_verification.json`. `check-code` alone checks integrity and Python syntax; it does not run an experiment.

## Run the wider experiment suite

```sh
python revision_2026/reproduce.py standard --workspace reproduction_output
```

This includes the four quick analyses, E1 allocation, E3 synthetic masking in both regions, E4 conditional sampling/nested summaries, E5 population/land-cover comparisons, E6 Fresno allocation, crossed population/sampling sensitivity, paired support, direct polygon intersections, actual-LOS masking, and all six additional S25-S30 analyses. The two synthetic experiments each generate 21,600 replicate rows, actual-LOS masking generates 57,600 rows, and the published-method benchmark generates 4,800 rows. S20 provider-download modes and the Census/DWR polygon mode remain separate commands below.

Run individual modes with the same `--workspace` argument:

| Analysis | Mode | Reproduction scope |
| --- | --- | --- |
| E1, 512 matched scale/origin/method cases | `allocation` | Recompute from native allocation/support arrays |
| E3, Bangkok artificial missingness | `synthetic-masking` | Regenerate masks and synthetic reference fields with fixed seed |
| E3, Fresno artificial missingness | `dwr-masking` | Same design on frozen DWR footprint/population arrays |
| E4 | `sampling-analysis`, `sampling-spatial` | Recompute nested/temporal/conditional summaries and spatial rankings from frozen per-pair integrals and nested mean fields |
| E5 | `population-analysis` | Recompute totals, coverage, epoch and land-cover summaries from aligned model arrays |
| E6 | `dwr-allocation` | Recompute allocation using the frozen union-of-footprints field |
| Population × pair-sample sensitivity | `cross-sensitivity` | Fixed common domain, conserved cross-moments and alternate integration resolutions |
| S13 | `paired`, `direct-geometry` | Reproject and directly intersect the frozen matched-grid arrays |
| S14 | `observed-los-masking` | Regenerate controlled masks on originally valid LOS pixels |
| S15 | `followup` | Recompute the covariance/factorial and median-rule diagnostics; the bundled WorldPop matched-grid array is frozen before normalization/outcome calculation |
| S18 | `spatial` | Recompute the Bangkok spatial decomposition on the fixed matched processing grid |
| S19 | `population-benchmark` | Redo Fresno Census-block geometry/raster intersections; the frozen 6,999-block query includes the 6,788 complete-block domain |
| S20 | `fine-population`, `spatial-dependence` | Download upstream Census/CA-POP files first |
| S21 | `partial-identification` | Recompute deterministic LOS-class allocation bounds |
| S22 | `published-method-benchmark` | Recompute the three estimator rules under fixed controlled masks |
| S23 | `reporting-tolerance` | Recompute declared tolerance decisions from the frozen S20 scale/origin table |
| S24 | `acs-population-benchmark` | Redo geometry overlay and exact 80-replicate aggregate ACS MOE |
| S25 | `additional-real-blocks` | Actual Census-block geometry and random/one-hole/four-hole masks; retain undefined results and partial-identification bounds |
| S26 | `additional-dependence` | Buffered training-only recalibration and deletion of all edges sharing each acquisition date |
| S27 | `additional-quality` | All 18 training-only post-filter rules on 44 excluded-edge residuals |
| S28 | `additional-regions` | Both fixed external regions and all 96 threshold/scale/origin designs |
| S29 | `additional-residential` | Two municipal residential definitions, six conserved allocations and incomplete common-domain coverage |
| S30 | `additional-gnss` | All validation sites, two excluded-control reference frames and three source-window sampling radii |

The standard suite creates sizable intermediate reprojection arrays. Reserve several GB of free memory and disk space. Outputs are written to the chosen workspace, leaving the released input/reference files unchanged. `compare` checks all available rerun tables:

```sh
python revision_2026/reproduce.py compare --workspace reproduction_output
```

Raw totals/areas are compared with relative tolerance `1e-8` and absolute tolerance `1e-6`; fraction/share/percentage columns use absolute tolerance `1e-9`. These cover small numerical-library/reprojection rounding differences. Support-distribution threshold comparisons explicitly use float64 so NumPy 1 and NumPy 2 scalar-promotion rules select the same pixels. Rows/columns, labels and missing-value positions must also match. Local measured differences and timings are retained in the current reproduction report.

## Download and reproduce S20 from provider files

CA-POP and the state/county Census source archives are acquired from their providers rather than bundled in Git:

```sh
python revision_2026/reproduce.py fetch-inputs --dataset s20 --workspace reproduction_output
python revision_2026/reproduce.py fine-population --workspace reproduction_output
python revision_2026/reproduce.py spatial-dependence --workspace reproduction_output
```

The downloader verifies all retained SHA-256 hashes, extracts the official TIGER and CA-POP archives, and fails on a changed artifact. The CA-POP URL is resolved from its actual Zenodo record, including the owner-prefixed archive name. S20 uses **6,787** complete official TIGER blocks; its common model domain is **5,713** blocks and **523,908** population units. The earlier S16–S19 Esri-hosted geometry has **6,788** complete blocks; both complete-block geometries sum to 523,953 people. Do not interchange the geometries silently.

## Full ACS source-file route

For convenience, the [v0.4.1 release](https://github.com/niubisile-wq/openinsar-observability-bias/releases/tag/v0.4.1) includes a checksum-verified companion archive with all six complete ACS source files plus the two official Census archives used in S20. The archive contains a short placement guide; the downloader below remains available if you prefer a direct provider download. CA-POP is intentionally obtained from Zenodo by the downloader and is not redistributed in the companion archive.

The compact ACS inputs contain all 637 Fresno County block groups, their official TIGER geometry, sequence-file estimates/MOEs and 80 variance replicates. The analysis still applies the fixed ROI criterion and selects 319 complete groups itself; neither raster-model predictions nor selected-outcome totals are embedded in that subset. Parent-file URLs and hashes are recorded under `data/workspace/inputs/acs_compact/PROVENANCE.json`.

To rerun from the complete provider archives instead:

```sh
python revision_2026/reproduce.py fetch-inputs --dataset acs --workspace reproduction_output
python revision_2026/reproduce.py acs-population-benchmark --workspace reproduction_output
```

When the full source files are present, the original sequence/geography join and variance-replicate parsing are used. The output comparison remains against the same frozen 319-group tables. Aggregate MOE uses the 80 replicate **sums**, not an independence-based sum of individual MOEs.

## Raw processing and remaining provider inputs

For the S16 Census/DWR footprint calculation, run `fetch-inputs --dataset dwr-polygons` and then `census-support`, with the same workspace. The downloader checks the 38,738 unique observation-cell codes before calculation. The small retained Census block query is already bundled; `fetch-inputs --dataset census-blocks` retrieves the same predeclared query anew and verifies 6,999 unique GEOIDs.

`INPUT_ACQUISITION.txt`, `input_manifests/`, the retained `core/fetch_*` scripts and the two original environment records describe the raw-data route. LiCSBAS2 is pinned to commit `2627e50a1c703becd4a0b18754b05c1cf8c564d5`; its processing environment uses NumPy 1.26.4, while the original main-experiment environment used NumPy 2.2.6. Keep those raw-processing environments separate. The unified reviewer environment is for replaying the released analysis arrays and is checked against the reported tables with stated floating-point tolerances.

The 437-pair full/training LiCSBAS stacks and all 224 pair-level fine support arrays are not in the small Git checkout. E4's fast route starts after pair integration; the original `core/experiment_sampling.py` without `--frozen-summaries` retains pixel-transfer code for a fully acquired input stack. DWR/TRE source velocity rasters and original observation polygons are described in the provider manifests; interpolated vertical-rate products are distinct from the Bangkok LOS product. Older prose-generation scripts are preserved for source history and should not overwrite the current manuscript.

## Locate the paper and evidence

The preceding S1-S24 publication audit covered 21 distinct analysis modes, 49 scientific-table comparisons across 22 data-based runs, and five geometry tests; full ACS and Census/DWR download routes were also executed. Its record remains in `provenance/current_reproduction_verification.json`. The new six-mode numerical rerun compares all 21 additional tables; GitHub Actions executes both S21-S24 and S25-S30 on Windows and Linux. Current audit records are retained separately, so the historical audit is not described as a run of the new manuscript. The current standard rerun completed 24 modes and 62 numerical table comparisons; see `provenance/current_standard_reproduction_v0.5.0.json`. These counts already include the six additional modes and their 21 tables; they are not added again.

- [`manuscript/manuscript.pdf`](manuscript/manuscript.pdf): current clean manuscript.
- [`manuscript/supplementary.pdf`](manuscript/supplementary.pdf): Supplementary S1–S30.
- [`manuscript/response.pdf`](manuscript/response.pdf): point-by-point response to all seven reviewer comments.
- [`REVIEWER_MAP.md`](REVIEWER_MAP.md): comments, changes, experiment entries and scope limits.
- `reference_outputs/`: retained scientific tables and full replicate tables.
- `data/`: frozen analysis inputs and public-data subsets with provenance.
- `MANIFEST.sha256.json`: release integrity for code, inputs, outputs and manuscripts.
- `provenance/`: environment/build records and the current reproduction verification.

The LaTeX source uses the included Springer Nature class/styles and prebuilt `.bbl` files. Compile from `manuscript/` with `pdflatex manuscript.tex` twice and likewise `supplementary.tex` and `response.tex`, or use `latexmk -pdf`. The shipped PDFs are the verified manuscript snapshot; regenerated layout depends on the TeX distribution.

## Interpretation and attribution

These analyses concern population-weighted support, allocation, aggregation and sampling. GHSL/WorldPop/CA-POP are allocation models. ACS is a coarse survey benchmark with uncertainty; CA-POP shares the Census block totals. Controlled masks hide originally observed data and do not validate genuinely unsupported locations. S30 provides an archived independent-GNSS comparison of external California LOS products, with conflicting source-date metadata, different temporal models and large negative residuals explicitly retained. It does not validate our Bangkok inversion; no matched independent GNSS/leveling accuracy validation of that inversion has been obtained. The earlier unvalidated exposed-person interpretation is withdrawn.

Author-owned code and aggregate outputs retain the release's CC BY 4.0 terms. Third-party data, derived arrays and LaTeX support files retain their applicable provider terms; see [`DATA_LICENSES.md`](DATA_LICENSES.md) and `THIRD_PARTY_NOTICES.txt`. Updated DOI publication and journal submission are separate from this GitHub release.
