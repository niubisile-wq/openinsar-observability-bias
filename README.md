# Open InSAR observability-bias manuscript repository

The current scientific revision is **[Spatial aggregation alters population-weighted estimates of InSAR observation support](revision_2026/manuscript/manuscript.pdf)**.

## Current revision and reviewer reproduction

Start with **[revision_2026/README.md](revision_2026/README.md)** for installation, a data-based quick reproduction, the full analysis suite, provider downloads, numerical comparisons, and the limits of each experiment.

```sh
git clone https://github.com/niubisile-wq/openinsar-observability-bias.git
cd openinsar-observability-bias
python -m venv .venv
# Linux/macOS:
.venv/bin/python -m pip install -r revision_2026/requirements_reviewer_lock.txt
.venv/bin/python revision_2026/reproduce.py quick --workspace reproduction_output
```

For Windows, use the explicit PowerShell commands in the [reproduction guide](revision_2026/README.md#start-here).

The [v0.4.1 release](https://github.com/niubisile-wq/openinsar-observability-bias/releases/tag/v0.4.1) also provides a checksum-verified companion archive of the complete upstream ACS source files and the official Census archives used by current SI S25 (historical S20). These files can still be fetched directly with the documented downloader.

Release **[v0.6.0](https://github.com/niubisile-wq/openinsar-observability-bias/releases/tag/v0.6.0)** corrects missing-unit estimation and DWR annualization, adds complete raw cropped DWR windows and class tables, and delivers the reorganized manuscript/SI and reviewer responses. Corrected error, availability and bounds are reported together. The 15 affected tables have been rerun from a fresh workspace. See the [correction record](revision_2026/REVIEW_CORRECTIONS_v0.6.0.txt) and [old/new SI crosswalk](revision_2026/provenance/section_crosswalk.csv). Earlier releases remain historical and are superseded where stated.

The current **[v0.6.1](https://github.com/niubisile-wq/openinsar-observability-bias/releases/tag/v0.6.1) document update** synchronizes the 2 October title, abstract, introduction and writing revisions. It also embeds references in the submission source ZIPs, orders the seven upload figures correctly, and supplies the final response locations. Scientific code, inputs and numerical reference outputs are unchanged from v0.6.0.

- [Current clean manuscript](revision_2026/manuscript/manuscript.pdf)
- [Numbered manuscript for response locations](revision_2026/manuscript/manuscript_numbered.pdf)
- [Manuscript marked against the retained local baseline](revision_2026/manuscript/manuscript_marked_local_baseline.pdf)
- [Supplementary Information, S1–S30](revision_2026/manuscript/supplementary.pdf)
- [Point-by-point response](revision_2026/manuscript/response.pdf)
- [Reviewer-to-experiment map](revision_2026/REVIEWER_MAP.md)
- [Scientific reference outputs](revision_2026/reference_outputs/)
- [Frozen analysis inputs and provenance](revision_2026/data/)

The quick route reruns the bounds, masking benchmark, tolerance and ACS analyses and checks scientific results against frozen tables. `python revision_2026/reproduce.py additional` reruns all six added analyses and checks every new table, including block, replicate, tile and station rows. The standard suite combines old and new analyses from processed inputs. Windows/Linux CI executes the quick, DWR-rate, full LOS-masking and additional suites. The retained raw-processing code, pinned LiCSBAS commit and acquisition manifests document the separate upstream-data route.

## Version history and DOI

The earlier manuscript in `14_esin_strengthened_v1/` and the historical folders remain available as source history. **DOI [10.5281/zenodo.21444768](https://doi.org/10.5281/zenodo.21444768) archives the earlier v0.3.0 release only.** It does not cover the 2026 revision. Publication of an updated Zenodo DOI remains pending; use the revision GitHub tag/commit to identify the current code in the meantime.

This repository release does not submit the manuscript to a journal. Numerical reproduction verifies computations for defined inputs and does not provide independent geodetic or household-population validation.
