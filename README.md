# Open InSAR observability-bias manuscript repository

The current scientific revision is **[Quantifying how spatial allocation and sampling alter population-weighted InSAR support summaries](revision_2026/manuscript/manuscript.pdf)**.

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

The [v0.4.1 release](https://github.com/niubisile-wq/openinsar-observability-bias/releases/tag/v0.4.1) also provides a checksum-verified companion archive of the complete upstream ACS source files and the official Census archives used by S20. These files can still be fetched directly with the documented downloader.

Release [v0.5.1](https://github.com/niubisile-wq/openinsar-observability-bias/releases/tag/v0.5.1) integrates six new experiments, figures, Supplementary S25-S30 and regenerated response locators. Its document asset includes six final PDFs and two self-contained source archives verified by independent recompilation. Science is fixed at [v0.5.0](https://github.com/niubisile-wq/openinsar-observability-bias/releases/tag/v0.5.0), commit `f70124f1cdc11e43848ce97c01f49fad6f7efb62`; this includes 21 new numerical output comparisons, bundled portable inputs and 44.7-MB/586.3-MB companion assets. The earlier S1-S24 baseline and complete ACS/Census source asset remain available at v0.4.1.

- [Current clean manuscript](revision_2026/manuscript/manuscript.pdf)
- [Numbered manuscript for response locations](revision_2026/manuscript/manuscript_numbered.pdf)
- [Manuscript marked against the retained local baseline](revision_2026/manuscript/manuscript_marked_local_baseline.pdf)
- [Supplementary Information, S1–S30](revision_2026/manuscript/supplementary.pdf)
- [Point-by-point response](revision_2026/manuscript/response.pdf)
- [Reviewer-to-experiment map](revision_2026/REVIEWER_MAP.md)
- [Scientific reference outputs](revision_2026/reference_outputs/)
- [Frozen analysis inputs and provenance](revision_2026/data/)

The quick route reruns S21–S24 and checks scientific results against frozen tables. `python revision_2026/reproduce.py additional` reruns all six S25-S30 analyses and checks every new table, including block, replicate, tile and station rows. The standard suite combines old and new analyses from processed inputs. Windows/Linux CI executes the quick and additional suites. The retained raw-processing code, pinned LiCSBAS commit and acquisition manifests document the separate upstream-data route.

## Version history and DOI

The earlier manuscript in `14_esin_strengthened_v1/` and the historical folders remain available as source history. **DOI [10.5281/zenodo.21444768](https://doi.org/10.5281/zenodo.21444768) archives the earlier v0.3.0 release only.** It does not cover the 2026 revision. Publication of an updated Zenodo DOI remains pending; use the revision GitHub tag/commit to identify the current code in the meantime.

This repository release does not submit the manuscript to a journal. Numerical reproduction verifies computations for defined inputs and does not provide independent geodetic or household-population validation.
