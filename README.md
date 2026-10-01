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

Release [v0.4.3](https://github.com/niubisile-wq/openinsar-observability-bias/releases/tag/v0.4.3) finalizes the manuscript and response wording, evidence limits, source attribution, actual line references and consistent manuscript variants. It makes no changes to scientific inputs, figures or results; [v0.4.1](https://github.com/niubisile-wq/openinsar-observability-bias/releases/tag/v0.4.1), commit `5bf98f19abbed981989f9c97ef9b118eea524c5b`, remains the fixed scientific reproduction package.

- [Current clean manuscript](revision_2026/manuscript/manuscript.pdf)
- [Numbered manuscript for response locations](revision_2026/manuscript/manuscript_numbered.pdf)
- [Manuscript marked against the retained local baseline](revision_2026/manuscript/manuscript_marked_local_baseline.pdf)
- [Supplementary Information, S1–S24](revision_2026/manuscript/supplementary.pdf)
- [Point-by-point response](revision_2026/manuscript/response.pdf)
- [Reviewer-to-experiment map](revision_2026/REVIEWER_MAP.md)
- [Scientific reference outputs](revision_2026/reference_outputs/)
- [Frozen analysis inputs and provenance](revision_2026/data/)

The quick route reruns S21–S24 and checks scientific results against frozen tables. The wider suite replays the main allocation, masking, sampling, population-model and paired/LOS analyses from processed inputs. The retained raw-processing code, pinned LiCSBAS commit and acquisition manifests document the separate upstream-data route.

## Version history and DOI

The earlier manuscript in `14_esin_strengthened_v1/` and the historical folders remain available as source history. **DOI [10.5281/zenodo.21444768](https://doi.org/10.5281/zenodo.21444768) archives the earlier v0.3.0 release only.** It does not cover the 2026 revision. Publication of an updated Zenodo DOI remains pending; use the revision GitHub tag/commit to identify the current code in the meantime.

This repository release does not submit the manuscript to a journal. Numerical reproduction verifies computations for defined inputs and does not provide independent geodetic or household-population validation.
