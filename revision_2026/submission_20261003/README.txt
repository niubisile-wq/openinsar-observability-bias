2026-10-03 submission snapshot

This directory versions the current submission documents and the two local audit archives alongside
the established scientific code and reference tables. The snapshot is separate from the manuscript/
folder's historical documents. The final submission_manifest.json identifies each delivered file by
SHA-256 and records PDF page counts after the final build; use those counts rather than earlier prose.

Contents
- documents/: the current clean manuscript, manuscript with line numbers, Supplementary Information,
  response and cover letter. The line-numbered manuscript is not a marked/highlighted comparison.
- source/LaTeX_Source.zip: the corresponding portable source archive for the five current documents.
- figures/Figure_1.pdf through figures/Figure_6.pdf: six main figures in current manuscript order.
- audits/Figure_Source_Data.zip: source data for the current figures with retained provenance.
- audits/Local_Validation_Audit.zip: CUSV/CUUT local GNSS feasibility, containing-pixel/3x3 series, date
  matching, source hashes and portable replay. It supplies no independent local accuracy score.
- audits/SBKK_Support_Audit.zip: the filtered LiCSBAS local input window, 104 acquisition dates, original
  support summary/time series, parent-file provenance and portable replay for the SBKK support check.
  This checks product-specific support and does not compare GNSS vertical rates directly with LOS.
- audits/Reproducibility_Records.zip: the section/path index, retained supporting records and a short
  path guide for delivery with the submission. This is a records package in addition to the two audits.
- supporting_records/: the retained numerical-check summary and four-row published leveling/PS fact
  table with its source audit; the source thesis PDF is not redistributed.
- section_output_index.json: S1-S31 to real public reference paths, historical path replacements,
  selected ledger scope, and the two additional audit archives.
- submission_manifest.json: the final file inventory, hashes and PDF page counts for this snapshot.

The main text has six figures and two tables. The SI has 31 sections, 20 figures and 29 tables.
The current five document PDFs total 101 pages: clean manuscript 27, line-numbered manuscript 27,
SI 37, response 9 and cover letter 1. These counts identify the supplied documents, not a claim
that all numerical experiments or every raw processing step were rerun for this document update.

Version scope
The numerical reference baseline is immutable v0.6.0, commit
0dd41575840917186aa504677d08996e34004b44. It corrected availability handling and annualized DWR intervals.
The immutable v0.6.1 snapshot (e3ae2c20c553ef8a26d2fc94577095fb5b2c659a) was a 2 October writing update;
its scientific code, inputs and numerical references stayed compatible with v0.6.0. The 3 October
submission snapshot adds the current documents, evidence-path index, retained supporting records and
local audit packages. Its release is v0.6.2; publication metadata records the final commit. Historical
tagged document snapshots must not be substituted for the current files in documents/.

Reproduction
For the established numerical results, follow ../README.md and run ../reproduce.py in its locked
Python 3.10 environment. Stable script names retain older section numbers; use section_output_index.json.
For S3 and the added local-feasibility paragraph of S17, extract the relevant audit ZIP into a fresh
directory and follow its own README and requirements. These audits are additional to the frozen
scientific-release reruns. PDF compilation, file integrity and numerical replay are distinct checks.
The complete original interferogram stack and independent common-datum geodetic accuracy validation
are not supplied by these small audit packages.

Source attribution
Third-party inputs retain provider terms and their included provenance. ../DATA_LICENSES.md,
../THIRD_PARTY_NOTICES.txt and the audit-specific notices describe their scope. The v0.4.1 and v0.5.0
companion input assets remain available at their historical releases; use the corrected current
reference tables, not superseded v0.5.0 derived class or zero-fill masking results.
