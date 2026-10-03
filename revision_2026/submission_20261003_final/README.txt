Final revision document snapshot v0.6.3, 3 October 2026

This directory identifies the final current manuscript, Supplementary Information, response and cover
letter. The five PDF files have 103 pages in total: clean manuscript 28, line-numbered
manuscript 28, SI 37, response 9 and cover letter 1.
The line-numbered manuscript is the same scientific text with line numbers, not a marked comparison.
The main article has six figures and two tables. SI contains S1-S31, 20 figures and 29 tables.

Contents
- documents/: five current PDFs.
- source/LaTeX_Source.zip: portable editable sources, bibliography, class/style files and required figures.
- figures/Figure_1.pdf through Figure_6.pdf: main figures in current manuscript order.
- audits/Local_Validation_Audit.zip: unchanged CUSV/CUUT local feasibility inputs and portable replay.
- audits/SBKK_Support_Audit.zip: unchanged filtered local LiCSBAS support window and portable replay.
- audits/Reproducibility_Records.zip: current path index and retained calculation/literature records.
- audits/Figure_Source_Data.zip: unchanged figure-source member bytes with their provenance.
- supporting_records/: retained numerical summary and four-row published Bangkok comparison/source audit.
- section_output_index.json: S1-S31 and historical code names mapped to released result locations.
- submission_manifest.json: SHA-256, file lengths and PDF page counts for this snapshot.

Version scope
The scientific code, frozen inputs and numerical references retain the corrected v0.6.0 baseline,
commit 0dd41575840917186aa504677d08996e34004b44. This is a document update, not a new numerical analysis.
The preceding v0.6.2 snapshot at ../submission_20261003/ and its tag remain unchanged; its commit is
e493b61624fdc08a89e7f7f9809709481196d043. The older manuscript/ directory is the v0.6.1 writing snapshot.

Reproduction
Follow ../README.md and its locked Python 3.10 environment for numerical reproduction. For the two
portable local audits, extract the respective ZIP and follow its included instructions. Source ZIP
compilation and file integrity checks are separate from numerical reproduction. Local support and
feasibility checks do not establish independent matched Bangkok geodetic accuracy or household truth.
The source thesis PDF, complete original interferogram stack and private author review logs are not
included. Source/provider terms remain as recorded in ../DATA_LICENSES.md and ../THIRD_PARTY_NOTICES.txt.
