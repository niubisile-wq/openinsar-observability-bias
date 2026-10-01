# Reviewer comments and reproducible evidence

The exact point-by-point responses are in [`manuscript/response.pdf`](manuscript/response.pdf) and the editable `response.tex`. Each response locates changes in the supplied [`manuscript/manuscript_numbered.pdf`](manuscript/manuscript_numbered.pdf); the precise locator text is also retained in [`provenance/final_document_locations.json`](provenance/final_document_locations.json). All original editor/reviewer comment blocks are retained verbatim. This table provides direct experiment entry points; it does not claim that a reviewer has accepted the answer.

| Comment | Revision response | Code / evidence |
| --- | --- | --- |
| R1.1: unvalidated hidden/exposed-population claim | Withdraw the unvalidated physical count; define a continuous allocation/support estimand; distinguish internal processing checks from external validation. | E1 `allocation`; S19 `paired` / `direct-geometry`; S13 `partial-identification`; S26 `acs-population-benchmark`; retained matched-processing summaries/manifests. |
| R1.2: coarse binary-majority support and consequences | Preserve fractional support and exposure–support cross-moments; report mask, scale, origin and counterexample sensitivities. | `allocation`, `synthetic-masking`, `dwr-masking`, `observed-los-masking`, `published-method-benchmark`. |
| R1.3: limitations of GHSL as population evidence | Add WorldPop, epoch/common-total/common-coverage analyses, Census/CA-POP within-block normalization, and an ACS administrative-unit benchmark with exact aggregate MOE. | `population-analysis`, `cross-sensitivity`, `fine-population`, `spatial-dependence`, `acs-population-benchmark`; S25/S26 limitations remain explicit. |
| R1.4: numerical, presentation and policy issues | Correct units, figures, source attribution and dates; remove unsupported risk/policy-priority interpretations. | Current manuscript/Supplementary Information, input manifests, output/hash verification. |
| R2.1: unsupported allegation of misuse in prior studies | Remove the allegation; distinguish estimands fairly and compare an actual published census-unit median-rate rule with within-unit prevalence and center sampling. | `published-method-benchmark` (S14), `followup`, Fresno DWR/Census reference tables and source acquisition notes. |
| R2.2: ESIN rationale and competing approaches | Explain allocation through the established covariance identity; retain counterexamples and evaluate alternatives on matched inputs. | E1/E3; S8 follow-up factorial outputs; S25 model sensitivity; S14 published-method comparison; S9 `reporting-tolerance`. |
| R2.3: Iran source/date confusion | Exclude the exploratory Iran branch from quantitative applications, correct the 2014–2020 product attribution, and make no completed new/old-product comparison claim. | Revised text and response; no invented Iran benchmark is supplied. |

Additional experiments strengthen several comments simultaneously. E4 reports variation conditional on the fixed 224-pair pool and preserves catalogue/network provenance; legacy-28 versus expanded-224 differences also change composition. The spatial Moran's I results are descriptive for the fixed Fresno domain. S9 tolerances are illustrative choices, not safety standards.

## Added S15, S20, S21, S27, S29, S30 evidence

| Added evidence | Relevant comments | What it adds and its limit |
| --- | --- | --- |
| S15 `additional-real-blocks` | R1.1, R1.2, R2.2 | Actual irregular Census blocks, 30 repeated masks, contiguous holes, estimator availability, common-computable errors and whole-known-domain bounds. Does not establish genuinely missing motion. |
| S20 `additional-dependence` | R1.2, R1.3, R2.2 | Buffered training-only spatial recalibration and 225 shared-date deletions. Released-grid recalibration is not original model-training CV; date deletions do not re-invert velocities. |
| S21 `additional-quality` | R1.1, R1.2, R2.2 | All 18 training-only quality rules, 44 excluded edges, coverage and internal weighted errors. These dependent displacement residuals are not external annual-rate accuracy. |
| S29 `additional-regions` | Editor contribution/scope, R2.1, R2.2 | Two regions fixed before outcome calculation; retain all 96 designs and small effects. Published averaged-coherence support differs from pair support and limits magnitude transfer. |
| S27 `additional-residential` | R1.3, R2.2 | Independent municipal land-use eligibility with strict/broad definitions. Household truth is unavailable; the common Census domain falls below the exploratory 90% coverage diagnostic. |
| S30 `additional-gnss` | R1.1, R2.1 | Archived independent observations compared with external LOS products; exclude all control sites and retain all radii, reference frames and large residuals. Source-date conflicts and temporal-model mismatch limit the comparison. |

**Evidence limit:** S30 supplies an archived external-product GNSS consistency comparison, including the negative Davis-Sacramento result. There is still no matched independent GNSS/leveling accuracy validation of our Bangkok inversion. Census/ACS and municipal-use benchmarks do not establish household locations or deformation in truly unsupported cells. These limits are explicit in the manuscript and replies.

**Editorial code-access item:** This GitHub revision supplies versioned code, inputs, tables and runnable instructions. A new DOI deposit and journal portal submission are separate remaining operations.
