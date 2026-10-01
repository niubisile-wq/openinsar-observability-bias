# Reviewer comments and reproducible evidence

The exact point-by-point responses are in [`manuscript/response.pdf`](manuscript/response.pdf) and the editable `response.tex`. Each response locates changes in the supplied [`manuscript/manuscript_numbered.pdf`](manuscript/manuscript_numbered.pdf); the precise locator text is also retained in [`provenance/reviewer_final_locators.json`](provenance/reviewer_final_locators.json). All original editor/reviewer comment blocks are retained verbatim. This table provides direct experiment entry points; it does not claim that a reviewer has accepted the answer.

| Comment | Revision response | Code / evidence |
| --- | --- | --- |
| R1.1: unvalidated hidden/exposed-population claim | Withdraw the unvalidated physical count; define a continuous allocation/support estimand; distinguish internal processing checks from external validation. | E1 `allocation`; S13 `paired` / `direct-geometry`; S21 `partial-identification`; S24 `acs-population-benchmark`; retained matched-processing summaries/manifests. |
| R1.2: coarse binary-majority support and consequences | Preserve fractional support and exposure–support cross-moments; report mask, scale, origin and counterexample sensitivities. | `allocation`, `synthetic-masking`, `dwr-masking`, `observed-los-masking`, `published-method-benchmark`. |
| R1.3: limitations of GHSL as population evidence | Add WorldPop, epoch/common-total/common-coverage analyses, Census/CA-POP within-block normalization, and a coarse ACS benchmark with exact aggregate MOE. | `population-analysis`, `cross-sensitivity`, `fine-population`, `spatial-dependence`, `acs-population-benchmark`; S20/S24 limitations remain explicit. |
| R1.4: numerical, presentation and policy issues | Correct units, figures, source attribution and dates; remove unsupported risk/policy-priority interpretations. | Current manuscript/Supplementary Information, input manifests, output/hash verification. |
| R2.1: unsupported allegation of misuse in prior studies | Remove the allegation; distinguish estimands fairly and compare an actual published census-unit median-rate rule with within-unit prevalence and center sampling. | `published-method-benchmark` (S22), `followup`, Fresno DWR/Census reference tables and source acquisition notes. |
| R2.2: ESIN rationale and competing approaches | Explain allocation through the established covariance identity; retain counterexamples and evaluate alternatives on matched inputs. | E1/E3; S15 follow-up factorial outputs; S20 model sensitivity; S22 published-method comparison; S23 `reporting-tolerance`. |
| R2.3: Iran source/date confusion | Exclude the exploratory Iran branch from quantitative applications, correct the 2014–2020 product attribution, and make no completed new/old-product comparison claim. | Revised text and response; no invented Iran benchmark is supplied. |

Additional experiments strengthen several comments simultaneously. E4 reports variation conditional on the fixed 224-pair pool and preserves catalogue/network provenance; legacy-28 versus expanded-224 differences also change composition. The spatial Moran's I results are descriptive for the fixed Fresno domain. S23 tolerances are illustrative choices, not safety standards.

**Evidence limit:** There is still no matched independent GNSS/leveling validation. Census/ACS benchmarks do not establish household locations or deformation in truly unsupported cells. These limits are part of the response, rather than missing outputs being described as completed validation.

**Editorial code-access item:** This GitHub revision supplies versioned code, inputs, tables and runnable instructions. A new DOI deposit and journal portal submission are separate remaining operations.
