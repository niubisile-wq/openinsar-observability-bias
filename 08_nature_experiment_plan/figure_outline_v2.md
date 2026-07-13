# Figure Outline v2
Date: 2026-07-09

This version turns the paper into an execution map. Every panel below has a single conclusion, an exact source file, and a concrete export target.

## Figure 1. Why this paper exists
### Panel A
Conclusion: Open InSAR observability is not random missingness; it is a structured censoring problem.
Source: `02_benchmark_v0_1/benchmark_v0_1_report.md`
Export target: concept diagram with three terms only: `observability failure`, `monitoring debt`, `risk underestimation factor`.

### Panel B
Conclusion: The benchmark spans multiple regions, not a single case study.
Source: `02_benchmark_v0_1/benchmark_region_evidence_v0_1.csv`
Export target: region map or region matrix with Po, Chao Phraya, Indus, Rhone, Brantas, Rhine.

### Panel C
Conclusion: The project boundary is benchmarkable, not just rhetorical.
Source: `02_benchmark_v0_1/benchmark_dataset_inventory_v0_1.csv`
Export target: source-data map with 3 tiers: primary benchmark, secondary benchmark, extension layer.

## Figure 2. Chao Phraya lead case
### Panel A
Conclusion: Strong subsidence dominates the lead case.
Source: `03_exposure_closure/chao_phraya_area_weighted_exposure_censoring/chao_phraya_area_weighted_exposure_summary.csv`
Export target: VLM summary heatmap or classification map from the cell table.

### Panel B
Conclusion: The same lead case has measurable observability censoring.
Source: `03_exposure_closure/chao_phraya_area_weighted_exposure_censoring/chao_phraya_area_weighted_exposure_summary.csv`
Export target: hidden-fraction map from `chao_phraya_area_weighted_exposure_cells.csv`.

### Panel C
Conclusion: The hidden exposure is not small in people or built-up area.
Source: `03_exposure_closure/chao_phraya_area_weighted_exposure_censoring/chao_phraya_area_weighted_exposure_summary.csv`
Export target: paired bar chart for total / strong / hidden population and built-up km2.

### Panel D
Conclusion: Transport infrastructure also falls into the hidden-exposure zone.
Source: `03_exposure_closure/chao_phraya_osm_exposure_censoring/chao_phraya_osm_exposure_censoring_summary.csv`
Export target: road + railway stacked exposure bars.

## Figure 3. Exposure-to-risk translation
### Panel A
Conclusion: Visible-only exposure materially undercounts strong-motion exposure.
Source: `03_exposure_closure/chao_phraya_area_weighted_exposure_censoring/chao_phraya_area_weighted_exposure_summary.csv`
Export target: visible vs full exposure waterfall chart.

### Panel B
Conclusion: The underestimation factor is non-trivial across regions.
Source: `02_benchmark_v0_1/benchmark_region_evidence_v0_1.csv`
Export target: regional bar chart of `bias_or_minus_one_proxy` or `landcover_adjusted_or`.

### Panel C
Conclusion: Population and built-up missed under strong motion are the policy-relevant quantities.
Source: `03_exposure_closure/chao_phraya_area_weighted_exposure_censoring/chao_phraya_area_weighted_exposure_summary.csv`
Export target: dual-axis or paired bar chart for missed population and missed built-up fraction.

### Panel D
Conclusion: Roads and rail are a separate infrastructure exposure layer.
Source: `03_exposure_closure/chao_phraya_osm_exposure_censoring/chao_phraya_osm_exposure_censoring_summary.csv`
Export target: transport exposure bars with hidden fraction.

## Figure 4. Mechanism and robustness
### Panel A
Conclusion: The censoring signal survives reasonable observability thresholds.
Source: `03_exposure_closure/chao_phraya_robustness_grid/chao_phraya_robustness_grid.csv`
Export target: threshold x strong-threshold heatmap of odds ratios.

### Panel B
Conclusion: The signal survives block bootstrap correction for spatial dependence.
Source: `03_exposure_closure/chao_phraya_robustness_grid/chao_phraya_robustness_grid.md`
Export target: interval plot for OR and not-majority fractions.

### Panel C
Conclusion: The result remains robust for moderate and strong deformation thresholds.
Source: `03_exposure_closure/chao_phraya_robustness_grid/chao_phraya_robustness_grid.csv`
Export target: line plot of OR vs strong threshold, faceted by observability threshold.

### Panel D
Conclusion: The robustness check is still a screening model, not a final spatial model.
Source: `08_nature_experiment_plan/claim_collision_matrix.md`
Export target: small text callout or boxed note in the figure caption, not a plotted panel.

## Figure 5. Controls and transfers
### Panel A
Conclusion: The benchmark is not a single dominant land-cover artifact.
Source: `02_benchmark_v0_1/benchmark_region_evidence_v0_1.csv`
Export target: land-cover composition vs OR scatter.

### Panel B
Conclusion: The signal is consistent with product-lineage bias, not one lucky frame.
Source: `01_innovation_reports/创新点落地进度_2026-07-09.md`
Export target: schematic of product lineage and observability failure.

### Panel C
Conclusion: Japan Niigata confirms the extension path to non-China/non-Europe public InSAR products.
Source: `04_japan_licsbas_probe/h5_velocity_summary.json`
Export target: small summary panel with sign/unit audit callout.

### Panel D
Conclusion: Iran confirms the extension path to a nationwide public InSAR product.
Source: `05_iran_insar_probe/iran_insar_probe_report.md`
Export target: small summary panel with sign/unit audit callout.

## Figure 6. Transfer and scope limits
### Panel A
Conclusion: The method transfers to another public deformation stack.
Source: `04_japan_licsbas_probe/h5_velocity_summary.json`
Export target: timeline or velocity summary.

### Panel B
Conclusion: The method also transfers to the Iran nationwide map.
Source: `05_iran_insar_probe/iran_insar_probe_report.md`
Export target: rate distribution / coverage summary.

### Panel C
Conclusion: GNSS anchors are sparse but adequate as independent support.
Source: `02_benchmark_v0_1/benchmark_region_evidence_v0_1.csv`
Export target: anchor availability matrix.

### Panel D
Conclusion: Europe remains the densest upgrade path once EGMS credentials are available.
Source: `06_egms_query_pack/EGMS_API突围查询包_2026-07-09.md`
Export target: no-credential upgrade schematic.

## Extended Data
- Threshold sensitivity table: `03_exposure_closure/chao_phraya_robustness_grid/chao_phraya_robustness_grid.csv`
- Block bootstrap table: same file
- Area-weighted cell table: `03_exposure_closure/chao_phraya_area_weighted_exposure_censoring/chao_phraya_area_weighted_exposure_cells.csv`
- Transport cell table: `03_exposure_closure/chao_phraya_osm_exposure_censoring/chao_phraya_osm_exposure_censoring_summary.csv`
- Dataset inventory: `02_benchmark_v0_1/benchmark_dataset_inventory_v0_1.csv`

## Build order
1. Lock Fig. 1 and Fig. 2 first.
2. Render Fig. 3 and Fig. 4 from the current numeric tables.
3. Use Fig. 5 and Fig. 6 as scope-control and transfer panels.
4. Put all bootstrap tables and audits into Extended Data.

## Review risk
- Do not overstate Fig. 2 as a global result; it is the lead case.
- Do not make Fig. 3 look like a formal spatial model; it is a robustness screen.
- Do not move transfer claims from Fig. 6 into the main causal claim.
