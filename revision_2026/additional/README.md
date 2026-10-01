All six additional experiments are executable from the public checkout. Frozen inputs are bundled in `../data/additional/`; expected tables are in `../reference_outputs/additional/`. No private path, API key, account, or GIS query is required for these numerical reruns.

From the repository root, install the locked reviewer environment, then run:

```text
python -m pip install -r revision_2026/requirements_reviewer_lock.txt
python revision_2026/reproduce.py additional --workspace reproduction_output
```

This executes every calculation and compares **all 21 new CSV tables**, including replicate, block, tile and station tables, against the archived results. Numerical tolerances are the same as the established reviewer comparator (tightened to 1e-9 absolute for proportions). It is not a checksum-only test.

Individual modes: `additional-real-blocks` (S25), `additional-dependence` (S26), `additional-quality` (S27), `additional-regions` (S28), `additional-residential` (S29), `additional-gnss` (S30). `standard` also includes these six modes. The existing `quick` suite retains S21-S24. Windows and Linux CI execute both suites.

To regenerate the scientific figures after a rerun:

```text
python revision_2026/additional/make_additional_figures.py --results reproduction_output/additional --output reproduction_output/additional_figures
```

`ASSETS.json` records the checksums and URLs of two companion release assets. The 44.7-MB analysis-input ZIP duplicates the bundled portable inputs. The 586.3-MB selected-source ZIP supplies complete A064/A137 velocity/coherence maps, completed provider notebooks, all selected municipal land-use query geometry, GHSL, DWR rate sources and metadata. Third-party terms and attribution are stated in its README and provenance records. The complete ACS/Census sources remain available with v0.4.1. Original interferogram stacks and other provider-controlled inputs follow the main reviewer guide; neither a new daily-GNSS fit nor every statewide ARIA time series is claimed.

The `provenance_scripts/` folder retains acquisition and preprocessing scripts as an audit trail. These optional upstream steps require their stated original sources and local environment; their historical default paths are not used by the portable six-mode launcher. `protocol_v1.json` freezes the new-region, masking and quality designs; spatial recalibration is explicitly exploratory, not original population-model cross-validation.

Interpret the added evidence at its actual level:

* S25 masks known published rates on actual irregular Census blocks. Undefined blocks are retained. Small conditional errors do not identify motion in missing blocks.
* S26 tests training-only spatial recalibration and deletion of shared acquisition dates. It does not refit velocities after every date deletion or assume independent pairs.
* S27 uses training-only quality filters and all 44 excluded-edge residuals. Residual millimeters are internal displacement errors, not mm/year external velocity accuracy.
* S28 retains all 96 designs across two fixed regions; their small effects limit transfer of the primary magnitude. Published averaged coherence differs from pair-level support.
* S29 uses independent municipal land-use eligibility, not household counts. The initial 90% coverage diagnostic fails: the joint domain covers 89.13-89.22% of complete Census population. Unmapped classification and non-2020 vintage remain explicit.
* S30 compares archived independent GNSS observations with a published external LOS field. Reference/control sites are excluded from validation; exact native source windows reproduce all three sampling radii. The archived GNSS fits are linear, while InSAR includes seasonal/step terms. Large residuals and the lack of improvement in Davis-Sacramento are retained. These data do not validate the Bangkok inversion or true missing pixels.

Sangha et al. (2026): DOI [10.1029/2026EA005214](https://doi.org/10.1029/2026EA005214), archived data [10.5281/zenodo.19493073](https://doi.org/10.5281/zenodo.19493073), CC BY 4.0. The April archived data version is fixed; dates come from completed notebook outputs and the publication, with stale configuration dates preserved as conflicting metadata. NGL source coordinates: `https://geodesy.unr.edu/NGLStationPages/DataHoldings.txt`. Municipal source: City of Fresno, public Existing Land Use FeatureServer layer 19. Archived MintPy implementation excerpts retain the accompanying GPL license, not the repository's code license.
