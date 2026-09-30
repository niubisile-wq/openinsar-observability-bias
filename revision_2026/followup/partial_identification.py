"""Finite-domain allocation bounds for LOS classes beneath the final mask."""
from pathlib import Path
import hashlib
import json
import os

import numpy as np
import pandas as pd


STUDY = Path(os.environ["INSAR_STRENGTHENING_ROOT"])
OUT = Path(os.environ.get("INSAR_PAIRED_ROOT", Path(__file__).resolve().parent)) / "partial_identification"
SRC = STUDY / "results" / "timeseries" / "matched_quality.npz"
THRESHOLDS = (-2.0, -5.0, -10.0)


def main():
    with np.load(SRC) as z:
        pop = z["population"].astype(np.float64)
        rate = z["full_vel_los"].astype(np.float64)
        final_valid = z["full_valid"].astype(bool)
    allocation_domain = np.isfinite(pop) & (pop >= 0)
    rate_valid = allocation_domain & final_valid & np.isfinite(rate)
    invalid = allocation_domain & ~rate_valid
    total = float(pop[allocation_domain].sum())
    visible_population = float(pop[rate_valid].sum())
    unknown_population = float(pop[invalid].sum())
    rows = []
    for threshold in THRESHOLDS:
        visible_class = float(pop[rate_valid & (rate <= threshold)].sum())
        lower = visible_class
        upper = min(total, visible_class + unknown_population)
        rows.append({
            "los_threshold_mm_yr": threshold,
            "population_allocation_units_in_fixed_domain": total,
            "visible_population_allocation_units": visible_population,
            "unsupported_population_allocation_units": unknown_population,
            "visible_class_allocation_units": visible_class,
            "finite_domain_class_lower_bound_units": lower,
            "finite_domain_class_upper_bound_units": upper,
            "class_share_lower_pct": 100.0 * lower / total if total else 0.0,
            "class_share_upper_pct": 100.0 * upper / total if total else 0.0,
            "identification_width_pp": 100.0 * (upper - lower) / total if total else 0.0,
        })
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "bounds.csv", index=False)
    metadata = {
        "source": str(SRC),
        "source_sha256": hashlib.sha256(SRC.read_bytes()).hexdigest(),
        "population_surface": "GHSL model allocation; units are not verified residents",
        "support_definition": "The retained final full-stack LiCSBAS validity mask",
        "observed_rate": "Relative LOS mm per year; not vertical velocity or a hazard threshold",
        "bound": "For unknown invalid-domain class share theta in [0, 1], total allocation is E_visible + theta * U_invalid.",
        "valid_population_units": visible_population,
        "unsupported_population_units": unknown_population,
        "domain_population_units": total,
        "limitations": ["Bounds are allocation sensitivity, not motion estimates in invalid cells.",
                        "They are conditional on one population surface, one final mask, and a finite geographic domain.",
                        "This provides no probabilistic coverage or external population validation."],
    }
    (OUT / "method.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
