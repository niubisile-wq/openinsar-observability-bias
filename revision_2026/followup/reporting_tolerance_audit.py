"""Evaluate declared reporting tolerances against retained covariance bounds."""
from pathlib import Path
import json
import os

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "reference_outputs" / "followup" / "fresno_population_allocation" / "scale_method_sensitivity.csv"
OUT = Path(os.environ.get("INSAR_PAIRED_ROOT", Path(__file__).resolve().parent)) / "reporting_tolerance"
TOLERANCES_PP = (0.5, 1.0, 2.0)


def main():
    frame = pd.read_csv(SOURCE)
    rows = []
    keys = ["population_model", "scale_m"]
    for (model, scale), group in frame.groupby(keys, sort=True):
        if len(group) != 4:
            raise ValueError(f"Expected four origins for {model} at {scale} m, got {len(group)}")
        for tolerance in TOLERANCES_PP:
            bounds = group["uniform_local_cauchy_bound_pp"].abs()
            local_l1 = group["uniform_local_L1_pp"].abs()
            net = group["uniform_minus_native_pp"].abs()
            rows.append({
                "population_model": model,
                "scale_m": int(scale),
                "declared_error_tolerance_pp": tolerance,
                "origins": len(group),
                "origins_certified_by_bound": int((bounds <= tolerance).sum()),
                "all_origins_certified": bool((bounds <= tolerance).all()),
                "max_bound_pp": float(bounds.max()),
                "max_local_L1_pp": float(local_l1.max()),
                "max_absolute_net_effect_pp": float(net.max()),
                "minimum_top_decile_overlap": float(group["uniform_top_decile_overlap"].min()),
            })
    OUT.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "tolerance_summary.csv", index=False)
    meta = {
        "input": str(SOURCE),
        "tolerances_pp": list(TOLERANCES_PP),
        "scope": "Finite Fresno domain, two 100-m allocation surfaces normalized to identical 2020 Census block totals, DWR footprint support, four grid origins at each scale.",
        "interpretation": "If the computed Cauchy-Schwarz upper bound is at or below the user-declared tolerance for all tested origins, the signed aggregate uniform-minus-native difference is bounded by that tolerance for the tested inputs. A failed bound is inconclusive; it does not establish that a coarse result is unacceptable.",
        "not_included": ["probability or confidence statements", "cross-region inference", "household location truth", "claims about physical deformation"],
    }
    (OUT / "method.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
