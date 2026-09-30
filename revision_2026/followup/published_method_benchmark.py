"""Same-reference comparison with the census-block median-rate rule.

All references are known LOS rates withheld artificially inside the final
valid domain. This is a conditional methodological benchmark, not validation
of genuinely invalid cells, demographic locations, or physical hazard.
"""
from pathlib import Path
import hashlib
import json
import os

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter


STUDY = Path(os.environ["INSAR_STRENGTHENING_ROOT"])
OUT = Path(os.environ.get("INSAR_PAIRED_ROOT", Path(__file__).resolve().parent)) / "published_method_benchmark"
SRC = STUDY / "results" / "timeseries" / "matched_quality.npz"
SEED = 20260930
REPLICATES = 100
MASK_FRACTION = 0.30
BLOCK_PIXELS = 10
THRESHOLD = -5.0


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with np.load(SRC) as z:
        pop = z["population"].astype(np.float64)
        rate = z["full_vel_los"].astype(np.float64)
        valid = z["full_valid"].astype(bool)
        spread = z["full_vstd"].astype(np.float64)
    valid &= np.isfinite(rate) & np.isfinite(spread) & (pop >= 0)
    target = rate <= THRESHOLD
    h, w = pop.shape
    yy, xx = np.indices((h, w))
    rr, cc = yy.ravel(), xx.ravel()
    rng = np.random.default_rng(SEED)
    patterns = {
        "random": lambda: rng.standard_normal(pop.shape),
        "spatial_cluster": lambda: gaussian_filter(rng.standard_normal(pop.shape), 8),
        "high_velocity_spread": lambda: spread + 0.25 * gaussian_filter(rng.standard_normal(pop.shape), 8),
        "low_velocity_spread": lambda: -spread + 0.25 * gaussian_filter(rng.standard_normal(pop.shape), 8),
    }
    valid_ids = np.flatnonzero(valid)
    full_pop = float(pop[valid].sum())
    rows = []

    # Fixed geometry, one mask fraction, and one threshold keep the method
    # comparison focused. Random-number order matches actual_los_masking.py.
    for rep in range(REPLICATES):
        for mechanism, make_score in patterns.items():
            score = make_score()
            order = valid_ids[np.argsort(score.ravel()[valid_ids], kind="stable")]
            count = int(round(len(order) * MASK_FRACTION))
            mask = np.zeros(pop.size, dtype=bool)
            mask[order[-count:]] = True
            mask = mask.reshape(pop.shape)
            observed = valid & ~mask
            for dy, dx in ((0, 0), (BLOCK_PIXELS // 2, 0), (0, BLOCK_PIXELS // 2),
                           (BLOCK_PIXELS // 2, BLOCK_PIXELS // 2)):
                gy, gx = (yy + dy) // BLOCK_PIXELS, (xx + dx) // BLOCK_PIXELS
                labels = (gy * int(gx.max() + 1) + gx).ravel()
                n = int(labels.max()) + 1

                def agg(field):
                    return np.bincount(labels, weights=field.ravel(), minlength=n)

                missed_population = agg(pop * mask)
                true_by_unit = agg(pop * mask * target)
                truth = float(true_by_unit.sum())
                observed_population = agg(pop * observed)
                observed_class = agg(pop * observed * target)
                prevalence = np.divide(observed_class, observed_population,
                                       out=np.zeros(n), where=observed_population > 0) * missed_population

                # Published census-block median-rate assignment, adapted here
                # to fixed 10-pixel test units and the same hidden allocation.
                observed_ids = np.flatnonzero(observed.ravel())
                order_by_unit = np.lexsort((rate.ravel()[observed_ids], labels[observed_ids]))
                unit_ids = labels[observed_ids[order_by_unit]]
                sorted_rates = rate.ravel()[observed_ids[order_by_unit]]
                unique_ids, first, counts = np.unique(unit_ids, return_index=True, return_counts=True)
                lower = first + (counts - 1) // 2
                upper = first + counts // 2
                medians = 0.5 * (sorted_rates[lower] + sorted_rates[upper])
                median_class = np.zeros(n, dtype=bool)
                median_class[unique_ids] = medians <= THRESHOLD
                median_estimate = missed_population * median_class

                # Nearest observed pixel to each geometric unit centre; tied
                # centre pixels are averaged as in the paired-mask experiment.
                rmin = np.full(n, h, dtype=np.int32)
                rmax = np.full(n, -1, dtype=np.int32)
                cmin = np.full(n, w, dtype=np.int32)
                cmax = np.full(n, -1, dtype=np.int32)
                np.minimum.at(rmin, labels, rr)
                np.maximum.at(rmax, labels, rr)
                np.minimum.at(cmin, labels, cc)
                np.maximum.at(cmax, labels, cc)
                cy = (rmin + rmax) / 2.0
                cx = (cmin + cmax) / 2.0
                distance2 = (rr - cy[labels]) ** 2 + (cc - cx[labels]) ** 2
                min_distance2 = np.full(n, np.inf)
                np.minimum.at(min_distance2, labels, np.where(observed.ravel(), distance2, np.inf))
                ties = observed.ravel() & (distance2 == min_distance2[labels])
                n_ties = np.bincount(labels, weights=ties.astype(np.float64), minlength=n)
                class_ties = np.bincount(labels, weights=(ties & target.ravel()).astype(np.float64), minlength=n)
                nearest = np.divide(class_ties, n_ties, out=np.zeros(n), where=n_ties > 0) * missed_population

                for method, estimate in (
                    ("within_block_observed_population_prevalence", prevalence),
                    ("census_unit_observed_median_rate", median_estimate),
                    ("nearest_observed_pixel_to_unit_center", nearest),
                ):
                    signed_error = float(estimate.sum() - truth)
                    abs_local_error = float(np.abs(estimate - true_by_unit).sum())
                    rows.append({
                        "replicate": rep,
                        "mechanism": mechanism,
                        "missing_fraction": MASK_FRACTION,
                        "los_threshold_mm_yr": THRESHOLD,
                        "block_pixels": BLOCK_PIXELS,
                        "origin_y_pixels": dy,
                        "origin_x_pixels": dx,
                        "method": method,
                        "observed_valid_pixels": int(valid.sum()),
                        "domain_population_allocation_units": full_pop,
                        "hidden_population_allocation_units": float(missed_population.sum()),
                        "true_hidden_class_population_allocation_units": truth,
                        "estimated_hidden_class_population_allocation_units": float(estimate.sum()),
                        "signed_error_population_allocation_units": signed_error,
                        "signed_error_percentage_points_of_domain": 100.0 * signed_error / full_pop,
                        "absolute_block_error_population_allocation_units": abs_local_error,
                        "absolute_block_error_percentage_points_of_domain": 100.0 * abs_local_error / full_pop,
                    })
        if (rep + 1) % 20 == 0:
            print(f"replicates {rep + 1}/{REPLICATES}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "replicates.csv", index=False)
    summary = (df.groupby(["mechanism", "method"], as_index=False)
        .agg(replicates=("replicate", "nunique"), origins=("origin_x_pixels", "count"),
             mean_truth=("true_hidden_class_population_allocation_units", "mean"),
             mean_estimate=("estimated_hidden_class_population_allocation_units", "mean"),
             mean_signed_error_pp=("signed_error_percentage_points_of_domain", "mean"),
             q025_signed_error_pp=("signed_error_percentage_points_of_domain", lambda x: x.quantile(.025)),
             q975_signed_error_pp=("signed_error_percentage_points_of_domain", lambda x: x.quantile(.975)),
             mean_local_absolute_error_pp=("absolute_block_error_percentage_points_of_domain", "mean")))
    summary.to_csv(OUT / "summary.csv", index=False)
    metadata = {
        "source": str(SRC),
        "source_sha256": hashlib.sha256(SRC.read_bytes()).hexdigest(),
        "seed": SEED,
        "replicates_per_mechanism": REPLICATES,
        "missing_fraction": MASK_FRACTION,
        "threshold_mm_per_year": THRESHOLD,
        "block_size_pixels": BLOCK_PIXELS,
        "grid_origins": "four zero/half-block origins, matching the existing S14 design",
        "methods": ["within-block observed-population class prevalence",
                    "observed median LOS rate assigned to whole census-like unit",
                    "nearest observed pixel to unit centre"],
        "truth": "Known full-rate GHSL-model allocation hidden by artificial masks inside the final valid LOS domain.",
        "uncertainty": "2.5th and 97.5th percentiles span seeded masks and the four predeclared grid origins; they are design-sensitivity ranges, not confidence intervals.",
        "limitations": ["The fixed square units approximate census units but do not reproduce Fresno census block boundaries.",
                        "The reference uses a population allocation surface, not enumerated household locations.",
                        "No inference about genuinely invalid or permanently unobserved pixels, vertical motion, or hazard is made."],
    }
    (OUT / "method.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
