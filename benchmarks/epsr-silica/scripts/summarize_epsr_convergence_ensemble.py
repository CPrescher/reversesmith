#!/usr/bin/env python3
"""Summarize preregistered fit-matched EPSR convergence-ensemble outputs."""

from __future__ import annotations

import argparse
import json
import math
import tomllib
from functools import lru_cache
from pathlib import Path

import numpy as np


case_root = Path(__file__).resolve().parents[1]
result_root = case_root / "results/epsr-convergence-ensemble"
protocol = tomllib.loads(
    (case_root / "expected/epsr-convergence-ensemble.toml").read_text()
)
raw = json.loads((result_root / "raw-score-summary.json").read_text())
parser = argparse.ArgumentParser()
parser.add_argument("--allow-incomplete", action="store_true")
args = parser.parse_args()


def fit_scalar(record: dict) -> float:
    reciprocal = record["common_reciprocal_distance_from_hidden_target"]
    return math.hypot(
        reciprocal["neutron_iq_rms"], reciprocal["xray_iq_rms"]
    ) / math.sqrt(2.0)


def checkpoint_number(name: str) -> int:
    return int(name.removeprefix("iter-"))


def ordered_fit_matches(native: dict, rsmith: dict, tolerance: float):
    native_points = sorted(
        ((fit_scalar(value), key, value) for key, value in native.items()),
        key=lambda item: (item[0], checkpoint_number(item[1])),
    )
    rsmith_points = sorted(
        ((fit_scalar(value), key, value) for key, value in rsmith.items()),
        key=lambda item: (item[0], checkpoint_number(item[1])),
    )

    @lru_cache(maxsize=None)
    def solve(native_index: int, rsmith_index: int):
        if native_index == len(native_points) or rsmith_index == len(rsmith_points):
            return 0, 0.0, ()
        candidates = [
            solve(native_index + 1, rsmith_index),
            solve(native_index, rsmith_index + 1),
        ]
        delta = abs(
            native_points[native_index][0] - rsmith_points[rsmith_index][0]
        )
        if delta <= tolerance:
            count, cost, pairs = solve(native_index + 1, rsmith_index + 1)
            candidates.append(
                (
                    count + 1,
                    cost + delta,
                    ((native_index, rsmith_index), *pairs),
                )
            )
        return min(candidates, key=lambda item: (-item[0], item[1], item[2]))

    _, _, indices = solve(0, 0)
    return [
        (native_points[native_index], rsmith_points[rsmith_index])
        for native_index, rsmith_index in indices
    ]


def pair_record(native_point, rsmith_point):
    native_fit, native_name, native = native_point
    rsmith_fit, rsmith_name, rsmith = rsmith_point
    native_structure = native["structural_distance_from_hidden_target"]
    rsmith_structure = rsmith["structural_distance_from_hidden_target"]
    native_reciprocal = native["common_reciprocal_distance_from_hidden_target"]
    rsmith_reciprocal = rsmith["common_reciprocal_distance_from_hidden_target"]
    return {
        "native_checkpoint": checkpoint_number(native_name),
        "rsmith_checkpoint": checkpoint_number(rsmith_name),
        "native_fit": native_fit,
        "rsmith_fit": rsmith_fit,
        "absolute_fit_difference": abs(native_fit - rsmith_fit),
        "fit_midpoint": (native_fit + rsmith_fit) / 2.0,
        "native_rdf": native_structure["mean_partial_rdf_rms"],
        "rsmith_rdf": rsmith_structure["mean_partial_rdf_rms"],
        "native_minus_rsmith_rdf": native_structure["mean_partial_rdf_rms"]
        - rsmith_structure["mean_partial_rdf_rms"],
        "native_partial_sq": native_reciprocal["mean_partial_sq_rms"],
        "rsmith_partial_sq": rsmith_reciprocal["mean_partial_sq_rms"],
        "native_minus_rsmith_partial_sq": native_reciprocal["mean_partial_sq_rms"]
        - rsmith_reciprocal["mean_partial_sq_rms"],
        "native_min_si_si_a": native["model_metrics"]["minimum_distance_a"][
            "Si-Si"
        ],
        "rsmith_min_si_si_a": rsmith["model_metrics"]["minimum_distance_a"][
            "Si-Si"
        ],
    }


def bootstrap_mean_ci(values, seed: int, samples: int = 20_000):
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(values), size=(samples, len(values)))
    means = values[indices].mean(axis=1)
    return [float(value) for value in np.quantile(means, [0.025, 0.975])]


def bootstrap_unpaired_mean_difference_ci(
    native_values, rsmith_values, seed: int, samples: int = 20_000
):
    native_values = np.asarray(native_values, dtype=float)
    rsmith_values = np.asarray(rsmith_values, dtype=float)
    rng = np.random.default_rng(seed)
    native_indices = rng.integers(
        0, len(native_values), size=(samples, len(native_values))
    )
    rsmith_indices = rng.integers(
        0, len(rsmith_values), size=(samples, len(rsmith_values))
    )
    differences = native_values[native_indices].mean(axis=1) - rsmith_values[
        rsmith_indices
    ].mean(axis=1)
    return [float(value) for value in np.quantile(differences, [0.025, 0.975])]


def distribution(values, seed: int):
    values = np.asarray(values, dtype=float)
    return {
        "values": [float(value) for value in values],
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "interquartile_range": [
            float(value) for value in np.quantile(values, [0.25, 0.75])
        ],
        "bootstrap_mean_95_ci": bootstrap_mean_ci(values, seed),
    }


expected_seeds = [f"seed-{seed}" for seed in protocol["design"]["seeds"]]
tolerance = float(protocol["outcomes"]["matching_tolerance"])
summary = {
    "status": "epsr_dense_matched_fit_ensemble_summarized",
    "protocol": "epsr-convergence-ensemble.toml",
    "matching": "maximum-cardinality monotone matching, then minimum total absolute fit difference",
    "cases": {},
}

for case_index, (case_name, methods) in enumerate(sorted(raw["cases"].items())):
    native_seeds = methods.get("native-epsr26", {})
    rsmith_seeds = methods.get("rsmith-epsr-unconstrained", {})
    available_seeds = sorted(set(native_seeds) & set(rsmith_seeds))
    if not args.allow_incomplete and available_seeds != expected_seeds:
        raise SystemExit(
            f"{case_name}: complete seeds {available_seeds} != expected {expected_seeds}"
        )
    case_summary = {"seeds": {}}
    primary_rdf_differences = []
    primary_sq_differences = []
    primary_fit_differences = []
    primary_native_rdf = []
    primary_rsmith_rdf = []
    best_fit_rdf_differences = []
    worst_fit_rdf_differences = []
    accepted_pair_counts = []
    native_failures = 0
    for seed_name in available_seeds:
        matches = ordered_fit_matches(
            native_seeds[seed_name], rsmith_seeds[seed_name], tolerance
        )
        records = [pair_record(native, rsmith) for native, rsmith in matches]
        records.sort(
            key=lambda item: (
                item["fit_midpoint"],
                item["absolute_fit_difference"],
                item["native_checkpoint"],
                item["rsmith_checkpoint"],
            )
        )
        if not records:
            raise SystemExit(f"{case_name}/{seed_name}: no fit-matched checkpoints")
        median_fit = float(np.median([record["fit_midpoint"] for record in records]))
        primary = min(
            records,
            key=lambda item: (
                abs(item["fit_midpoint"] - median_fit),
                item["absolute_fit_difference"],
                item["native_checkpoint"],
                item["rsmith_checkpoint"],
            ),
        )
        native_final = native_seeds[seed_name]["iter-100"]
        native_final_minimum = native_final["model_metrics"]["minimum_distance_a"][
            "Si-Si"
        ]
        failed = native_final_minimum < float(
            protocol["outcomes"]["native_failure_threshold_si_si_a"]
        )
        native_failures += int(failed)
        accepted_pair_counts.append(len(records))
        primary_rdf_differences.append(primary["native_minus_rsmith_rdf"])
        primary_sq_differences.append(primary["native_minus_rsmith_partial_sq"])
        primary_fit_differences.append(primary["absolute_fit_difference"])
        primary_native_rdf.append(primary["native_rdf"])
        primary_rsmith_rdf.append(primary["rsmith_rdf"])
        best_fit_rdf_differences.append(records[0]["native_minus_rsmith_rdf"])
        worst_fit_rdf_differences.append(records[-1]["native_minus_rsmith_rdf"])
        case_summary["seeds"][seed_name] = {
            "accepted_pair_count": len(records),
            "primary": primary,
            "best_fit_end": records[0],
            "worst_fit_end": records[-1],
            "all_pairs": records,
            "native_iteration_100_min_si_si_a": native_final_minimum,
            "native_iteration_100_below_threshold": failed,
        }
    rdf = np.asarray(primary_rdf_differences, dtype=float)
    case_summary["aggregate"] = {
        "matched_seed_count": len(available_seeds),
        "accepted_pair_count_range": [
            min(accepted_pair_counts),
            max(accepted_pair_counts),
        ],
        "rsmith_lower_primary_rdf_count": int(np.sum(rdf > 0.0)),
        "native_iteration_100_close_contact_count": native_failures,
        "primary_absolute_fit_difference": distribution(
            primary_fit_differences, 20260793 + case_index
        ),
        "native_minus_rsmith_primary_rdf": distribution(
            primary_rdf_differences, 20260803 + case_index
        ),
        "primary_rdf_unpaired_bootstrap_mean_difference_95_ci": bootstrap_unpaired_mean_difference_ci(
            primary_native_rdf,
            primary_rsmith_rdf,
            20260808 + case_index,
        ),
        "native_minus_rsmith_primary_partial_sq": distribution(
            primary_sq_differences, 20260813 + case_index
        ),
        "sensitivity": {
            "native_minus_rsmith_best_fit_end_rdf": distribution(
                best_fit_rdf_differences, 20260823 + case_index
            ),
            "native_minus_rsmith_worst_fit_end_rdf": distribution(
                worst_fit_rdf_differences, 20260833 + case_index
            ),
        },
    }
    summary["cases"][case_name] = case_summary

(result_root / "matched-fit-summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n"
)
for case_name, case in summary["cases"].items():
    aggregate = case["aggregate"]
    rdf = aggregate["native_minus_rsmith_primary_rdf"]
    print(
        f"{case_name}: {aggregate['matched_seed_count']} seeds, "
        f"rsmith lower RDF in {aggregate['rsmith_lower_primary_rdf_count']}, "
        f"native close contacts in "
        f"{aggregate['native_iteration_100_close_contact_count']}; "
        f"mean native-rsmith RDF {rdf['mean']:.6f}, "
        f"95% CI [{rdf['bootstrap_mean_95_ci'][0]:.6f}, "
        f"{rdf['bootstrap_mean_95_ci'][1]:.6f}]"
    )
