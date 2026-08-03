#!/usr/bin/env python3
"""Summarize the preregistered EPSR control/start sensitivity matrix."""

from __future__ import annotations

import argparse
import json
import math
import tomllib
from functools import lru_cache
from pathlib import Path

import numpy as np


case_root = Path(__file__).resolve().parents[1]
result_root = case_root / "results/epsr-control-start-sensitivity"
protocol = tomllib.loads(
    (case_root / "expected/epsr-control-start-sensitivity.toml").read_text()
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
        difference = abs(
            native_points[native_index][0] - rsmith_points[rsmith_index][0]
        )
        if difference <= tolerance:
            count, cost, pairs = solve(native_index + 1, rsmith_index + 1)
            candidates.append(
                (count + 1, cost + difference, ((native_index, rsmith_index), *pairs))
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
    }


def bootstrap_mean_ci(values, seed: int, samples: int = 20_000):
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(values), size=(samples, len(values)))
    return [
        float(value)
        for value in np.quantile(values[indices].mean(axis=1), [0.025, 0.975])
    ]


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


expected_seeds = [
    f"seed-{seed}" for seed in protocol["design"]["refinement_seeds"]
]
expected_arms = [arm["name"] for arm in protocol["arms"]]
tolerance = float(protocol["matching"]["tolerance"])
summary = {
    "status": "epsr_control_start_sensitivity_summarized",
    "protocol": "epsr-control-start-sensitivity.toml",
    "matching": protocol["matching"]["algorithm"],
    "cases": {},
}

for case_index, (case_name, arms) in enumerate(sorted(raw["cases"].items())):
    if not args.allow_incomplete and sorted(arms) != sorted(expected_arms):
        raise SystemExit(f"{case_name}: incomplete sensitivity-arm set")
    case_summary = {"arms": {}}
    for arm_index, arm_name in enumerate(expected_arms):
        if arm_name not in arms:
            continue
        methods = arms[arm_name]
        native_seeds = methods.get("native-epsr26", {})
        rsmith_seeds = methods.get("rsmith-epsr-unconstrained", {})
        available_seeds = sorted(set(native_seeds) & set(rsmith_seeds))
        if not args.allow_incomplete and available_seeds != expected_seeds:
            raise SystemExit(
                f"{case_name}/{arm_name}: complete seeds {available_seeds} "
                f"!= expected {expected_seeds}"
            )
        arm_summary = {"seeds": {}}
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
                raise SystemExit(
                    f"{case_name}/{arm_name}/{seed_name}: no fit-matched checkpoints"
                )
            median_fit = float(
                np.median([record["fit_midpoint"] for record in records])
            )
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
            native_final_minimum = native_final["model_metrics"][
                "minimum_distance_a"
            ]["Si-Si"]
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
            arm_summary["seeds"][seed_name] = {
                "accepted_pair_count": len(records),
                "primary": primary,
                "best_fit_end": records[0],
                "worst_fit_end": records[-1],
                "all_pairs": records,
                "native_iteration_100_min_si_si_a": native_final_minimum,
                "native_iteration_100_below_threshold": failed,
            }
        if not available_seeds:
            continue
        random_seed = 20260803 + 100 * case_index + 10 * arm_index
        rdf = np.asarray(primary_rdf_differences, dtype=float)
        paired_rdf = distribution(primary_rdf_differences, random_seed)
        win_count = int(np.sum(rdf > 0.0))
        complete = len(available_seeds) == len(expected_seeds)
        robust = (
            paired_rdf["median"] > 0.0 and win_count >= 4 if complete else None
        )
        strong = (
            robust and paired_rdf["bootstrap_mean_95_ci"][0] > 0.0
            if complete
            else None
        )
        sensitive = (
            paired_rdf["mean"] <= 0.0 or win_count < 3 if complete else None
        )
        arm_summary["aggregate"] = {
            "matched_seed_count": len(available_seeds),
            "accepted_pair_count_range": [
                min(accepted_pair_counts),
                max(accepted_pair_counts),
            ],
            "rsmith_lower_primary_rdf_count": win_count,
            "native_iteration_100_close_contact_count": native_failures,
            "primary_absolute_fit_difference": distribution(
                primary_fit_differences, random_seed + 1
            ),
            "native_minus_rsmith_primary_rdf": paired_rdf,
            "primary_rdf_unpaired_bootstrap_mean_difference_95_ci": (
                bootstrap_unpaired_mean_difference_ci(
                    primary_native_rdf, primary_rsmith_rdf, random_seed + 2
                )
            ),
            "native_minus_rsmith_primary_partial_sq": distribution(
                primary_sq_differences, random_seed + 3
            ),
            "sensitivity": {
                "native_minus_rsmith_best_fit_end_rdf": distribution(
                    best_fit_rdf_differences, random_seed + 4
                ),
                "native_minus_rsmith_worst_fit_end_rdf": distribution(
                    worst_fit_rdf_differences, random_seed + 5
                ),
            },
            "decision": {
                "robust": robust,
                "strong": strong,
                "sensitive": sensitive,
            },
        }
        case_summary["arms"][arm_name] = arm_summary
    summary["cases"][case_name] = case_summary

complete_decisions = [
    arm["aggregate"]["decision"]
    for case in summary["cases"].values()
    for arm in case["arms"].values()
    if arm["aggregate"]["decision"]["sensitive"] is not None
]
summary["overall_decision"] = {
    "evaluated_case_arms": len(complete_decisions),
    "retain_bounded_superiority_claim": bool(complete_decisions)
    and not any(decision["sensitive"] for decision in complete_decisions),
    "all_case_arms_robust": bool(complete_decisions)
    and all(decision["robust"] for decision in complete_decisions),
    "all_case_arms_strong": bool(complete_decisions)
    and all(decision["strong"] for decision in complete_decisions),
}

(result_root / "matched-fit-summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n"
)
for case_name, case in summary["cases"].items():
    for arm_name, arm in case["arms"].items():
        aggregate = arm["aggregate"]
        rdf = aggregate["native_minus_rsmith_primary_rdf"]
        print(
            f"{case_name}/{arm_name}: rsmith lower RDF in "
            f"{aggregate['rsmith_lower_primary_rdf_count']}/"
            f"{aggregate['matched_seed_count']}; mean native-rsmith "
            f"{rdf['mean']:.6f}, 95% CI "
            f"[{rdf['bootstrap_mean_95_ci'][0]:.6f}, "
            f"{rdf['bootstrap_mean_95_ci'][1]:.6f}]; "
            f"decision={aggregate['decision']}"
        )
