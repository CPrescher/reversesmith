#!/usr/bin/env python3
"""Summarize fixed-budget and achieved-fit-overlap SiO2 ensembles."""

from __future__ import annotations

import json
import math
import random
import statistics
import tomllib
from pathlib import Path


BOOTSTRAP_REPLICATES = 20_000
BOOTSTRAP_SEED = 20260803


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def distribution(values: list[float], rng: random.Random) -> dict:
    medians = []
    for _ in range(BOOTSTRAP_REPLICATES):
        sample = [values[rng.randrange(len(values))] for _ in values]
        medians.append(statistics.median(sample))
    return {
        "n": len(values),
        "median": statistics.median(values),
        "q1": percentile(values, 0.25),
        "q3": percentile(values, 0.75),
        "minimum": min(values),
        "maximum": max(values),
        "bootstrap_median_ci95": [
            percentile(medians, 0.025),
            percentile(medians, 0.975),
        ],
    }


def metrics(run: dict) -> dict[str, float]:
    reciprocal = run["common_reciprocal_distance_from_hidden_target"]
    structural = run["structural_distance_from_hidden_target"]
    fit = run["fit"]
    minimum_distance = run["model_metrics"]["minimum_distance_a"]
    neutron = reciprocal["neutron_iq_rms"]
    xray = reciprocal["xray_iq_rms"]
    result = {
        "combined_total_rms": math.sqrt((neutron * neutron + xray * xray) / 2.0),
        "neutron_rms": neutron,
        "xray_rms": xray,
        "mean_partial_sq_rms": reciprocal["mean_partial_sq_rms"],
        "mean_partial_rdf_rms": structural["mean_partial_rdf_rms"],
        "si_defect_fraction_error": structural["si_defect_fraction_error"],
        "o_defect_fraction_error": structural["o_defect_fraction_error"],
        "o_si_o_mean_error_degrees": structural["o_si_o_mean_error_degrees"],
        "si_o_si_mean_error_degrees": structural["si_o_si_mean_error_degrees"],
        "ring_total_variation": structural["ring_total_variation"],
        "minimum_oo_distance_a": minimum_distance["O-O"],
        "minimum_si_o_distance_a": minimum_distance["Si-O"],
        "minimum_si_si_distance_a": minimum_distance["Si-Si"],
        "wall_seconds": fit.get("sampling_wall_seconds", fit.get("wall_seconds")),
    }
    if fit.get("accepted_moves") is not None:
        result["accepted_moves"] = float(fit["accepted_moves"])
    return result


case_root = Path(__file__).resolve().parents[1]
root = case_root / "results/multiseed-comparison"
raw = json.loads((root / "raw-score-summary.json").read_text())
protocol = tomllib.loads(
    (case_root / "expected/multiseed-comparison.toml").read_text()
)
tolerance = protocol["outcomes"]["matched_fit"][
    "maximum_absolute_scalar_difference"
]
minimum_pairs = protocol["outcomes"]["matched_fit"]["minimum_pairs_for_claim"]
rng = random.Random(BOOTSTRAP_SEED)
summary = {
    "status": "completed_multiseed_cross_program_summary",
    "bootstrap_replicates": BOOTSTRAP_REPLICATES,
    "bootstrap_seed": BOOTSTRAP_SEED,
    "cases": {},
}

for case_name, methods in raw["cases"].items():
    method_metrics = {
        method: {seed: metrics(run) for seed, run in seeds.items()}
        for method, seeds in methods.items()
    }
    case_summary = {"fixed_budget": {}, "paired_same_seed": {}, "matched_fit": {}}
    for method, seeds in method_metrics.items():
        metric_names = sorted({name for values in seeds.values() for name in values})
        case_summary["fixed_budget"][method] = {
            name: distribution(
                [values[name] for values in seeds.values() if name in values], rng
            )
            for name in metric_names
        }

    for first, second in (
        ("rsmith-pace-w3", "rsmith-rmc"),
        ("rsmith-pace-w3", "rsmith-pedone-w3"),
    ):
        shared = sorted(set(method_metrics[first]) & set(method_metrics[second]))
        comparison = {}
        for name in (
            "combined_total_rms",
            "mean_partial_sq_rms",
            "mean_partial_rdf_rms",
            "ring_total_variation",
        ):
            differences = [
                method_metrics[first][seed][name] - method_metrics[second][seed][name]
                for seed in shared
            ]
            comparison[name + "_first_minus_second"] = distribution(differences, rng)
        case_summary["paired_same_seed"][f"{first}_vs_{second}"] = comparison

    method_names = sorted(method_metrics)
    for first_index, first in enumerate(method_names):
        first_ranked = sorted(
            method_metrics[first].items(), key=lambda item: item[1]["combined_total_rms"]
        )
        for second in method_names[first_index + 1 :]:
            second_ranked = sorted(
                method_metrics[second].items(),
                key=lambda item: item[1]["combined_total_rms"],
            )
            pairs = [
                (first_item, second_item)
                for first_item, second_item in zip(first_ranked, second_ranked)
                if abs(
                    first_item[1]["combined_total_rms"]
                    - second_item[1]["combined_total_rms"]
                )
                <= tolerance
            ]
            key = f"{first}_vs_{second}"
            matched = {
                "pair_count": len(pairs),
                "claim_supported": len(pairs) >= minimum_pairs,
            }
            if pairs:
                matched["median_absolute_fit_difference"] = statistics.median(
                    abs(a[1]["combined_total_rms"] - b[1]["combined_total_rms"])
                    for a, b in pairs
                )
                for name in ("mean_partial_sq_rms", "mean_partial_rdf_rms"):
                    matched[name + "_first_minus_second"] = distribution(
                        [a[1][name] - b[1][name] for a, b in pairs], rng
                    )
            case_summary["matched_fit"][key] = matched
    summary["cases"][case_name] = case_summary

(root / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
for case_name, case in summary["cases"].items():
    print(case_name)
    for method, values in case["fixed_budget"].items():
        print(
            f"  {method:20s} fit={values['combined_total_rms']['median']:.6f} "
            f"rdf={values['mean_partial_rdf_rms']['median']:.6f} "
            f"wall={values['wall_seconds']['median']:.3f}s"
        )
    supported = [
        name
        for name, values in case["matched_fit"].items()
        if values["claim_supported"]
    ]
    print("  fit-overlap:", ", ".join(supported) if supported else "none")
