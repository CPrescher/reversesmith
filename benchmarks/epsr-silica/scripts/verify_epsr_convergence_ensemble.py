#!/usr/bin/env python3
"""Verify the committed dense EPSR convergence-ensemble observations."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
result_root = case_root / "results/epsr-convergence-ensemble"
raw_path = result_root / "raw-score-summary.json"
summary_path = result_root / "matched-fit-summary.json"
protocol = tomllib.loads(
    (case_root / "expected/epsr-convergence-ensemble.toml").read_text()
)
observed = tomllib.loads(
    (case_root / "expected/epsr-convergence-ensemble-observed.toml").read_text()
)
raw = json.loads(raw_path.read_text())
summary = json.loads(summary_path.read_text())
parser = argparse.ArgumentParser()
parser.add_argument("--strict-provenance", action="store_true")
args = parser.parse_args()


def close(label: str, actual: float, expected: float) -> None:
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=5e-12):
        raise SystemExit(f"{label}: {actual} != {expected}")


if args.strict_provenance:
    for path, key in (
        (raw_path, "raw_score_sha256"),
        (summary_path, "matched_fit_summary_sha256"),
    ):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != observed[key]:
            raise SystemExit(f"{path.name}: hash {digest} != {observed[key]}")

expected_seeds = [f"seed-{seed}" for seed in protocol["design"]["seeds"]]
expected_checkpoints = {
    "native-epsr26": {
        f"iter-{iteration:03d}"
        for iteration in protocol["methods"]["native_epsr26"]["checkpoints"]
    },
    "rsmith-epsr-unconstrained": {
        f"iter-{iteration:03d}"
        for iteration in protocol["methods"]["rsmith_epsr_unconstrained"][
            "checkpoints"
        ]
    },
}

for case_name, expected_case in observed["cases"].items():
    raw_case = raw["cases"][case_name]
    for method, checkpoints in expected_checkpoints.items():
        if sorted(raw_case[method]) != expected_seeds:
            raise SystemExit(f"{case_name}/{method}: incomplete seed set")
        for seed_name, runs in raw_case[method].items():
            if set(runs) != checkpoints:
                raise SystemExit(f"{case_name}/{method}/{seed_name}: incomplete checkpoints")

    aggregate = summary["cases"][case_name]["aggregate"]
    if aggregate["matched_seed_count"] != expected_case["matched_seed_count"]:
        raise SystemExit(f"{case_name}: wrong matched-seed count")
    if aggregate["accepted_pair_count_range"] != [
        expected_case["accepted_pair_count"],
        expected_case["accepted_pair_count"],
    ]:
        raise SystemExit(f"{case_name}: wrong matched-pair count")
    if (
        aggregate["rsmith_lower_primary_rdf_count"]
        != expected_case["rsmith_lower_primary_rdf_count"]
    ):
        raise SystemExit(f"{case_name}: wrong primary RDF win count")
    if (
        aggregate["native_iteration_100_close_contact_count"]
        != expected_case["native_iteration_100_close_contact_count"]
    ):
        raise SystemExit(f"{case_name}: wrong native close-contact count")

    primary_fit = aggregate["primary_absolute_fit_difference"]
    close(
        f"{case_name}/maximum primary fit difference",
        max(primary_fit["values"]),
        expected_case["maximum_primary_absolute_fit_difference"],
    )
    if max(primary_fit["values"]) > protocol["outcomes"]["matching_tolerance"]:
        raise SystemExit(f"{case_name}: primary pair exceeds fit tolerance")

    checks = (
        (
            aggregate["native_minus_rsmith_primary_rdf"],
            "primary_rdf_mean_native_minus_rsmith",
            "primary_rdf_median_native_minus_rsmith",
            "primary_rdf_bootstrap_mean_95_ci",
        ),
        (
            aggregate["native_minus_rsmith_primary_partial_sq"],
            "primary_partial_sq_mean_native_minus_rsmith",
            None,
            "primary_partial_sq_bootstrap_mean_95_ci",
        ),
        (
            aggregate["sensitivity"]["native_minus_rsmith_best_fit_end_rdf"],
            "best_fit_end_rdf_mean_native_minus_rsmith",
            None,
            "best_fit_end_rdf_bootstrap_mean_95_ci",
        ),
        (
            aggregate["sensitivity"]["native_minus_rsmith_worst_fit_end_rdf"],
            "worst_fit_end_rdf_mean_native_minus_rsmith",
            None,
            "worst_fit_end_rdf_bootstrap_mean_95_ci",
        ),
    )
    for distribution, mean_key, median_key, interval_key in checks:
        close(f"{case_name}/{mean_key}", distribution["mean"], expected_case[mean_key])
        if median_key:
            close(
                f"{case_name}/{median_key}",
                distribution["median"],
                expected_case[median_key],
            )
        for actual, expected in zip(
            distribution["bootstrap_mean_95_ci"], expected_case[interval_key]
        ):
            close(f"{case_name}/{interval_key}", actual, expected)

    for actual, expected in zip(
        aggregate["primary_rdf_unpaired_bootstrap_mean_difference_95_ci"],
        expected_case["primary_rdf_unpaired_bootstrap_mean_difference_95_ci"],
    ):
        close(f"{case_name}/unpaired RDF bootstrap interval", actual, expected)

    recorded_minima = [
        seed["native_iteration_100_min_si_si_a"]
        for seed in summary["cases"][case_name]["seeds"].values()
    ]
    for actual, expected in zip(
        recorded_minima, expected_case["native_iteration_100_min_si_si_a"]
    ):
        close(f"{case_name}/native iteration-100 minimum", actual, expected)

print(
    "EPSR convergence ensemble: PASS "
    "(10 seeds, matched fit, hidden structure and stability reproduced)"
)
