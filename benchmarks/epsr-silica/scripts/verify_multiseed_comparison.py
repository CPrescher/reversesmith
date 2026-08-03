#!/usr/bin/env python3
"""Verify the committed multi-seed SiO2 summary against local run products."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
result_root = case_root / "results/multiseed-comparison"
observed_path = case_root / "expected/multiseed-comparison-observed.toml"
observed = tomllib.loads(observed_path.read_text())
summary_path = result_root / "summary.json"
raw_path = result_root / "raw-score-summary.json"
parser = argparse.ArgumentParser()
parser.add_argument(
    "--strict-provenance",
    action="store_true",
    help="also require hashes and machine-specific wall times from the recorded run",
)
args = parser.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(label: str, actual: float, expected: float) -> None:
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=5e-10):
        raise SystemExit(f"{label}: {actual} != {expected}")


if args.strict_provenance:
    if sha256(raw_path) != observed["raw_score_sha256"]:
        raise SystemExit("raw score summary hash mismatch")
    if sha256(summary_path) != observed["summary_sha256"]:
        raise SystemExit("statistical summary hash mismatch")

summary = json.loads(summary_path.read_text())
metric_map = {
    "fit": "combined_total_rms",
    "neutron": "neutron_rms",
    "xray": "xray_rms",
    "partial_sq": "mean_partial_sq_rms",
    "rdf": "mean_partial_rdf_rms",
}
for case_name, expected_case in observed["cases"].items():
    actual_case = summary["cases"][case_name]
    for method, expected_method in expected_case["methods"].items():
        actual_method = actual_case["fixed_budget"][method]
        for short_name, summary_name in metric_map.items():
            close(
                f"{case_name}/{method}/{short_name}",
                actual_method[summary_name]["median"],
                expected_method[short_name],
            )
        if actual_method["wall_seconds"]["median"] <= 0.0:
            raise SystemExit(f"{case_name}/{method}: invalid wall time")
        if args.strict_provenance:
            close(
                f"{case_name}/{method}/wall_s",
                actual_method["wall_seconds"]["median"],
                expected_method["wall_s"],
            )
        for short_name, summary_name in (("fit", "combined_total_rms"), ("rdf", "mean_partial_rdf_rms")):
            close(
                f"{case_name}/{method}/{short_name}_q1",
                actual_method[summary_name]["q1"],
                expected_method[f"{short_name}_q1"],
            )
            close(
                f"{case_name}/{method}/{short_name}_q3",
                actual_method[summary_name]["q3"],
                expected_method[f"{short_name}_q3"],
            )

    paired = actual_case["paired_same_seed"]
    expected_paired = expected_case["paired"]
    for comparator, prefix in (
        ("rsmith-pace-w3_vs_rsmith-rmc", "pace_minus_rmc"),
        ("rsmith-pace-w3_vs_rsmith-pedone-w3", "pace_minus_pedone"),
    ):
        actual = paired[comparator]["mean_partial_rdf_rms_first_minus_second"]
        close(f"{case_name}/{prefix}/median", actual["median"], expected_paired[f"{prefix}_rdf_median"])
        for index, endpoint in enumerate(expected_paired[f"{prefix}_rdf_ci95"]):
            close(f"{case_name}/{prefix}/ci95/{index}", actual["bootstrap_median_ci95"][index], endpoint)

for case in summary["cases"].values():
    for method in case["fixed_budget"].values():
        for key in ("minimum_si_o_distance_a", "minimum_oo_distance_a", "minimum_si_si_distance_a"):
            if method[key]["minimum"] <= 0.0:
                raise SystemExit(f"invalid minimum distance: {key}")

print("multi-seed SiO2 comparison: PASS (120 endpoints, committed medians and paired intervals reproduced)")
