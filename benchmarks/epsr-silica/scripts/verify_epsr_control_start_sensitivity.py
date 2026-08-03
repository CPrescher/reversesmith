#!/usr/bin/env python3
"""Verify committed EPSR control/start sensitivity observations."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
result_root = case_root / "results/epsr-control-start-sensitivity"
raw_path = result_root / "raw-score-summary.json"
summary_path = result_root / "matched-fit-summary.json"
protocol = tomllib.loads(
    (case_root / "expected/epsr-control-start-sensitivity.toml").read_text()
)
observed = tomllib.loads(
    (case_root / "expected/epsr-control-start-sensitivity-observed.toml").read_text()
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

expected_seeds = [
    f"seed-{seed}" for seed in protocol["design"]["refinement_seeds"]
]
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
expected_arms = {arm["name"] for arm in protocol["arms"]}
endpoint_counts = {method: 0 for method in expected_checkpoints}
for case_name, arms in raw["cases"].items():
    if set(arms) != expected_arms:
        raise SystemExit(f"{case_name}: incomplete arm set")
    for arm_name, methods in arms.items():
        for method, checkpoints in expected_checkpoints.items():
            seeds = methods.get(method, {})
            if sorted(seeds) != expected_seeds:
                raise SystemExit(f"{case_name}/{arm_name}/{method}: incomplete seeds")
            for seed_name, runs in seeds.items():
                if set(runs) != checkpoints:
                    raise SystemExit(
                        f"{case_name}/{arm_name}/{method}/{seed_name}: "
                        "incomplete checkpoints"
                    )
                endpoint_counts[method] += len(runs)

if endpoint_counts["native-epsr26"] != observed["native_endpoint_count"]:
    raise SystemExit("wrong native endpoint count")
if endpoint_counts["rsmith-epsr-unconstrained"] != observed["rsmith_endpoint_count"]:
    raise SystemExit("wrong rsmith endpoint count")

strong_count = 0
sensitive_count = 0
for expected in observed["observations"]:
    label = f"{expected['case']}/{expected['arm']}"
    aggregate = summary["cases"][expected["case"]]["arms"][expected["arm"]][
        "aggregate"
    ]
    if aggregate["matched_seed_count"] != len(expected_seeds):
        raise SystemExit(f"{label}: wrong matched-seed count")
    if aggregate["accepted_pair_count_range"] != expected[
        "accepted_pair_count_range"
    ]:
        raise SystemExit(f"{label}: wrong accepted-pair range")
    if aggregate["rsmith_lower_primary_rdf_count"] != expected[
        "rsmith_lower_primary_rdf_count"
    ]:
        raise SystemExit(f"{label}: wrong RDF win count")
    if aggregate["native_iteration_100_close_contact_count"] != expected[
        "native_close_contact_count"
    ]:
        raise SystemExit(f"{label}: wrong close-contact count")
    maximum_fit_difference = max(
        aggregate["primary_absolute_fit_difference"]["values"]
    )
    close(
        f"{label}/maximum primary fit difference",
        maximum_fit_difference,
        expected["maximum_primary_fit_difference"],
    )
    if maximum_fit_difference > protocol["matching"]["tolerance"]:
        raise SystemExit(f"{label}: primary pair exceeds fit tolerance")
    rdf = aggregate["native_minus_rsmith_primary_rdf"]
    close(
        f"{label}/primary RDF mean",
        rdf["mean"],
        expected["primary_rdf_mean_native_minus_rsmith"],
    )
    close(
        f"{label}/primary RDF median",
        rdf["median"],
        expected["primary_rdf_median_native_minus_rsmith"],
    )
    for actual, wanted in zip(
        rdf["bootstrap_mean_95_ci"], expected["primary_rdf_bootstrap_mean_95_ci"]
    ):
        close(f"{label}/paired interval", actual, wanted)
    for actual, wanted in zip(
        aggregate["primary_rdf_unpaired_bootstrap_mean_difference_95_ci"],
        expected["primary_rdf_unpaired_bootstrap_mean_difference_95_ci"],
    ):
        close(f"{label}/unpaired interval", actual, wanted)
    for decision in ("robust", "strong", "sensitive"):
        if aggregate["decision"][decision] is not expected[decision]:
            raise SystemExit(f"{label}: wrong {decision} decision")
    strong_count += int(aggregate["decision"]["strong"])
    sensitive_count += int(aggregate["decision"]["sensitive"])

overall = summary["overall_decision"]
if overall["evaluated_case_arms"] != observed["evaluated_case_arms"]:
    raise SystemExit("wrong evaluated case-arm count")
if strong_count != observed["strong_case_arms"]:
    raise SystemExit("wrong strong case-arm count")
if sensitive_count != observed["sensitive_case_arms"]:
    raise SystemExit("wrong sensitive case-arm count")
if overall["retain_bounded_superiority_claim"] is not observed[
    "retain_bounded_superiority_claim"
]:
    raise SystemExit("wrong overall claim decision")

print(
    "EPSR control/start sensitivity: PASS "
    "(1,260 endpoints, 12 matched case-arms, frozen decisions reproduced)"
)
