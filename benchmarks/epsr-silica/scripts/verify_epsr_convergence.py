#!/usr/bin/env python3
"""Verify the committed SiO2 EPSR convergence-pilot observations."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
result_root = case_root / "results/epsr-convergence-pilot"
raw_path = result_root / "raw-score-summary.json"
observed = tomllib.loads(
    (case_root / "expected/epsr-convergence-pilot-observed.toml").read_text()
)
raw = json.loads(raw_path.read_text())
parser = argparse.ArgumentParser()
parser.add_argument(
    "--strict-provenance",
    action="store_true",
    help="also require the recorded raw hash and machine-specific wall times",
)
args = parser.parse_args()


def close(label: str, actual: float, expected: float) -> None:
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=5e-10):
        raise SystemExit(f"{label}: {actual} != {expected}")


if args.strict_provenance:
    if hashlib.sha256(raw_path.read_bytes()).hexdigest() != observed["raw_score_sha256"]:
        raise SystemExit("EPSR convergence raw-score hash mismatch")

for case_name, expected_case in observed["cases"].items():
    actual_case = raw["cases"][case_name]
    for method, expected_method in expected_case.items():
        checkpoints = actual_case[method]
        for index, iteration in enumerate(observed["iterations"]):
            run = checkpoints[f"iter-{iteration:03d}"]
            reciprocal = run["common_reciprocal_distance_from_hidden_target"]
            actual = {
                "fit": math.hypot(
                    reciprocal["neutron_iq_rms"], reciprocal["xray_iq_rms"]
                )
                / math.sqrt(2.0),
                "rdf": run["structural_distance_from_hidden_target"][
                    "mean_partial_rdf_rms"
                ],
                "partial_sq": reciprocal["mean_partial_sq_rms"],
                "min_si_si_a": run["model_metrics"]["minimum_distance_a"]["Si-Si"],
            }
            for metric, value in actual.items():
                close(
                    f"{case_name}/{method}/{iteration}/{metric}",
                    value,
                    expected_method[metric][index],
                )
            if run["fit"]["wall_seconds"] <= 0.0:
                raise SystemExit(f"{case_name}/{method}/{iteration}: invalid wall time")
            if args.strict_provenance:
                close(
                    f"{case_name}/{method}/{iteration}/wall_s",
                    run["fit"]["wall_seconds"],
                    expected_method["wall_s"][index],
                )

for case_name in observed["cases"]:
    for iteration in observed["iterations"][:-1]:
        checkpoint = f"iter-{iteration:03d}"
        guarded = (
            result_root
            / case_name
            / "rsmith-epsr"
            / "seed-20260802"
            / checkpoint
            / "refined.xyz"
        )
        unguarded = (
            result_root
            / case_name
            / "rsmith-epsr-unconstrained"
            / "seed-20260802"
            / checkpoint
            / "refined.xyz"
        )
        if guarded.read_bytes() != unguarded.read_bytes():
            raise SystemExit(f"constraint sensitivity changed coordinates: {case_name}/{checkpoint}")

for case_name in observed["cases"]:
    pcof = (
        result_root
        / case_name
        / "native-epsr26-rminex-control"
        / "seed-20260802"
        / "iter-100"
        / "DTBsilica.pcof"
    )
    lines = pcof.read_text().splitlines()
    retained = {}
    for index, line in enumerate(lines[:-1]):
        fields = line.split()
        if len(fields) >= 2 and (fields[0], fields[1]) in {
            ("Si", "Si"),
            ("Si", "O"),
            ("O", "O"),
        }:
            retained[(fields[0], fields[1])] = float(lines[index + 1].split()[0])
    for pair, expected_value in {
        ("Si", "Si"): 2.0,
        ("Si", "O"): 1.35,
        ("O", "O"): 2.0,
    }.items():
        if not math.isclose(retained[pair], expected_value, rel_tol=0.0, abs_tol=5e-7):
            raise SystemExit(
                f"retained rminex {case_name}/{pair}: "
                f"{retained[pair]} != {expected_value}"
            )

    for iteration in observed["iterations"]:
        checkpoint = f"iter-{iteration:03d}"
        default = (
            result_root
            / case_name
            / "native-epsr26"
            / "seed-20260802"
            / checkpoint
            / "Cross.ato"
        )
        rminex = (
            result_root
            / case_name
            / "native-epsr26-rminex-control"
            / "seed-20260802"
            / checkpoint
            / "Cross.ato"
        )
        if default.read_bytes() != rminex.read_bytes():
            raise SystemExit(f"native rminex control changed coordinates: {case_name}/{checkpoint}")

for pair in observed["matched_fit"]:
    case = raw["cases"][pair["case"]]
    native = case["native-epsr26"][f"iter-{pair['native_iteration']:03d}"]
    rsmith = case["rsmith-epsr-unconstrained"][f"iter-{pair['rsmith_iteration']:03d}"]
    native_recip = native["common_reciprocal_distance_from_hidden_target"]
    rsmith_recip = rsmith["common_reciprocal_distance_from_hidden_target"]
    native_fit = math.hypot(native_recip["neutron_iq_rms"], native_recip["xray_iq_rms"]) / math.sqrt(2.0)
    rsmith_fit = math.hypot(rsmith_recip["neutron_iq_rms"], rsmith_recip["xray_iq_rms"]) / math.sqrt(2.0)
    close("matched fit difference", abs(native_fit - rsmith_fit), pair["absolute_fit_difference"])
    native_rdf = native["structural_distance_from_hidden_target"]["mean_partial_rdf_rms"]
    rsmith_rdf = rsmith["structural_distance_from_hidden_target"]["mean_partial_rdf_rms"]
    close("matched RDF difference", native_rdf - rsmith_rdf, pair["native_minus_rsmith_rdf"])
    wall_ratio = native["fit"]["wall_seconds"] / rsmith["fit"]["wall_seconds"]
    if wall_ratio <= 1.0:
        raise SystemExit("matched-fit rsmith run is not faster")
    if args.strict_provenance:
        close("matched wall ratio", wall_ratio, pair["native_over_rsmith_wall"])
    if rsmith_rdf >= native_rdf:
        raise SystemExit("matched-fit rsmith RDF recovery is not better")

print("EPSR convergence pilot: PASS (unconstrained parity, inert rminex control, four matched-fit pairs, iteration-100 stop reproduced)")
