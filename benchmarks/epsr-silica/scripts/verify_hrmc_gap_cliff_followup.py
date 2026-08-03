#!/usr/bin/env python3
"""Verify the bounded GAP acceptance-cliff follow-up."""

from __future__ import annotations

import json
import math
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
expected = tomllib.loads(
    (case_root / "expected/hrmc-gap-cliff-followup-smoke.toml").read_text()
)
score = json.loads(
    (case_root / "results/hrmc-weight-sweep/score-summary.json").read_text()
)
guards = expected["guards"]
failures: list[str] = []
first_failing_has_case_failure = False


def label(weight: float) -> str:
    return f"gap-w{weight:g}".replace(".", "p")


def progress(run: dict, kind: str, pure: dict) -> float:
    start = run["fit"][f"start_{kind}_rms"]
    final = run["fit"][f"refined_{kind}_rms"]
    pure_final = pure["fit"][f"refined_{kind}_rms"]
    return (start - final) / (start - pure_final)


for case_name, runs in score["cases"].items():
    pure = runs["pure-rmc"]
    for weight in guards["required_weights"]:
        run_name = label(weight)
        if run_name not in runs:
            failures.append(f"{case_name}: missing {run_name}")
            continue
        run = runs[run_name]
        if run["fit"].get("attempted_moves") != guards["attempted_moves"]:
            failures.append(f"{case_name}/{run_name}: wrong move count")
        wall = run["fit"].get("wall_seconds", 0.0)
        if not math.isfinite(wall) or wall <= 0.0:
            failures.append(f"{case_name}/{run_name}: invalid wall time")
        for pair, threshold in guards["minimum_distance_a"].items():
            value = run["model_metrics"]["minimum_distance_a"][pair]
            if value + 1.0e-9 < threshold:
                failures.append(
                    f"{case_name}/{run_name}: {pair} minimum {value} < {threshold}"
                )

    passing = runs[label(guards["highest_passing_weight"])]
    for kind in ("neutron", "xray"):
        value = progress(passing, kind, pure)
        if value < guards["weight_0p4_minimum_fit_progress_fraction"]:
            failures.append(f"{case_name}/gap-w0p4: {kind} progress {value} too small")
    baseline_energy = runs["gap-w0p001"]["fit"]["energy_change_ev_per_atom"]
    passing_energy = passing["fit"]["energy_change_ev_per_atom"]
    reduction = 1.0 - passing_energy / baseline_energy
    if (
        reduction
        < guards["weight_0p4_minimum_energy_drift_reduction_fraction_vs_0p001"]
    ):
        failures.append(f"{case_name}/gap-w0p4: energy-drift reduction too small")

    failing = runs[label(guards["first_failing_weight"])]
    if any(
        progress(failing, kind, pure) < guards["minimum_fit_progress_fraction"]
        for kind in ("neutron", "xray")
    ):
        first_failing_has_case_failure = True

if not first_failing_has_case_failure:
    failures.append("gap-w0p5: expected a progress failure in at least one direction")

if failures:
    raise SystemExit(
        "GAP cliff-follow-up verification failed:\n- " + "\n- ".join(failures)
    )
print("passed bounded GAP acceptance-cliff follow-up checks")
