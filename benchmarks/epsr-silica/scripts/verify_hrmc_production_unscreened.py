#!/usr/bin/env python3
"""Verify the single-seed joint-acceptance Pedone/GAP/PACE comparison."""

from __future__ import annotations

import json
import math
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
expected = tomllib.loads(
    (case_root / "expected/hrmc-production-unscreened-smoke.toml").read_text()
)
score = json.loads(
    (case_root / "results/hrmc-production-unscreened/score-summary.json").read_text()
)
pure_score = json.loads(
    (case_root / "results/cross-recovery/score-summary.json").read_text()
)
guards = expected["guards"]
failures: list[str] = []


def progress(run: dict, kind: str, pure_fit: dict) -> float:
    start = run["fit"][f"start_{kind}_rms"]
    final = run["fit"][f"refined_{kind}_rms"]
    pure_final = pure_fit[f"refined_{kind}_rms"]
    return (start - final) / (start - pure_final)


for case_name, runs in score["cases"].items():
    pure_case = pure_score["cases"][case_name]
    pure_fit = pure_case["fits"]["rsmith-joint"]
    pure_rdf = pure_case["structural_distance_from_hidden_target"]["rsmith_joint"]
    for run_name in guards["required_runs"]:
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

    pace = runs[f"pace-w{guards['pace_selected_weight']:g}"]
    for kind in ("neutron", "xray"):
        value = progress(pace, kind, pure_fit)
        if value < 1.0 - 1.0e-3:
            failures.append(f"{case_name}/pace-w3: {kind} progress {value} < 0.999")
    if (
        pace["structural_distance_from_hidden_target"]["mean_partial_rdf_rms"]
        >= pure_rdf["mean_partial_rdf_rms"]
    ):
        failures.append(f"{case_name}/pace-w3: hidden-target RDF did not improve")

    for gap_name, pace_name in (
        ("gap-w0p1", "pace-w3"),
        ("gap-w0p3", "pace-w10"),
        ("gap-w0p4", "pace-w30"),
    ):
        speedup = runs[gap_name]["fit"]["wall_seconds"] / runs[pace_name]["fit"][
            "wall_seconds"
        ]
        if speedup < guards["pace_minimum_speedup_over_gap"]:
            failures.append(
                f"{case_name}/{pace_name}: PACE/GAP speedup {speedup} below guard"
            )

if failures:
    raise SystemExit("HRMC production-smoke verification failed:\n- " + "\n- ".join(failures))
print("passed single-seed joint-acceptance Pedone/GAP/PACE checks")
