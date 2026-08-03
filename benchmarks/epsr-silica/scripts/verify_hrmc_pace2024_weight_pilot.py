#!/usr/bin/env python3
"""Verify the bounded Erhard-2024 PACE weight calibration."""

from __future__ import annotations

import json
import math
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
expected = tomllib.loads(
    (case_root / "expected/hrmc-pace2024-weight-pilot-smoke.toml").read_text()
)
score = json.loads(
    (case_root / "results/hrmc-pace2024-weight-pilot/score-summary.json").read_text()
)
pure_score = json.loads(
    (case_root / "results/hrmc-weight-sweep/score-summary.json").read_text()
)
guards = expected["guards"]
failures: list[str] = []


def label(weight: float) -> str:
    return f"pace-w{weight:g}".replace(".", "p")


def progress(run: dict, kind: str, pure: dict) -> float:
    start = run["fit"][f"start_{kind}_rms"]
    final = run["fit"][f"refined_{kind}_rms"]
    pure_final = pure["fit"][f"refined_{kind}_rms"]
    return (start - final) / (start - pure_final)


first_failing_has_failure = False
for case_name, runs in score["cases"].items():
    pure = pure_score["cases"][case_name]["pure-rmc"]
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
        if value < guards["minimum_fit_progress_fraction"]:
            failures.append(
                f"{case_name}/{label(guards['highest_passing_weight'])}: "
                f"{kind} progress {value} too small"
            )
    failing = runs[label(guards["first_failing_weight"])]
    if any(
        progress(failing, kind, pure) < guards["minimum_fit_progress_fraction"]
        for kind in ("neutron", "xray")
    ):
        first_failing_has_failure = True

if not first_failing_has_failure:
    failures.append("pace-w100: expected a progress failure in at least one direction")
if failures:
    raise SystemExit("PACE weight-pilot verification failed:\n- " + "\n- ".join(failures))
print("passed bounded Erhard-2024 PACE weight-calibration checks")
