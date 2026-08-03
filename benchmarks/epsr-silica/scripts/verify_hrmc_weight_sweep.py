#!/usr/bin/env python3
"""Verify the bounded HRMC weight pilot against its committed guards."""

from __future__ import annotations

import json
import math
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
expected = tomllib.loads(
    (case_root / "expected/hrmc-weight-sweep-smoke.toml").read_text()
)
score = json.loads(
    (case_root / "results/hrmc-weight-sweep/score-summary.json").read_text()
)
guards = expected["guards"]
failures: list[str] = []


def label(model: str, weight: float) -> str:
    return f"{model}-w{weight:g}".replace(".", "p")


def progress(run: dict, kind: str, pure: dict) -> float:
    start = run["fit"][f"start_{kind}_rms"]
    final = run["fit"][f"refined_{kind}_rms"]
    pure_final = pure["fit"][f"refined_{kind}_rms"]
    return (start - final) / (start - pure_final)


for case_name, runs in score["cases"].items():
    pure = runs.get("pure-rmc")
    if pure is None:
        failures.append(f"{case_name}: missing pure-RMC control")
        continue
    for model in ("pedone", "gap"):
        weights = guards[f"required_{model}_weights"]
        for weight in weights:
            run_name = label(model, weight)
            if run_name not in runs:
                failures.append(f"{case_name}: missing {run_name}")
                continue
            run = runs[run_name]
            fit = run["fit"]
            if fit.get("attempted_moves") != guards["attempted_moves"]:
                failures.append(f"{case_name}/{run_name}: wrong move count")
            wall = fit.get("wall_seconds", 0.0)
            if not math.isfinite(wall) or wall <= 0.0:
                failures.append(f"{case_name}/{run_name}: invalid wall time")
            minima = run["model_metrics"]["minimum_distance_a"]
            for pair, threshold in guards["minimum_distance_a"].items():
                if minima[pair] + 1.0e-9 < threshold:
                    failures.append(
                        f"{case_name}/{run_name}: {pair} minimum {minima[pair]} < {threshold}"
                    )

    for weight in expected["selection"]["pedone_production_bracket"]:
        run_name = label("pedone", weight)
        for kind in ("neutron", "xray"):
            if (
                progress(runs[run_name], kind, pure)
                < guards["minimum_fit_progress_fraction"]
            ):
                failures.append(f"{case_name}/{run_name}: {kind} progress guard failed")

    failing_pedone = label("pedone", guards["pedone_first_guard_failing_weight"])
    if all(
        progress(runs[failing_pedone], kind, pure)
        >= guards["minimum_fit_progress_fraction"]
        for kind in ("neutron", "xray")
    ):
        failures.append(
            f"{case_name}/{failing_pedone}: expected a progress guard failure"
        )

    gap_passing = label("gap", guards["gap_highest_guard_passing_weight"])
    for kind in ("neutron", "xray"):
        if (
            progress(runs[gap_passing], kind, pure)
            < guards["gap_weight_0p3_minimum_fit_progress_fraction"]
        ):
            failures.append(f"{case_name}/{gap_passing}: {kind} progress too small")
    gap_failing = label("gap", guards["gap_first_guard_failing_weight"])
    if (
        runs[gap_failing]["fit"]["accepted_moves"]
        > guards["gap_weight_1_maximum_accepted_moves"]
    ):
        failures.append(f"{case_name}/{gap_failing}: acceptance cliff not reproduced")

if failures:
    raise SystemExit(
        "HRMC weight-pilot verification failed:\n- " + "\n- ".join(failures)
    )
print("passed bounded Pedone/GAP HRMC energy-weight pilot checks")
