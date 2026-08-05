#!/usr/bin/env python3
"""Independently score the incremental high-pressure recovery steps."""

from __future__ import annotations

import argparse
import json
import math
import re
import tomllib
from pathlib import Path

import numpy as np

from pressure_common import PAIR_LABELS, reciprocal_curves, structure_metrics


parser = argparse.ArgumentParser()
parser.add_argument("--quiet", action="store_true")
args = parser.parse_args()
case_root = Path(__file__).resolve().parents[1]
root = case_root / "results/pressure-steps"
protocol = tomllib.loads((case_root / "expected/pressure-steps.toml").read_text())
fixture = protocol["fixture"]


def read_curve(path: Path):
    return np.asarray(
        [
            [float(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
    )


def coordination_tv(target, model):
    keys = set(target) | set(model)
    target_total, model_total = sum(target.values()), sum(model.values())
    return 0.5 * sum(
        abs(target.get(key, 0) / target_total - model.get(key, 0) / model_total)
        for key in keys
    )


def score_model(model, target, curves):
    rdf_rms = {
        pair: float(np.sqrt(np.mean((model["rdf"][pair] - target["rdf"][pair]) ** 2)))
        for pair in target["rdf"]
    }
    q_values = curves["neutron"][:, 0]
    _, neutron, xray = reciprocal_curves(
        model, q_values, fixture["rdf_bin_width_a"]
    )
    neutron_rms = float(np.sqrt(np.mean((neutron - curves["neutron"][:, 1]) ** 2)))
    xray_rms = float(np.sqrt(np.mean((xray - curves["xray"][:, 1]) ** 2)))
    quantile_errors, counts = {}, {}
    for label in PAIR_LABELS.values():
        quantile_errors[label], counts[label] = {}, {}
        target_values = target["_lower_tail_values"][label]
        model_values = model["_lower_tail_values"][label]
        for level in fixture["lower_tail_quantiles"]:
            key = f"{float(level):.3f}"
            threshold = float(target["lower_tail"][label]["quantiles_a"][key])
            quantile_errors[label][key] = abs(
                float(model["lower_tail"][label]["quantiles_a"][key]) - threshold
            )
            counts[label][key] = {
                "model": int(np.sum(model_values < threshold)),
                "model_fraction": float(np.mean(model_values < threshold)),
                "target": int(np.sum(target_values < threshold)),
                "target_fraction": float(np.mean(target_values < threshold)),
                "threshold_a": threshold,
            }
    all_errors = [value for pair in quantile_errors.values() for value in pair.values()]
    return {
        "counts_below_target_lower_tail": counts,
        "joint_iq_rms": math.sqrt((neutron_rms**2 + xray_rms**2) / 2.0),
        "lower_tail_quantile_error_a": quantile_errors,
        "mean_lower_tail_quantile_error_a": sum(all_errors) / len(all_errors),
        "mean_partial_rdf_rms": sum(rdf_rms.values()) / len(rdf_rms),
        "minimum_distance_a": {
            label: model["lower_tail"][label]["minimum_a"]
            for label in PAIR_LABELS.values()
        },
        "neutron_iq_rms": neutron_rms,
        "partial_rdf_rms": rdf_rms,
        "si_coordination_total_variation": coordination_tv(
            target["si_coordination"], model["si_coordination"]
        ),
        "xray_iq_rms": xray_rms,
    }


def run_metadata(run: Path):
    result = {
        "wall_seconds": float((run / "wall-seconds.txt").read_text())
        if (run / "wall-seconds.txt").is_file()
        else None
    }
    log = (run / "rsmith.log").read_text() if (run / "rsmith.log").is_file() else ""
    accepted = re.search(r"accepted (\d+)/(\d+)", log)
    if accepted:
        result["accepted_moves"] = int(accepted.group(1))
        result["attempted_moves"] = int(accepted.group(2))
    initial = re.search(r"Initial PACE/native energy =\s*([+\-0-9.eE]+)", log)
    final = re.findall(r"\[E:\s*([+\-0-9.eE]+)\]", log)
    if initial and final:
        result["energy_change_ev_per_atom"] = (
            float(final[-1]) - float(initial.group(1))
        ) / fixture["atoms"]
    return result


scores = {
    "status": "hp_sio2_incremental_pressure_steps_scored",
    "steps": {},
}
for step in protocol["steps"]:
    name = step["name"]
    step_root = root / name
    curves = {
        kind: read_curve(step_root / f"inputs/target-{kind}-iq.dat")
        for kind in ("neutron", "xray")
    }
    target = structure_metrics(step_root / "inputs/hidden-target.data", fixture)
    start = score_model(
        structure_metrics(step_root / "inputs/cross-start.data", fixture),
        target,
        curves,
    )
    result = {"start": start, "methods": {}, "native_forward": {}}
    for method in protocol["design"]["methods"]:
        run = step_root / "runs" / method
        structure = run / ("Cross.ato" if method == "native-epsr26" else "refined.xyz")
        if not structure.is_file():
            continue
        score = score_model(structure_metrics(structure, fixture), target, curves)
        score.update(run_metadata(run))
        result["methods"][method] = score
    forward = step_root / "native-forward-summary.json"
    if forward.is_file():
        result["native_forward"] = json.loads(forward.read_text())
    result["informative_step"] = (
        start["mean_partial_rdf_rms"] >= 0.02
        or start["si_coordination_total_variation"] >= 0.02
    )
    scores["steps"][name] = result

complete = all(
    set(step["methods"]) == set(protocol["design"]["methods"])
    for step in scores["steps"].values()
)
if complete:
    pace_better = all(
        step["methods"]["rsmith-pace-w30"]["mean_partial_rdf_rms"]
        < step["methods"]["native-epsr26"]["mean_partial_rdf_rms"]
        for step in scores["steps"].values()
    )
    lower_tail_guard = all(
        step["methods"]["rsmith-pace-w30"]["mean_lower_tail_quantile_error_a"]
        <= step["methods"]["native-epsr26"]["mean_lower_tail_quantile_error_a"]
        for step in scores["steps"].values()
    )
    scores["decision"] = {
        "all_steps_informative": all(
            step["informative_step"] for step in scores["steps"].values()
        ),
        "lower_tail_guard": lower_tail_guard,
        "pace_recovery_advantage": pace_better,
        "wait_for_20gpa_before_replication": True,
    }
(root / "score-summary.json").write_text(
    json.dumps(scores, indent=2, sort_keys=True) + "\n"
)
if not args.quiet:
    print(json.dumps(scores, indent=2, sort_keys=True))
