#!/usr/bin/env python3
"""Independently score all four methods over the 0-to-70 GPa series."""

from __future__ import annotations

import argparse
import csv
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
root = case_root / "results/pressure-series-70"
protocol = tomllib.loads((case_root / "expected/pressure-series-70.toml").read_text())
fixture = protocol["fixture"]


def read_curve(path: Path):
    return np.asarray(
        [
            [float(value) for value in line.split()[:2]]
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
    _, neutron, xray = reciprocal_curves(model, q_values, fixture["rdf_bin_width_a"])
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


def run_metadata(run: Path, method: str):
    result = {
        "wall_seconds": float((run / "wall-seconds.txt").read_text())
        if (run / "wall-seconds.txt").is_file()
        else None
    }
    if method == "native-rmcprofile":
        audit = run / "move-audit.json"
        if audit.is_file():
            data = json.loads(audit.read_text())
            result["accepted_moves"] = data["total_accepted"]
            result["attempted_moves"] = data["total_moves"]
        return result
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


structures = {
    "rsmith-rmc": "refined.xyz",
    "rsmith-pace-w30": "refined.xyz",
    "native-epsr26": "Cross.ato",
    "native-rmcprofile": "cross.rmc6f",
}
scores = {"status": "hp_sio2_pressure_series_70_scored", "steps": {}}
csv_rows = []
for step in protocol["steps"]:
    name = step["name"]
    step_root = root / name
    curves = {
        kind: read_curve(step_root / f"inputs/target-{kind}-iq.dat")
        for kind in ("neutron", "xray")
    }
    target = structure_metrics(step_root / "inputs/hidden-target.data", fixture)
    start = score_model(
        structure_metrics(step_root / "inputs/cross-start.data", fixture), target, curves
    )
    result = {"start": start, "methods": {}, "native_forward": {}}
    for method in protocol["design"]["methods"]:
        run = step_root / "runs" / method
        structure = run / structures[method]
        if not structure.is_file():
            continue
        score = score_model(structure_metrics(structure, fixture), target, curves)
        score.update(run_metadata(run, method))
        result["methods"][method] = score
        csv_rows.append(
            {
                "step": name,
                "target_gpa": step["target_gpa"],
                "method": method,
                "joint_iq_rms": score["joint_iq_rms"],
                "mean_partial_rdf_rms": score["mean_partial_rdf_rms"],
                "si_coordination_tv": score["si_coordination_total_variation"],
                "mean_lower_tail_error_a": score["mean_lower_tail_quantile_error_a"],
                "wall_seconds": score["wall_seconds"],
            }
        )
    for label, filename in (
        ("native-epsr26", "native-forward-summary.json"),
        ("native-rmcprofile", "native-rmcprofile-forward-summary.json"),
    ):
        path = step_root / filename
        if path.is_file():
            result["native_forward"][label] = json.loads(path.read_text())
    result["informative_step"] = (
        start["mean_partial_rdf_rms"] >= 0.02
        or start["si_coordination_total_variation"] >= 0.02
    )
    if result["methods"]:
        result["lowest_hidden_rdf_method"] = min(
            result["methods"],
            key=lambda method: result["methods"][method]["mean_partial_rdf_rms"],
        )
    scores["steps"][name] = result

complete = all(
    set(step["methods"]) == set(protocol["design"]["methods"])
    for step in scores["steps"].values()
)
scores["complete"] = complete
if complete:
    scores["series"] = {
        "all_steps_informative": all(
            step["informative_step"] for step in scores["steps"].values()
        ),
        "lowest_hidden_rdf_counts": {
            method: sum(
                step["lowest_hidden_rdf_method"] == method
                for step in scores["steps"].values()
            )
            for method in protocol["design"]["methods"]
        },
    }
(root / "score-summary.json").write_text(
    json.dumps(scores, indent=2, sort_keys=True) + "\n"
)
if csv_rows:
    with (root / "plot-data.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=csv_rows[0])
        writer.writeheader()
        writer.writerows(csv_rows)
if not args.quiet:
    print(json.dumps(scores, indent=2, sort_keys=True))
