#!/usr/bin/env python3
"""Score the frozen reciprocal/local-real/joint PACE comparison."""

from __future__ import annotations

import json
import math
import re
import tomllib
from pathlib import Path

import numpy as np

from pressure_common import reciprocal_curves, structure_metrics, total_gr_from_iq


case_root = Path(__file__).resolve().parents[1]
root = case_root / "results/sq-gr-domain"
protocol = tomllib.loads((case_root / "expected/sq-gr-domain.toml").read_text())
fixture = protocol["fixture"]


def read_curve(path):
    return np.asarray(
        [
            [float(value) for value in line.split()[:2]]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
    )


def coordination_tv(target, model):
    keys = set(target) | set(model)
    nt, nm = sum(target.values()), sum(model.values())
    return 0.5 * sum(
        abs(target.get(key, 0) / nt - model.get(key, 0) / nm) for key in keys
    )


q_target = read_curve(root / "inputs/target-xray-iq.dat")
gr_target = read_curve(root / "inputs/target-xray-gr.dat")
q_values, r_values = q_target[:, 0], gr_target[:, 0]
r_local = (r_values >= fixture["gr_fit_min_a"]) & (
    r_values <= fixture["gr_fit_max_a"]
)
rdf_radius = (
    np.arange(round(fixture["rdf_cutoff_a"] / fixture["rdf_bin_width_a"])) + 0.5
) * fixture["rdf_bin_width_a"]
rdf_local = (rdf_radius >= fixture["gr_fit_min_a"]) & (
    rdf_radius <= fixture["gr_fit_max_a"]
)
target = structure_metrics(root / "inputs/hidden-target.data", fixture)


def score(path):
    model = structure_metrics(path, fixture)
    rdf_full, rdf_local_values = [], []
    for pair in target["rdf"]:
        delta = model["rdf"][pair] - target["rdf"][pair]
        rdf_full.append(float(np.sqrt(np.mean(delta**2))))
        rdf_local_values.append(float(np.sqrt(np.mean(delta[rdf_local] ** 2))))
    _, _, xray = reciprocal_curves(model, q_values, fixture["rdf_bin_width_a"])
    density = fixture["atoms"] / math.prod(model["box_a"])
    gr = total_gr_from_iq(
        q_values,
        xray,
        r_values,
        density,
        fixture["gr_qmax_a_inverse"],
        fixture["gr_lorch"],
    )
    tail_errors = []
    for pair in target["lower_tail"]:
        for key, value in target["lower_tail"][pair]["quantiles_a"].items():
            tail_errors.append(abs(model["lower_tail"][pair]["quantiles_a"][key] - value))
    return {
        "joint_hidden_partial_rdf_rms": float(np.mean(rdf_full)),
        "local_hidden_partial_rdf_rms": float(np.mean(rdf_local_values)),
        "mean_lower_tail_quantile_error_a": float(np.mean(tail_errors)),
        "minimum_distance_a": {
            pair: model["lower_tail"][pair]["minimum_a"]
            for pair in target["lower_tail"]
        },
        "si_coordination_total_variation": coordination_tv(
            target["si_coordination"], model["si_coordination"]
        ),
        "xray_gr_local_rms": float(
            np.sqrt(np.mean((gr[r_local] - gr_target[r_local, 1]) ** 2))
        ),
        "xray_iq_rms": float(np.sqrt(np.mean((xray - q_target[:, 1]) ** 2))),
    }


result = {
    "status": "hp_sio2_sq_gr_domain_scored",
    "start": score(root / "inputs/cross-start.data"),
    "methods": {},
}
for method in protocol["design"]["methods"]:
    run = root / "runs" / method
    if not (run / "refined.xyz").is_file():
        continue
    values = score(run / "refined.xyz")
    values["wall_seconds"] = float((run / "wall-seconds.txt").read_text())
    log = (run / "rsmith.log").read_text()
    accepted = re.search(r"accepted (\d+)/(\d+)", log)
    if accepted:
        values["accepted_moves"] = int(accepted.group(1))
        values["attempted_moves"] = int(accepted.group(2))
    initial = re.search(r"Initial PACE/native energy =\s*([+\-0-9.eE]+)", log)
    final = re.findall(r"\[E:\s*([+\-0-9.eE]+)\]", log)
    if initial and final:
        values["energy_change_ev_per_atom"] = (
            float(final[-1]) - float(initial.group(1))
        ) / fixture["atoms"]
    result["methods"][method] = values
result["complete"] = set(result["methods"]) == set(protocol["design"]["methods"])
(root / "score-summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(json.dumps(result, indent=2, sort_keys=True))
