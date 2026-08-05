#!/usr/bin/env python3
"""Prepare the frozen 20-to-30 GPa S(Q)/g(r) domain comparison."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import tomllib
from pathlib import Path

import numpy as np

from pressure_common import (
    map_to_box,
    read_lammps,
    reciprocal_curves,
    structure_metrics,
    total_gr_from_iq,
    write_curve,
    write_lammps,
)


parser = argparse.ArgumentParser()
parser.add_argument("--force", action="store_true")
parser.add_argument("--pace-model", type=Path)
args = parser.parse_args()
case_root = Path(__file__).resolve().parents[1]
protocol_path = case_root / "expected/sq-gr-domain.toml"
protocol = tomllib.loads(protocol_path.read_text())
fixture, design = protocol["fixture"], protocol["design"]
source_root = case_root / "reference/local/pressure-series-70"
model = args.pace_model or (
    case_root.parent
    / "epsr-silica/reference/local/public-ace2024/SiOx_potential.yace"
)
start_source = source_root / "020gpa.data"
target_source = source_root / "030gpa.data"
for path, expected in (
    (start_source, protocol["start_sha256"]),
    (target_source, protocol["target_sha256"]),
    (model, protocol["ace_model_sha256"]),
):
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise SystemExit(f"missing or changed frozen input: {path}")
root = case_root / "results/sq-gr-domain"
if root.exists():
    if not args.force:
        raise SystemExit(f"output exists: {root} (pass --force)")
    shutil.rmtree(root)
inputs = root / "inputs"
inputs.mkdir(parents=True)

start_positions, start_types, start_box, start_origin = read_lammps(start_source)
target_positions, target_types, target_box, target_origin = read_lammps(target_source)
if not np.array_equal(start_types, target_types):
    raise SystemExit("20 and 30 GPa atom types differ")
start_mapped = map_to_box(start_positions, start_box, start_origin, target_box)
hidden_positions = map_to_box(target_positions, target_box, target_origin, target_box)
start_path, target_path = inputs / "cross-start.data", inputs / "hidden-target.data"
write_lammps(start_path, "20 GPa start mapped to 30 GPa box", start_mapped, start_types, target_box)
write_lammps(target_path, "hidden 30 GPa target", hidden_positions, target_types, target_box)

q_values = np.arange(
    fixture["q_min_a_inverse"],
    fixture["q_max_a_inverse"],
    fixture["q_step_a_inverse"],
)
r_values = (
    np.arange(round(fixture["rdf_cutoff_a"] / fixture["rdf_bin_width_a"])) + 0.5
) * fixture["rdf_bin_width_a"]
start_metrics = structure_metrics(start_path, fixture)
target_metrics = structure_metrics(target_path, fixture)
_, _, start_xray = reciprocal_curves(
    start_metrics, q_values, fixture["rdf_bin_width_a"]
)
_, _, target_xray = reciprocal_curves(
    target_metrics, q_values, fixture["rdf_bin_width_a"]
)
density = fixture["atoms"] / math.prod(target_metrics["box_a"])
start_gr = total_gr_from_iq(
    q_values,
    start_xray,
    r_values,
    density,
    fixture["gr_qmax_a_inverse"],
    fixture["gr_lorch"],
)
target_gr = total_gr_from_iq(
    q_values,
    target_xray,
    r_values,
    density,
    fixture["gr_qmax_a_inverse"],
    fixture["gr_lorch"],
)
r_fit = (r_values >= fixture["gr_fit_min_a"]) & (
    r_values <= fixture["gr_fit_max_a"]
)
sq_start_chi2 = float(
    np.sum(((start_xray - target_xray) / fixture["sq_sigma_iq"]) ** 2)
)
gr_sigma = math.sqrt(float(np.sum((start_gr[r_fit] - target_gr[r_fit]) ** 2)) / sq_start_chi2)
sq_path, gr_path = inputs / "target-xray-iq.dat", inputs / "target-xray-gr.dat"
write_curve(sq_path, q_values, target_xray, fixture["sq_sigma_iq"])
gr_path.write_text(
    "# r g_X(r) sigma\n"
    + "".join(
        f"{r:.12g} {value:.12g} {gr_sigma:.12g}\n"
        for r, value in zip(r_values, target_gr)
    )
)


def config_text(method):
    weights = {
        "sq-only": design["sq_only_weights"],
        "gr-only": design["gr_only_weights"],
        "sq-plus-gr": design["joint_weights"],
    }[method]
    text = f'''[system]
structure = "{start_path}"
format = "lammps"
types = {{ "1" = "Si", "2" = "O" }}

[data]
'''
    if weights["sq"]:
        text += f'''[data.xray_sq]
file = "{sq_path}"
sigma_column = 3
convention = "iq"
weight = {weights["sq"]}
fit_min = {fixture["q_min_a_inverse"]}
fit_max = {fixture["q_max_a_inverse"]}
'''
    if weights["gr"]:
        text += f'''[data.xray_gr]
file = "{gr_path}"
sigma_column = 3
weight = {weights["gr"]}
fit_min = {fixture["gr_fit_min_a"]}
fit_max = {fixture["gr_fit_max_a"]}
qmax = {fixture["gr_qmax_a_inverse"]}
lorch = {str(fixture["gr_lorch"]).lower()}
qdamp = {fixture["gr_qdamp_a_inverse"]}
'''
    text += f'''
[rmc]
max_moves = {design["moves"]}
max_step = 0.1
checkpoint_every = 1000000000
seed = {design["seed"]}
print_every = {design["moves"]}
target_acceptance = 0.25
adjust_step_every = 1000

[sq]
qmin = {fixture["q_min_a_inverse"]}
qmax = {fixture["q_max_a_inverse"]}
nq = {len(q_values)}
lorch = false
rdf_cutoff = {fixture["rdf_cutoff_a"]}
rdf_nbins = {len(r_values)}

[ml_potential]
backend = "pace_native"
model = "{model}"
weight = {design["pace_weight"]}
'''
    return text


for method in design["methods"]:
    run = root / "runs" / method
    run.mkdir(parents=True)
    (run / "run.toml").write_text(config_text(method))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


manifest = {
    "files": {
        path.name: sha256(path) for path in (start_path, target_path, sq_path, gr_path)
    },
    "gr_sigma": gr_sigma,
    "gr_start_chi2": float(np.sum(((start_gr[r_fit] - target_gr[r_fit]) / gr_sigma) ** 2)),
    "model_sha256": sha256(model),
    "protocol_sha256": sha256(protocol_path),
    "sq_start_chi2": sq_start_chi2,
    "start_box_a": start_box.tolist(),
    "status": "hp_sio2_sq_gr_domain_prepared",
    "target_box_a": target_box.tolist(),
}
(root / "input-manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
print(json.dumps(manifest, indent=2, sort_keys=True))
