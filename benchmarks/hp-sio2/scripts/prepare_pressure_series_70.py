#!/usr/bin/env python3
"""Prepare the frozen 0-to-70 GPa series in 10 GPa increments."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tomllib
from pathlib import Path

import numpy as np

from pressure_common import (
    PAIR_TYPES,
    map_to_box,
    read_lammps,
    reciprocal_curves,
    structure_metrics,
    write_curve,
    write_lammps,
    write_rmc6f,
)


parser = argparse.ArgumentParser()
parser.add_argument("--force", action="store_true")
parser.add_argument("--pace-model", type=Path)
args = parser.parse_args()
case_root = Path(__file__).resolve().parents[1]
protocol_path = case_root / "expected/pressure-series-70.toml"
protocol = tomllib.loads(protocol_path.read_text())
source_root = case_root / "reference/local/pressure-series-70"
model = args.pace_model or (
    case_root.parent
    / "epsr-silica/reference/local/public-ace2024/SiOx_potential.yace"
)
if not (source_root / "IMPORT.json").is_file() or not model.is_file():
    raise SystemExit("import the pressure series and public ACE model first")
if hashlib.sha256(model.read_bytes()).hexdigest() != protocol["ace_model_sha256"]:
    raise SystemExit("PACE model hash differs from the frozen protocol")
root = case_root / "results/pressure-series-70"
if root.exists():
    if not args.force:
        raise SystemExit(f"output exists: {root} (pass --force)")
    shutil.rmtree(root)
root.mkdir(parents=True)


def sha256(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_config(path, structure, neutron, xray, method):
    fixture, design = protocol["fixture"], protocol["design"]
    nq = int(
        round(
            (fixture["q_max_a_inverse"] - fixture["q_min_a_inverse"])
            / fixture["q_step_a_inverse"]
        )
    )
    nbins = int(round(fixture["rdf_cutoff_a"] / fixture["rdf_bin_width_a"]))
    text = f'''[system]
structure = "{structure}"
format = "lammps"
types = {{ "1" = "Si", "2" = "O" }}

[data]
[data.neutron_sq]
file = "{neutron}"
sigma_column = 3
convention = "iq"
fit_min = {fixture["q_min_a_inverse"]}
fit_max = {fixture["q_max_a_inverse"]}
scattering_lengths = {{ Si = 4.1491, O = 5.803 }}
[data.xray_sq]
file = "{xray}"
sigma_column = 3
convention = "iq"
fit_min = {fixture["q_min_a_inverse"]}
fit_max = {fixture["q_max_a_inverse"]}

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
nq = {nq}
lorch = false
rdf_cutoff = {fixture["rdf_cutoff_a"]}
rdf_nbins = {nbins}
'''
    if method == "rsmith-pace-w30":
        text += f'''
[ml_potential]
backend = "pace_native"
model = "{model}"
weight = {design["pace_weight"]}
'''
    path.write_text(text)


fixture = protocol["fixture"]
q_values = np.arange(
    fixture["q_min_a_inverse"],
    fixture["q_max_a_inverse"],
    fixture["q_step_a_inverse"],
)
manifest = {
    "model_sha256": sha256(model),
    "protocol_sha256": sha256(protocol_path),
    "status": "hp_sio2_pressure_series_70_prepared",
    "steps": {},
}
for step in protocol["steps"]:
    name = step["name"]
    step_root, inputs = root / name, root / name / "inputs"
    inputs.mkdir(parents=True)
    start_positions, start_types, start_box, start_origin = read_lammps(
        source_root / f"{int(step['start_gpa']):03d}gpa.data"
    )
    target_positions, target_types, target_box, target_origin = read_lammps(
        source_root / f"{int(step['target_gpa']):03d}gpa.data"
    )
    if not np.array_equal(start_types, target_types):
        raise SystemExit(f"atom types differ in {name}")
    hidden_positions = map_to_box(
        target_positions, target_box, target_origin, target_box
    )
    start_mapped = map_to_box(start_positions, start_box, start_origin, target_box)
    hidden_path, start_path = inputs / "hidden-target.data", inputs / "cross-start.data"
    hidden_rmc, start_rmc = inputs / "hidden-target.rmc6f", inputs / "cross-start.rmc6f"
    write_lammps(hidden_path, f"hidden {step['target_gpa']} GPa target", hidden_positions, target_types, target_box)
    write_lammps(start_path, f"{step['start_gpa']} GPa start mapped to target box", start_mapped, start_types, target_box)
    write_rmc6f(hidden_rmc, f"hidden {step['target_gpa']} GPa target", hidden_positions, target_types, target_box)
    write_rmc6f(start_rmc, f"{step['start_gpa']} GPa start mapped to target box", start_mapped, start_types, target_box)
    target_metrics = structure_metrics(hidden_path, fixture)
    partials, neutron, xray = reciprocal_curves(
        target_metrics, q_values, fixture["rdf_bin_width_a"]
    )
    neutron_path, xray_path = inputs / "target-neutron-iq.dat", inputs / "target-xray-iq.dat"
    write_curve(neutron_path, q_values, neutron, fixture["sigma_iq"])
    write_curve(xray_path, q_values, xray, fixture["sigma_iq"])
    (inputs / "hidden-partial-sq.dat").write_text(
        "# Q S_SiSi S_SiO S_OO\n"
        + "".join(
            f"{q:.12g} "
            + " ".join(f"{partials[pair][index]:.12g}" for pair in PAIR_TYPES)
            + "\n"
            for index, q in enumerate(q_values)
        )
    )
    for method in ("rsmith-rmc", "rsmith-pace-w30"):
        run = step_root / "runs" / method
        run.mkdir(parents=True)
        write_config(run / "run.toml", start_path, neutron_path, xray_path, method)
    files = [hidden_path, start_path, hidden_rmc, start_rmc, neutron_path, xray_path, inputs / "hidden-partial-sq.dat"]
    manifest["steps"][name] = {
        "files": {path.name: sha256(path) for path in files},
        "start_box_a": start_box.tolist(),
        "target_box_a": target_box.tolist(),
    }
(root / "input-manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
print(json.dumps(manifest, indent=2, sort_keys=True))
