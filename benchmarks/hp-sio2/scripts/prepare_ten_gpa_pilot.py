#!/usr/bin/env python3
"""Prepare the frozen 0-to-10 GPa synthetic SiO2 recovery pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import tomllib
from collections import Counter
from pathlib import Path


TYPE_NAMES = {1: "Si", 2: "O"}
PAIR_TYPES = ((1, 1), (1, 2), (2, 2))
PAIR_LABELS = {(1, 1): "Si-Si", (1, 2): "Si-O", (2, 2): "O-O"}
NEUTRON_B = {1: 4.1491, 2: 5.803}
XRAY = {
    1: (
        [6.2915, 3.0353, 1.9891, 1.541],
        [2.4386, 32.3337, 0.6785, 81.6937],
        1.1407,
    ),
    2: (
        [3.0485, 2.2868, 1.5463, 0.867],
        [13.2771, 5.7011, 0.3239, 32.9089],
        0.2508,
    ),
}

parser = argparse.ArgumentParser()
parser.add_argument("--force", action="store_true")
parser.add_argument("--pace-model", type=Path)
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
repo_root = case_root.parents[1]
protocol_path = case_root / "expected/ten-gpa-pilot.toml"
protocol = tomllib.loads(protocol_path.read_text())
source_root = case_root / "reference/local/ten-gpa"
model = args.pace_model or (
    repo_root
    / "benchmarks/epsr-silica/reference/local/public-ace2024/SiOx_potential.yace"
)
for required in (
    source_root / "IMPORT.json",
    source_root / "ambient-0gpa.data",
    source_root / "target-10gpa.data",
    model,
):
    if not required.is_file():
        raise SystemExit(f"missing prerequisite: {required}")
if hashlib.sha256(model.read_bytes()).hexdigest() != protocol["ace"]["model_sha256"]:
    raise SystemExit("PACE model does not match frozen Erhard-2024 hash")

result_root = case_root / "results/ten-gpa-pilot"
if result_root.exists():
    if not args.force:
        raise SystemExit(f"output exists: {result_root} (pass --force)")
    shutil.rmtree(result_root)
input_root = result_root / "inputs"
run_root = result_root / "runs"
input_root.mkdir(parents=True)
run_root.mkdir()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_lammps(path: Path):
    lines = path.read_text().splitlines()
    bounds = {}
    atoms = []
    atom_style = None
    in_atoms = False
    for line in lines:
        fields = line.split()
        if len(fields) >= 4 and fields[-2:] in (
            ["xlo", "xhi"],
            ["ylo", "yhi"],
            ["zlo", "zhi"],
        ):
            bounds[fields[-2][0]] = (float(fields[0]), float(fields[1]))
            continue
        if line.strip().startswith("Atoms"):
            atom_style = "charge" if "charge" in line else "atomic"
            in_atoms = True
            continue
        if not in_atoms or not line.strip():
            continue
        if line.lstrip()[0].isalpha():
            in_atoms = False
            continue
        fields = line.split()
        offset = 3 if atom_style == "charge" else 2
        atoms.append(
            (
                int(fields[0]),
                int(fields[1]),
                tuple(float(value) for value in fields[offset : offset + 3]),
            )
        )
    atoms.sort()
    box = tuple(bounds[axis][1] - bounds[axis][0] for axis in "xyz")
    origin = tuple(bounds[axis][0] for axis in "xyz")
    if len(atoms) != protocol["source"]["atoms"] or set(bounds) != set("xyz"):
        raise ValueError(f"unexpected structure content in {path}")
    return atoms, box, origin


def map_to_box(atoms, old_box, old_origin, new_box):
    return [
        (
            atom_id,
            atom_type,
            tuple(
                (((position[axis] - old_origin[axis]) / old_box[axis]) % 1.0)
                * new_box[axis]
                for axis in range(3)
            ),
        )
        for atom_id, atom_type, position in atoms
    ]


def write_lammps(path: Path, title: str, atoms, box):
    path.write_text(
        f"{title}\n\n{len(atoms)} atoms\n2 atom types\n\n"
        + "".join(
            f"0.0 {length:.15g} {axis}lo {axis}hi\n"
            for axis, length in zip("xyz", box)
        )
        + "\nMasses\n\n1 28.0855\n2 15.9994\n\nAtoms # charge\n\n"
        + "".join(
            f"{atom_id} {atom_type} 0.0 "
            + " ".join(f"{value:.15g}" for value in position)
            + "\n"
            for atom_id, atom_type, position in atoms
        )
    )


def displacement(first, second, box):
    return tuple(
        first[axis]
        - second[axis]
        - box[axis] * round((first[axis] - second[axis]) / box[axis])
        for axis in range(3)
    )


def partial_rdf(atoms, box, dr, rmax):
    counts = Counter(atom_type for _, atom_type, _ in atoms)
    positions = [position for _, _, position in atoms]
    types = [atom_type for _, atom_type, _ in atoms]
    nbins = int(round(rmax / dr))
    histograms = {pair: [0] * nbins for pair in PAIR_TYPES}
    for first in range(len(atoms) - 1):
        for second in range(first + 1, len(atoms)):
            vector = displacement(positions[second], positions[first], box)
            distance = math.sqrt(sum(value * value for value in vector))
            if distance < rmax:
                pair = tuple(sorted((types[first], types[second])))
                histograms[pair][int(distance / dr)] += 1
    volume = math.prod(box)
    radius = [(index + 0.5) * dr for index in range(nbins)]
    result = {}
    for pair in PAIR_TYPES:
        first_type, second_type = pair
        values = []
        for index, count in enumerate(histograms[pair]):
            inner, outer = index * dr, (index + 1) * dr
            shell = 4.0 * math.pi * (outer**3 - inner**3) / 3.0
            denominator = counts[first_type] * (counts[second_type] / volume) * shell
            factor = 2.0 if first_type == second_type else 1.0
            values.append(factor * count / denominator)
        result[pair] = values
    return radius, result


def partial_sq(radius, rdf, number_density, q_values, dr):
    prefactor = 4.0 * math.pi * number_density * dr
    return {
        pair: [
            1.0
            + prefactor
            * sum(
                r * (g - 1.0) * math.sin(q * r)
                for r, g in zip(radius, values)
            )
            / q
            for q in q_values
        ]
        for pair, values in rdf.items()
    }


def xray_form_factor(atom_type, q):
    amplitudes, decays, constant = XRAY[atom_type]
    s2 = (q / (4.0 * math.pi)) ** 2
    return constant + sum(
        amplitude * math.exp(-decay * s2)
        for amplitude, decay in zip(amplitudes, decays)
    )


def totals(partials, q_values):
    concentrations = {1: 1.0 / 3.0, 2: 2.0 / 3.0}
    mean_b = sum(concentrations[k] * NEUTRON_B[k] for k in concentrations)
    neutron_weights = {
        pair: (1.0 if pair[0] == pair[1] else 2.0)
        * concentrations[pair[0]]
        * concentrations[pair[1]]
        * NEUTRON_B[pair[0]]
        * NEUTRON_B[pair[1]]
        / mean_b**2
        for pair in PAIR_TYPES
    }
    neutron = [
        sum(neutron_weights[pair] * partials[pair][index] for pair in PAIR_TYPES)
        for index in range(len(q_values))
    ]
    xray = []
    for index, q in enumerate(q_values):
        form = {atom_type: xray_form_factor(atom_type, q) for atom_type in (1, 2)}
        mean_f = sum(concentrations[k] * form[k] for k in concentrations)
        xray.append(
            sum(
                (1.0 if pair[0] == pair[1] else 2.0)
                * concentrations[pair[0]]
                * concentrations[pair[1]]
                * form[pair[0]]
                * form[pair[1]]
                * partials[pair][index]
                / mean_f**2
                for pair in PAIR_TYPES
            )
        )
    return neutron, xray


def write_curve(path, q_values, values, sigma):
    path.write_text(
        "# Q i(Q)=S(Q)-1 sigma\n"
        + "".join(
            f"{q:.12g} {value - 1.0:.12g} {sigma:.12g}\n"
            for q, value in zip(q_values, values)
        )
    )


def write_config(path: Path, structure: Path, neutron: Path, xray: Path, moves: int, seed: int, pace: bool):
    fixture = protocol["fixture"]
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
max_moves = {moves}
max_step = 0.1
checkpoint_every = 1000000000
seed = {seed}
print_every = {moves}
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
    if pace:
        text += f'''
[ml_potential]
backend = "pace_native"
model = "{model}"
weight = {protocol["pilot"]["pace_weight"]}
'''
    path.write_text(text)


ambient_atoms, ambient_box, ambient_origin = read_lammps(
    source_root / "ambient-0gpa.data"
)
target_atoms, target_box, target_origin = read_lammps(
    source_root / "target-10gpa.data"
)
if Counter(atom_type for _, atom_type, _ in target_atoms) != Counter(
    {1: protocol["source"]["silicon_atoms"], 2: protocol["source"]["oxygen_atoms"]}
):
    raise SystemExit("unexpected target composition")
hidden_atoms = map_to_box(target_atoms, target_box, target_origin, target_box)
start_atoms = map_to_box(ambient_atoms, ambient_box, ambient_origin, target_box)
hidden_path = input_root / "hidden-target.data"
start_path = input_root / "cross-start.data"
write_lammps(hidden_path, "hidden ACE 10 GPa target", hidden_atoms, target_box)
write_lammps(start_path, "0 GPa ACE glass affinely mapped to 10 GPa box", start_atoms, target_box)

fixture = protocol["fixture"]
q_values = [
    fixture["q_min_a_inverse"] + index * fixture["q_step_a_inverse"]
    for index in range(
        int(
            round(
                (fixture["q_max_a_inverse"] - fixture["q_min_a_inverse"])
                / fixture["q_step_a_inverse"]
            )
        )
    )
]
radius, rdf = partial_rdf(
    hidden_atoms,
    target_box,
    fixture["rdf_bin_width_a"],
    fixture["rdf_cutoff_a"],
)
partials = partial_sq(
    radius,
    rdf,
    len(hidden_atoms) / math.prod(target_box),
    q_values,
    fixture["rdf_bin_width_a"],
)
neutron, xray = totals(partials, q_values)
neutron_path = input_root / "target-neutron-iq.dat"
xray_path = input_root / "target-xray-iq.dat"
write_curve(neutron_path, q_values, neutron, fixture["sigma_iq"])
write_curve(xray_path, q_values, xray, fixture["sigma_iq"])
(input_root / "hidden-partial-sq.dat").write_text(
    "# Q S_SiSi S_SiO S_OO\n"
    + "".join(
        f"{q:.12g} "
        + " ".join(f"{partials[pair][index]:.12g}" for pair in PAIR_TYPES)
        + "\n"
        for index, q in enumerate(q_values)
    )
)

for method in protocol["pilot"]["methods"]:
    for moves in protocol["pilot"]["checkpoints_moves"]:
        run = run_root / method / f"moves-{int(moves):06d}"
        run.mkdir(parents=True)
        write_config(
            run / "run.toml",
            start_path,
            neutron_path,
            xray_path,
            int(moves),
            int(protocol["pilot"]["seed"]),
            method == "rsmith-pace-w3",
        )

summary = {
    "status": "hp_sio2_10gpa_pilot_prepared",
    "protocol_sha256": sha256(protocol_path),
    "pace_model": str(model),
    "pace_model_sha256": sha256(model),
    "ambient_source_sha256": sha256(source_root / "ambient-0gpa.data"),
    "target_source_sha256": sha256(source_root / "target-10gpa.data"),
    "ambient_box_a": ambient_box,
    "target_box_a": target_box,
    "affine_scale": [target_box[i] / ambient_box[i] for i in range(3)],
    "files": {
        path.name: sha256(path)
        for path in (
            hidden_path,
            start_path,
            neutron_path,
            xray_path,
            input_root / "hidden-partial-sq.dat",
        )
    },
}
(input_root / "fixture-summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(summary, indent=2, sort_keys=True))
