#!/usr/bin/env python3
"""Validate the DTBsilicaNX reference potential and run one joint EPSR epoch.

The licensed EPSR26 tutorial files remain local and ignored.  This script
parses their atom and potential-control records, reconstructs EPSR's reference
pair potentials independently, checks them against a fresh native ``.o01``,
and writes rsmith tabulated potentials in eV.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import time
from pathlib import Path


KJ_PER_MOL_PER_EV = 96.48533212
COULOMB_KJ_A_PER_MOL = 1390.0
CURVE_LIMIT = 2.0e-7
ENERGY_LIMIT = 0.01


def numeric_rows(path: Path, minimum_columns: int = 2):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        try:
            values = tuple(float(value) for value in line.split())
        except ValueError:
            continue
        if len(values) >= minimum_columns:
            rows.append(values)
    return rows


def atom_parameters(path: Path, expected_species: str):
    lines = path.read_text().splitlines()
    for index, line in enumerate(lines[:-1]):
        words = line.split()
        if len(words) >= 2 and words[0] == expected_species and words[1] == expected_species:
            values = tuple(float(value) for value in lines[index + 1].split())
            if len(values) < 5:
                break
            return {
                "epsilon_kj_mol": values[0],
                "sigma_a": values[1],
                "mass": values[2],
                "charge_e": values[3],
                "charge_radius_a": values[4],
            }
    raise ValueError(f"could not find {expected_species} pair parameters in {path}")


def pcof_setting(path: Path, name: str):
    pattern = re.compile(rf"^\s*{re.escape(name)}\s+([^\s]+)", re.MULTILINE)
    match = pattern.search(path.read_text())
    if match is None:
        raise ValueError(f"setting {name!r} not found in {path}")
    return float(match.group(1))


def cosine_cutoff(r: float, onset: float, cutoff: float):
    if r <= onset:
        return 1.0
    if r >= cutoff:
        return 0.0
    return 0.5 * (1.0 + math.cos(math.pi * (r - onset) / (cutoff - onset)))


def hummer_cutoff(r: float, cutoff: float):
    if r >= cutoff:
        return 0.0
    ratio = r / cutoff
    return (1.0 - ratio) ** 4 * (1.0 + ratio * (1.6 + 0.4 * ratio))


def reference_pair_value(r: float, pair: dict, controls: dict):
    # EPSR replaces its first r=0 table entry with the next finite-grid value.
    radius = max(r, 0.0075)
    power = controls["power"]
    if power <= 6.0:
        raise ValueError("EPSR modified Lennard-Jones power must exceed 6")
    ratio = 6.0 / power
    root = math.sqrt(ratio)
    prefactor = 1.0 / ((1.0 / root - root) * root ** ((power + 6.0) / (power - 6.0)))
    diameter_ratio = pair["sigma_a"] / radius
    short_range = prefactor * pair["epsilon_kj_mol"] * (
        diameter_ratio**power - diameter_ratio**6
    )
    short_range *= cosine_cutoff(radius, controls["rminpt"], controls["rmaxpt"])

    # DTBsilica has zero charge-cloud radii, so EPSR's chsph() reduces to 1/r.
    if pair["charge_radius_a"] != (0.0, 0.0):
        raise ValueError("non-zero EPSR charge-cloud radii are not implemented by this benchmark")
    coulomb = (
        COULOMB_KJ_A_PER_MOL
        * pair["charge_product"]
        / radius
        * hummer_cutoff(radius, controls["rmaxpt"])
    )
    return controls["refpotfac"] * (short_range + coulomb)


def pair_parameters(atoms: dict):
    pairs = []
    for first, second in (("Si", "Si"), ("Si", "O"), ("O", "O")):
        a, b = atoms[first], atoms[second]
        pairs.append(
            {
                "label": f"{first}-{second}",
                "epsilon_kj_mol": math.sqrt(a["epsilon_kj_mol"] * b["epsilon_kj_mol"]),
                "sigma_a": 0.5 * (a["sigma_a"] + b["sigma_a"]),
                "charge_product": a["charge_e"] * b["charge_e"],
                "charge_radius_a": (a["charge_radius_a"], b["charge_radius_a"]),
            }
        )
    return pairs


def curve_metrics(native_o01: Path, pairs: list[dict], controls: dict):
    rows = numeric_rows(native_o01, 7)
    metrics = {}
    for pair_index, pair in enumerate(pairs):
        reference = []
        residual = []
        for row in rows:
            radius = row[0]
            if not 0.5 <= radius < controls["rmaxpt"]:
                continue
            expected = row[1 + 2 * pair_index]
            calculated = reference_pair_value(radius, pair, controls)
            reference.append(expected)
            residual.append(calculated - expected)
        rms = math.sqrt(sum(value * value for value in residual) / len(residual))
        dynamic_range = max(reference) - min(reference)
        metrics[pair["label"]] = {
            "n_points": len(reference),
            "rms_kj_mol": rms,
            "max_absolute_kj_mol": max(abs(value) for value in residual),
            "rms_over_dynamic_range": rms / dynamic_range,
        }
    return metrics


def write_potential_tables(output: Path, pairs: list[dict], controls: dict):
    output.mkdir(parents=True, exist_ok=True)
    paths = {}
    for pair in pairs:
        path = output / f"epsr-reference-{pair['label']}.dat"
        with path.open("w") as stream:
            stream.write("# r(A) V_ref(eV); independently reconstructed from local EPSR26 inputs\n")
            for index in range(int(round(controls["rmaxpt"] / 0.001)) + 1):
                radius = index * 0.001
                value = reference_pair_value(radius, pair, controls) / KJ_PER_MOL_PER_EV
                stream.write(f"{radius:.3f} {value:.15g}\n")
        paths[pair["label"]] = path
    return paths


def xray_parameters(path: Path):
    rows = [row for row in numeric_rows(path, 11) if len(row) == 11]
    if len(rows) != 2:
        raise ValueError(f"expected two X-ray form-factor rows in {path}")
    return rows


def xray_form_factor(values, q: float):
    s2 = (q / (4.0 * math.pi)) ** 2
    return values[0] + sum(
        amplitude * math.exp(-decay * s2)
        for amplitude, decay in zip(values[1::2], values[2::2])
    )


def write_faber_ziman_targets(source: Path, native: Path, output: Path):
    # EPSR .u01 stores the calculated total scattering.  The similarly shaped
    # .v01 file stores data-minus-model residuals and must not be used as a
    # refinement target.
    rows = numeric_rows(native / "DTBsilica.EPSR.u01", 5)
    neutron_weights = [row[2] for row in numeric_rows(source / "DTBsilica.NWTStot.wts", 3) if len(row) == 3]
    if len(neutron_weights) != 3:
        raise ValueError("expected three neutron pair weights")
    neutron_scale = 1.0 / sum(neutron_weights)
    xray_rows = xray_parameters(source / "DTBsilica.XWTS.wts")

    neutron_path = output / "target-neutron-faber-ziman-iq.dat"
    xray_path = output / "target-xray-faber-ziman-iq.dat"
    with neutron_path.open("w") as neutron, xray_path.open("w") as xray:
        neutron.write("# Q i_FZ(Q) sigma_FZ\n")
        xray.write("# Q i_FZ(Q) sigma_FZ\n")
        for row in rows:
            q = row[0]
            if row[1] != 0.0 and row[2] > 0.0:
                neutron.write(
                    f"{q:.12g} {row[1] * neutron_scale:.12g} {row[2] * neutron_scale:.12g}\n"
                )
            if row[3] != 0.0 and row[4] > 0.0:
                f_si = xray_form_factor(xray_rows[0], q)
                f_o = xray_form_factor(xray_rows[1], q)
                mean_f = (f_si + 2.0 * f_o) / 3.0
                mean_f2 = (f_si * f_si + 2.0 * f_o * f_o) / 3.0
                scale = mean_f2 / (mean_f * mean_f)
                xray.write(f"{q:.12g} {row[3] * scale:.12g} {row[4] * scale:.12g}\n")
    return neutron_path, xray_path


def write_config(output: Path, structure: Path, targets, tables, moves: int, seed: int):
    config = output / "run.toml"
    config.write_text(
        f'''[system]
structure = "{structure}"
format = "lammps"
types = {{ "1" = "Si", "2" = "O" }}

[data]
[data.neutron_sq]
file = "{targets[0]}"
sigma_column = 3
convention = "iq"
fit_min = 0.1
fit_max = 29.95
scattering_lengths = {{ Si = 4.1491, O = 5.803 }}

[data.xray_sq]
file = "{targets[1]}"
sigma_column = 3
convention = "iq"
fit_min = 0.7
fit_max = 19.95

[rmc]
max_moves = {moves}
max_step = 0.1
checkpoint_every = 1000000000
seed = {seed}
print_every = {moves}
target_acceptance = 0.25
adjust_step_every = 1000

[sq]
qmin = 0.0
qmax = 30.0
nq = 600
lorch = false
rdf_cutoff = 22.41
rdf_nbins = 747

[potential]
cutoff = 12.0
weight = 1.0

[[potential.tabulated]]
pair = "Si-Si"
file = "{tables['Si-Si']}"

[[potential.tabulated]]
pair = "Si-O"
file = "{tables['Si-O']}"

[[potential.tabulated]]
pair = "O-O"
file = "{tables['O-O']}"

[epsr]
mode = "pure"
iterations = 1
feedback = 0.9
smooth_sigma = 0.03
moves_per_iteration = {moves}
temperature = 0.025852
min_r = 1.0
'''
    )
    return config


def extract_float(pattern: str, text: str, label: str):
    match = re.search(pattern, text)
    if match is None:
        raise ValueError(f"could not extract {label} from rsmith log")
    return float(match.group(1))


def extract_dataset_residuals(text: str):
    matches = re.findall(
        r"Post-MC dataset (\d+): squared residual =\s*([+\-0-9.eE]+), "
        r"RMS =\s*([+\-0-9.eE]+) \((\d+) Q points\)",
        text,
    )
    if len(matches) != 2:
        raise ValueError("expected two post-MC dataset residual records")
    labels = {"1": "xray", "2": "neutron"}
    return {
        labels[index]: {
            "squared_residual": float(squared),
            "rms": float(rms),
            "q_points": int(points),
        }
        for index, squared, rms, points in matches
    }


parser = argparse.ArgumentParser()
parser.add_argument("--binary", type=Path)
parser.add_argument("--moves", type=int, default=6000)
parser.add_argument("--seed", type=int, default=20260802)
parser.add_argument("--force", action="store_true")
args = parser.parse_args()
if args.moves <= 0:
    raise SystemExit("--moves must be positive")

case_root = Path(__file__).resolve().parents[1]
repo_root = case_root.parents[1]
source = case_root / "reference/local/upstream"
native = case_root / "results/native-zero-move"
structure = case_root / "results/rsmith-zero-move/dtbsilica.data"
for required in (source / "si.ato", source / "o.ato", source / "DTBsilica.pcof", native / "DTBsilica.EPSR.o01", native / "DTBsilica.EPSR.u01", structure):
    if not required.is_file():
        raise SystemExit(f"missing prerequisite {required}; run the local import and zero-move gates first")

output = case_root / "results/reference-potential-smoke"
if output.exists():
    if not args.force:
        raise SystemExit(f"output exists: {output} (pass --force to replace it)")
    shutil.rmtree(output)
output.mkdir(parents=True)

atoms = {
    "Si": atom_parameters(source / "si.ato", "Si"),
    "O": atom_parameters(source / "o.ato", "O"),
}
controls = {
    "power": pcof_setting(source / "DTBsilica.pcof", "power"),
    "rminpt": pcof_setting(source / "DTBsilica.pcof", "rminpt"),
    "rmaxpt": pcof_setting(source / "DTBsilica.pcof", "rmaxpt"),
    "refpotfac": pcof_setting(source / "DTBsilica.pcof", "refpotfac"),
}
pairs = pair_parameters(atoms)
curves = curve_metrics(native / "DTBsilica.EPSR.o01", pairs, controls)
for label, metric in curves.items():
    if metric["rms_over_dynamic_range"] > CURVE_LIMIT:
        raise SystemExit(f"{label} reference-potential parity failed: {metric}")

tables = write_potential_tables(output, pairs, controls)
targets = write_faber_ziman_targets(source, native, output)
config = write_config(output, structure, targets, tables, args.moves, args.seed)
binary = args.binary or repo_root / "target/release/rsmith"
if not binary.is_file():
    raise SystemExit(f"rsmith release binary not found: {binary}")

started = time.perf_counter()
completed = subprocess.run(
    [str(binary), str(config), "--output-dir", str(output), "--quiet"],
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    env={**os.environ, "RAYON_NUM_THREADS": "1"},
    check=False,
)
wall_seconds = time.perf_counter() - started
(output / "driver.log").write_text(completed.stdout)
if completed.returncode != 0:
    raise SystemExit(f"rsmith failed with exit code {completed.returncode}; see {output / 'driver.log'}")

log = (output / "rsmith.log").read_text()
initial_energy_ev = extract_float(r"Initial potential energy =\s*([+\-0-9.eE]+) eV", log, "initial energy")
native_report = (native / "DTBsilica.EPSR.out").read_text()
native_energy_kj_per_mol = extract_float(
    r"Intermolecular energy and U_pol\s*([+\-0-9.eE]+)", native_report, "native energy"
)
rsmith_energy_kj_per_mol = initial_energy_ev * KJ_PER_MOL_PER_EV / 6000.0
energy_relative_difference = abs(rsmith_energy_kj_per_mol - native_energy_kj_per_mol) / abs(native_energy_kj_per_mol)
if energy_relative_difference > ENERGY_LIMIT:
    raise SystemExit(
        f"reference-potential energy parity failed: native={native_energy_kj_per_mol}, "
        f"rsmith={rsmith_energy_kj_per_mol}, relative={energy_relative_difference}"
    )

summary = {
    "status": "passed_reference_potential_and_joint_epsr_smoke",
    "source": "local licensed EPSR26 DTBsilicaNX tutorial",
    "mixing": {
        "epsilon": "geometric",
        "sigma": "arithmetic",
        "charge": "product",
    },
    "controls": controls,
    "pairs": pairs,
    "native_curve_parity": curves,
    "energy_parity": {
        "native_kj_mol_per_epsr_molecule": native_energy_kj_per_mol,
        "rsmith_kj_mol_per_atom": rsmith_energy_kj_per_mol,
        "relative_difference": energy_relative_difference,
    },
    "smoke": {
        "seed": args.seed,
        "iterations": 1,
        "moves": args.moves,
        "threads": 1,
        "wall_seconds": wall_seconds,
        "initial_energy_ev_total": initial_energy_ev,
        "post_mc_residuals": extract_dataset_residuals(log),
        "summed_squared_residual": extract_float(
            r"EPSR iter 1: chi2 =\s*([+\-0-9.eE]+)", log, "summed squared residual"
        ),
        "acceptance_percent": extract_float(r"acceptance =\s*([+\-0-9.eE]+)%", log, "acceptance"),
        "max_ep_update_ev": extract_float(r"max \|.EP\| =\s*([+\-0-9.eE]+) eV", log, "EP update"),
    },
}
(output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(json.dumps(summary, indent=2, sort_keys=True))
print(f"Reference-potential smoke output: {output}")
