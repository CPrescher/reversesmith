#!/usr/bin/env python3
"""Run rsmith on the exact SrTiO3 tutorial configuration and compare outputs."""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
import re
import subprocess
import tomllib
from pathlib import Path


SCATTERING_LENGTHS = {"O": 5.803, "Ti": -3.438, "Sr": 7.02}
CONCENTRATIONS = {"O": 0.6, "Ti": 0.2, "Sr": 0.2}


def read_space_curve(path: Path):
    rows = []
    for raw in path.read_text().splitlines():
        if raw.strip() and not raw.lstrip().startswith("#"):
            fields = raw.split()
            try:
                row = tuple(float(value) for value in fields)
            except ValueError:
                continue
            if len(row) >= 2:
                rows.append(row)
    return rows


def read_csv(path: Path):
    with path.open(newline="") as stream:
        return [tuple(float(value) for value in row) for row in csv.reader(stream) if row and not row[0].startswith("#")]


def interpolate(curve, x, column=1):
    xs = [row[0] for row in curve]
    hi = bisect.bisect_left(xs, x)
    if hi < len(curve) and abs(curve[hi][0] - x) < 1e-12:
        return curve[hi][column]
    if hi == 0 or hi == len(curve):
        raise ValueError(f"{x} is outside the candidate grid")
    lo = hi - 1
    fraction = (x - curve[lo][0]) / (curve[hi][0] - curve[lo][0])
    return curve[lo][column] + fraction * (curve[hi][column] - curve[lo][column])


def metrics(reference, candidate, ref_column, fit_min, fit_max, candidate_column=1):
    selected = [row for row in reference if fit_min <= row[0] <= fit_max]
    expected = [row[ref_column] for row in selected]
    residual = [interpolate(candidate, row[0], candidate_column) - row[ref_column] for row in selected]
    dynamic_range = max(expected) - min(expected)
    rms = math.sqrt(sum(value * value for value in residual) / len(residual))
    maximum = max(abs(value) for value in residual)
    return {
        "n_points": len(selected),
        "rms_error": rms,
        "max_absolute_error": maximum,
        "rms_over_dynamic_range": rms / dynamic_range,
        "max_over_dynamic_range": maximum / dynamic_range,
    }


def rebin(curve, factor):
    rebinned = []
    for start in range(0, len(curve) - factor + 1, factor):
        block = curve[start : start + factor]
        n_columns = len(block[0])
        rebinned.append(tuple(sum(row[column] for row in block) / factor for column in range(n_columns)))
    return rebinned


def convert_rmc7(path: Path, output: Path):
    lines = path.read_text().splitlines()
    atom_count = int(next(line for line in lines if line.startswith("Number of atoms:")).split()[-1])
    cell_line = next(line for line in lines if line.startswith("Cell (Ang/deg):"))
    box = [float(value) for value in cell_line.split(":", 1)[1].split()[:3]]
    atom_start = lines.index("Atoms:") + 1
    atoms = []
    for line in lines[atom_start : atom_start + atom_count]:
        fields = line.split()
        atoms.append((fields[1], tuple(float(value) for value in fields[3:6])))
    species = []
    for symbol, _ in atoms:
        if symbol not in species:
            species.append(symbol)
    type_ids = {symbol: index + 1 for index, symbol in enumerate(species)}
    with output.open("w") as stream:
        stream.write("official RMCProfile SrTiO3 tutorial configuration\n\n")
        stream.write(f"{atom_count} atoms\n{len(species)} atom types\n\n")
        stream.write(f"0.0 {box[0]:.12f} xlo xhi\n0.0 {box[1]:.12f} ylo yhi\n0.0 {box[2]:.12f} zlo zhi\n\n")
        stream.write("Atoms # charge\n\n")
        for atom_id, (symbol, fractional) in enumerate(atoms, 1):
            xyz = [fractional[axis] * box[axis] for axis in range(3)]
            stream.write(f"{atom_id} {type_ids[symbol]} 0.0 {xyz[0]:.12f} {xyz[1]:.12f} {xyz[2]:.12f}\n")
    return type_ids


parser = argparse.ArgumentParser()
parser.add_argument("--native-run-dir", type=Path, required=True)
parser.add_argument("--binary", type=Path)
args = parser.parse_args()
case_root = Path(__file__).resolve().parents[1]
repo_root = case_root.parents[1]
binary = args.binary or repo_root / "target/release/rsmith"
results = case_root / "results/rsmith"
results.mkdir(parents=True, exist_ok=True)

type_ids = convert_rmc7(args.native_run_dir / "SRTIO3_5K.rmc7", results / "srtio3.data")
input_fq = read_space_curve(args.native_run_dir / "srtio3_wk_5k_rmc.fq")
fq_path = results / "srtio3-fq.dat"
fq_path.write_text("# Q F(Q) sigma\n" + "".join(f"{q:.12g} {value:.12g} 1.0\n" for q, value in input_fq))
types_toml = ", ".join(f'"{type_id}" = "{symbol}"' for symbol, type_id in type_ids.items())
config = results / "forward.toml"
config.write_text(
    f'''[system]
structure = "{results / 'srtio3.data'}"
format = "lammps"
types = {{ {types_toml} }}

[data]
[data.neutron_sq]
file = "{fq_path}"
sigma_column = 3
convention = "iq"
fit_min = 0.691
fit_max = 39.931

[rmc]
max_moves = 0
seed = 20260802

[sq]
qmin = 0.02
qmax = 39.96
nq = 1997
lorch = false
rdf_cutoff = 15.6
rdf_nbins = 780
'''
)
subprocess.run([str(binary), str(config), "--output-dir", str(results), "--quiet"], check=True)

native_fq = read_csv(args.native_run_dir / "srtio3_wk_5k_rmc.fq_1.csv")
rsmith_fq = read_space_curve(results / "start_neutron_sq.dat")
b_average = sum(CONCENTRATIONS[name] * SCATTERING_LENGTHS[name] for name in CONCENTRATIONS)
# RMCProfile's F(Q) is 0.01*<b>^2*[S(Q)-1], and its calculated curve is
# sampled on the 0.02 A^-1 model grid that is offset by +0.009 A^-1 from
# the experimental x values written to the comparison CSV.
rsmith_rmc_fq = [(row[0] - 0.009, 0.01 * b_average * b_average * row[1]) for row in rsmith_fq]
summary = {
    "native_version": "RMCProfile7b.35",
    "atom_count": 2880,
    "fq": metrics(native_fq, rsmith_rmc_fq, 2, 0.691, 39.931),
}

native_gr = read_csv(args.native_run_dir / "SRTIO3_5K.rmc7_grpartials.csv")
rsmith_gr_raw = read_space_curve(results / "start_gr.dat")
# RMCProfile labels histogram bin i by its upper edge (i+1)*dr, whereas
# rsmith labels the identical bin by its centre (i+0.5)*dr.
rsmith_gr = [(row[0] + 0.01, *row[1:]) for row in rsmith_gr_raw]
pair_names = ("O-O", "O-Ti", "O-Sr", "Ti-Ti", "Ti-Sr", "Sr-Sr")
summary["partial_gr"] = {
    pair: metrics(native_gr, rsmith_gr, column, 1.5, 15.59, column)
    for column, pair in enumerate(pair_names, 1)
}
native_gr_rebinned = rebin(native_gr, 5)
rsmith_gr_rebinned = rebin(rsmith_gr, 5)
summary["partial_gr_rebinned_0.1A"] = {
    pair: metrics(native_gr_rebinned, rsmith_gr_rebinned, column, 1.5, 15.5, column)
    for column, pair in enumerate(pair_names, 1)
}

weights = []
species = ("O", "Ti", "Sr")
for first_index, first in enumerate(species):
    for second in species[first_index:]:
        factor = 1.0 if first == second else 2.0
        weights.append(
            factor
            * CONCENTRATIONS[first]
            * CONCENTRATIONS[second]
            * SCATTERING_LENGTHS[first]
            * SCATTERING_LENGTHS[second]
            / (b_average * b_average)
        )
rsmith_total_g = [
    (row[0], 0.01 * b_average * b_average * (sum(weight * value for weight, value in zip(weights, row[1:])) - 1.0))
    for row in rsmith_gr
]
native_total_g = read_csv(args.native_run_dir / "srtio3_wk_5k_rmc.gr_1.csv")
summary["total_gr"] = metrics(native_total_g, rsmith_total_g, 2, 1.5, 15.59)
summary["total_gr_rebinned_0.1A"] = metrics(
    rebin(native_total_g, 5), rebin(rsmith_total_g, 5), 2, 1.5, 15.5
)

expected_path = case_root / "expected/native-forward.toml"
expected = tomllib.loads(expected_path.read_text())
if summary["fq"]["rms_over_dynamic_range"] > expected["fq_rms_over_dynamic_range_max"]:
    raise AssertionError("native RMCProfile F(Q) RMS regression limit exceeded")
if summary["fq"]["max_over_dynamic_range"] > expected["fq_max_over_dynamic_range_max"]:
    raise AssertionError("native RMCProfile F(Q) maximum-error regression limit exceeded")
total_rebinned = summary["total_gr_rebinned_0.1A"]
if total_rebinned["rms_over_dynamic_range"] > expected["total_gr_rebinned_rms_over_dynamic_range_max"]:
    raise AssertionError("native RMCProfile rebinned total G(r) RMS regression limit exceeded")
if total_rebinned["max_over_dynamic_range"] > expected["total_gr_rebinned_max_over_dynamic_range_max"]:
    raise AssertionError("native RMCProfile rebinned total G(r) maximum-error limit exceeded")
for pair, pair_metrics in summary["partial_gr_rebinned_0.1A"].items():
    if pair_metrics["rms_over_dynamic_range"] > expected["partial_gr_rebinned_rms_over_dynamic_range_max"]:
        raise AssertionError(f"native RMCProfile {pair} rebinned RDF RMS regression limit exceeded")
    if pair_metrics["max_over_dynamic_range"] > expected["partial_gr_rebinned_max_over_dynamic_range_max"]:
        raise AssertionError(f"native RMCProfile {pair} rebinned RDF maximum-error limit exceeded")

(results / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(json.dumps(summary, indent=2, sort_keys=True))
print(f"passed frozen regression checks from {expected_path}")
