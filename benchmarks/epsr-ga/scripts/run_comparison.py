#!/usr/bin/env python3
"""Convert the EPSR Ga configuration, run rsmith at zero moves, and compare."""

from __future__ import annotations

import argparse
import bisect
import json
import math
import subprocess
import tomllib
from pathlib import Path

import numpy as np


def read_curve(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        try:
            row = tuple(float(value) for value in line.split())
        except ValueError:
            continue
        if len(row) >= 2:
            rows.append(row)
    return rows


def interpolate(curve, x, column=1):
    xs = [row[0] for row in curve]
    hi = bisect.bisect_left(xs, x)
    if hi < len(curve) and abs(xs[hi] - x) < 1.0e-10:
        return curve[hi][column]
    if hi == 0 or hi == len(curve):
        raise ValueError(f"{x} outside candidate grid")
    lo = hi - 1
    fraction = (x - xs[lo]) / (xs[hi] - xs[lo])
    return curve[lo][column] + fraction * (curve[hi][column] - curve[lo][column])


def metrics(reference, candidate, fit_min, fit_max):
    selected = [row for row in reference if fit_min <= row[0] <= fit_max]
    expected = [row[1] for row in selected]
    residual = [interpolate(candidate, row[0]) - row[1] for row in selected]
    dynamic_range = max(expected) - min(expected)
    rms = math.sqrt(sum(value * value for value in residual) / len(residual))
    maximum = max(abs(value) for value in residual)
    return {
        "n_points": len(selected),
        "rms_error": rms,
        "max_absolute_error": maximum,
        "dynamic_range": dynamic_range,
        "rms_over_dynamic_range": rms / dynamic_range,
        "max_over_dynamic_range": maximum / dynamic_range,
    }


def block_average(curve, width):
    start = math.ceil(curve[0][0] / width) * width
    stop = math.floor(curve[-1][0] / width) * width
    points = []
    index = 0
    x = start
    while x <= stop + 1.0e-10:
        upper = x + width
        values = []
        while index < len(curve) and curve[index][0] < upper - 1.0e-10:
            if curve[index][0] >= x - 1.0e-10:
                values.append(curve[index][1])
            index += 1
        if values:
            points.append((x + 0.5 * width, sum(values) / len(values)))
        x = upper
    return points


def convert_ato(path: Path, output: Path):
    lines = path.read_text().splitlines()
    header = lines[0].split()
    atom_count = int(header[0])
    box = float(header[1])
    records = lines[2:]
    if len(records) < 5 * atom_count:
        raise ValueError("truncated EPSR .ato file")
    atoms = []
    for index in range(atom_count):
        fields = records[5 * index].split()
        species = records[5 * index + 1].split()[0]
        xyz = tuple((float(fields[axis]) + 0.5 * box) % box for axis in (1, 2, 3))
        atoms.append((species, xyz))
    if {species for species, _ in atoms} != {"Ga"}:
        raise ValueError("this converter is intentionally specific to monatomic LiquidGa50C")
    with output.open("w") as stream:
        stream.write("EPSR26 LiquidGa50C zero-move configuration\n\n")
        stream.write(f"{atom_count} atoms\n1 atom types\n\n")
        for axis in "xyz":
            stream.write(f"0.0 {box:.12f} {axis}lo {axis}hi\n")
        stream.write("\nAtoms # charge\n\n")
        for atom_id, (_, xyz) in enumerate(atoms, 1):
            stream.write(f"{atom_id} 1 0.0 {xyz[0]:.12f} {xyz[1]:.12f} {xyz[2]:.12f}\n")
    return atom_count, box, np.asarray([xyz for _, xyz in atoms])


def direct_histogram_oracle(positions, box, cutoff, nbins):
    """Independent NumPy O(N^2) minimum-image pair histogram."""
    histogram = np.zeros(nbins, dtype=np.int64)
    dr = cutoff / nbins
    for index in range(len(positions) - 1):
        displacement = positions[index + 1 :] - positions[index]
        displacement -= box * np.rint(displacement / box)
        distance = np.sqrt(np.einsum("ij,ij->i", displacement, displacement))
        bins = np.floor(distance[distance < cutoff] / dr).astype(np.int64)
        histogram += np.bincount(bins, minlength=nbins)
    radius = (np.arange(nbins) + 0.5) * dr
    number_density = len(positions) / box**3
    normalisation = len(positions) * number_density * 4.0 * math.pi * radius**2 * dr
    gr = 2.0 * histogram / normalisation
    return list(zip(radius.tolist(), gr.tolist()))


parser = argparse.ArgumentParser()
parser.add_argument("--binary", type=Path)
parser.add_argument("--native-run-dir", type=Path)
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
repo_root = case_root.parents[1]
native = args.native_run_dir or case_root / "results/native-zero-move"
binary = args.binary or repo_root / "target/release/rsmith"
if not binary.is_file():
    raise SystemExit(f"rsmith binary not found: {binary}; run cargo build --release")

manifest = tomllib.loads((case_root / "manifest.toml").read_text())
if manifest["tolerances"]["status"] != "preregistered_before_rsmith_result":
    raise SystemExit("forward tolerances are not marked as preregistered")

results = case_root / "results/rsmith-zero-move"
results.mkdir(parents=True, exist_ok=True)
source_ato = case_root / "reference/local/upstream/LiquidGa50C.ato"
if not source_ato.is_file():
    raise SystemExit("local upstream configuration is missing; rerun import_local_reference.sh")
atom_count, box, positions = convert_ato(source_ato, results / "liquid-ga.data")

native_iq = read_curve(native / "LiquidGa50C.EPSR.f01")
data_path = results / "native-partial-iq.dat"
data_path.write_text(
    "# Q S(Q)-1 sigma\n"
    + "".join(f"{row[0]:.12g} {row[1]:.12g} {max(row[2], 1.0e-6):.12g}\n" for row in native_iq)
)

rdf_cutoff = math.floor((0.5 * box - 1.0e-8) / 0.03) * 0.03
rdf_nbins = round(rdf_cutoff / 0.03)
config = results / "forward.toml"
config.write_text(
    f'''[system]
structure = "{results / 'liquid-ga.data'}"
format = "lammps"
types = {{ "1" = "Ga" }}

[data]
[data.neutron_sq]
file = "{data_path}"
sigma_column = 3
convention = "iq"
fit_min = 0.05
fit_max = 29.95

[rmc]
max_moves = 0
seed = 20260802

[sq]
qmin = 0.0
qmax = 30.0
nq = 600
lorch = false
rdf_cutoff = {rdf_cutoff:.12g}
rdf_nbins = {rdf_nbins}
'''
)
subprocess.run([str(binary), str(config), "--output-dir", str(results), "--quiet"], check=True)

rsmith_iq = read_curve(results / "start_neutron_sq.dat")
native_gr = read_curve(native / "LiquidGa50C.EPSR.g01")
rsmith_gr = read_curve(results / "start_gr.dat")
oracle_gr = direct_histogram_oracle(positions, box, rdf_cutoff, rdf_nbins)
native_gr_rebinned = block_average(native_gr, 0.12)
rsmith_gr_rebinned = block_average(rsmith_gr, 0.12)

summary = {
    "native_version": manifest["reference_version"],
    "atom_count": atom_count,
    "box_angstrom": box,
    "rdf_cutoff_angstrom": rdf_cutoff,
    "rdf_nbins": rdf_nbins,
    "oracle_gr": metrics(oracle_gr, rsmith_gr, 0.5, rdf_cutoff - 0.03),
    "gr_rebinned_0.12A": metrics(native_gr_rebinned, rsmith_gr_rebinned, 1.5, rdf_cutoff - 0.12),
    "iq": metrics(native_iq, rsmith_iq, 0.5, 29.95),
}

tolerances = manifest["tolerances"]
checks = (
    ("gr_rebinned_0.12A", "rms_over_dynamic_range", "gr_rebinned_rms_over_dynamic_range_max"),
    ("gr_rebinned_0.12A", "max_over_dynamic_range", "gr_rebinned_max_over_dynamic_range_max"),
    ("iq", "rms_over_dynamic_range", "iq_rms_over_dynamic_range_max"),
    ("iq", "max_over_dynamic_range", "iq_max_over_dynamic_range_max"),
    ("oracle_gr", "rms_over_dynamic_range", "oracle_gr_rms_over_dynamic_range_max"),
    ("oracle_gr", "max_over_dynamic_range", "oracle_gr_max_over_dynamic_range_max"),
)
failures = []
for section, metric, limit_name in checks:
    observed = summary[section][metric]
    limit = tolerances[limit_name]
    if observed > limit:
        failures.append(f"{section}.{metric}={observed:.6g} exceeds {limit_name}={limit:.6g}")

(results / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(json.dumps(summary, indent=2, sort_keys=True))
if failures:
    raise SystemExit("frozen forward checks failed:\n  " + "\n  ".join(failures))
print("passed frozen pre-result EPSR26 forward checks")
