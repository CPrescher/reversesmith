#!/usr/bin/env python3
"""Run and score rsmith against a zero-move EPSR26 DTBsilicaNX calculation."""

from __future__ import annotations

import argparse
import bisect
import hashlib
import json
import math
import subprocess
from pathlib import Path

import numpy as np


PAIR_NAMES = ("Si-Si", "Si-O", "O-O")
NATIVE_PAIR_COLUMNS = (1, 3, 5)


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
    if not rows:
        raise ValueError(f"no numeric curve rows in {path}")
    return rows


def interpolate(curve, x, column):
    xs = [row[0] for row in curve]
    hi = bisect.bisect_left(xs, x)
    if hi < len(curve) and abs(xs[hi] - x) < 1.0e-9:
        return curve[hi][column]
    if hi == 0 or hi == len(curve):
        raise ValueError(f"{x} outside candidate grid")
    lo = hi - 1
    fraction = (x - xs[lo]) / (xs[hi] - xs[lo])
    return curve[lo][column] + fraction * (curve[hi][column] - curve[lo][column])


def metrics(reference, candidate, fit_min, fit_max, ref_column=1, candidate_column=1):
    selected = [row for row in reference if fit_min <= row[0] <= fit_max]
    expected = np.asarray([row[ref_column] for row in selected])
    calculated = np.asarray(
        [interpolate(candidate, row[0], candidate_column) for row in selected]
    )
    residual = calculated - expected
    dynamic_range = float(np.ptp(expected))
    rms = float(np.sqrt(np.mean(residual * residual)))
    maximum = float(np.max(np.abs(residual)))

    design = np.column_stack((calculated, np.ones_like(calculated)))
    scale, offset = np.linalg.lstsq(design, expected, rcond=None)[0]
    affine_residual = scale * calculated + offset - expected
    affine_rms = float(np.sqrt(np.mean(affine_residual * affine_residual)))
    scale_only = float(np.dot(calculated, expected) / np.dot(calculated, calculated))
    scale_only_residual = scale_only * calculated - expected
    scale_only_rms = float(np.sqrt(np.mean(scale_only_residual * scale_only_residual)))
    denominator = dynamic_range if dynamic_range > 1.0e-30 else 1.0
    return {
        "n_points": len(selected),
        "rms_error": rms,
        "max_absolute_error": maximum,
        "dynamic_range": dynamic_range,
        "rms_over_dynamic_range": rms / denominator,
        "max_over_dynamic_range": maximum / denominator,
        "best_affine_scale": float(scale),
        "best_affine_offset": float(offset),
        "affine_rms_over_dynamic_range": affine_rms / denominator,
        "best_scale_only": scale_only,
        "scale_only_rms_over_dynamic_range": scale_only_rms / denominator,
    }


def block_average(curve, width, columns):
    start = math.ceil(curve[0][0] / width) * width
    stop = math.floor(curve[-1][0] / width) * width
    points = []
    index = 0
    x = start
    while x <= stop + 1.0e-10:
        upper = x + width
        rows = []
        while index < len(curve) and curve[index][0] < upper - 1.0e-10:
            if curve[index][0] >= x - 1.0e-10:
                rows.append(curve[index])
            index += 1
        if rows:
            points.append(
                (
                    x + 0.5 * width,
                    *(sum(row[c] for row in rows) / len(rows) for c in columns),
                )
            )
        x = upper
    return points


def convert_ato(path: Path, output: Path):
    lines = path.read_text().splitlines()
    atom_count, box = int(lines[0].split()[0]), float(lines[0].split()[1])
    records = lines[2:]
    if len(records) < 5 * atom_count:
        raise ValueError("truncated EPSR .ato file")
    atoms = []
    for index in range(atom_count):
        coordinate_fields = records[5 * index].split()
        species = records[5 * index + 1].split()[0]
        position = tuple(
            (float(coordinate_fields[axis]) + 0.5 * box) % box for axis in (1, 2, 3)
        )
        atoms.append((species, position))
    species_counts = {
        name: sum(atom[0] == name for atom in atoms) for name in ("Si", "O")
    }
    if species_counts != {"Si": 2000, "O": 4000}:
        raise ValueError(f"unexpected DTBsilicaNX composition: {species_counts}")

    with output.open("w") as stream:
        stream.write("EPSR26 DTBsilicaNX zero-move configuration\n\n")
        stream.write(f"{atom_count} atoms\n2 atom types\n\n")
        for axis in "xyz":
            stream.write(f"0.0 {box:.12f} {axis}lo {axis}hi\n")
        stream.write("\nAtoms # charge\n\n")
        for atom_id, (species, position) in enumerate(atoms, 1):
            type_id = 1 if species == "Si" else 2
            stream.write(
                f"{atom_id} {type_id} 0.0 "
                f"{position[0]:.12f} {position[1]:.12f} {position[2]:.12f}\n"
            )
    positions = np.asarray([position for _, position in atoms], dtype=np.float64)
    type_ids = np.asarray(
        [0 if species == "Si" else 1 for species, _ in atoms], dtype=np.int8
    )
    return atom_count, box, species_counts, positions, type_ids


def direct_partial_rdf_oracle(positions, type_ids, box, cutoff, bins):
    histograms = np.zeros((3, bins), dtype=np.int64)
    inverse_dr = bins / cutoff
    for index in range(len(positions) - 1):
        displacement = positions[index + 1 :] - positions[index]
        displacement -= box * np.rint(displacement / box)
        distance = np.sqrt(np.einsum("ij,ij->i", displacement, displacement))
        selected = distance < cutoff
        if not np.any(selected):
            continue
        other_types = type_ids[index + 1 :][selected]
        pair_indices = (
            np.zeros_like(other_types)
            if type_ids[index] == 0
            else np.where(other_types == 0, 1, 2)
        )
        if type_ids[index] == 0:
            pair_indices = np.where(other_types == 0, 0, 1)
        bin_indices = np.floor(distance[selected] * inverse_dr).astype(np.int64)
        for pair_index in range(3):
            histograms[pair_index] += np.bincount(
                bin_indices[pair_indices == pair_index], minlength=bins
            )

    dr = cutoff / bins
    radius = (np.arange(bins) + 0.5) * dr
    volume = box**3
    counts = (int(np.sum(type_ids == 0)), int(np.sum(type_ids == 1)))
    pair_types = ((0, 0), (0, 1), (1, 1))
    partials = []
    shell = 4.0 * math.pi * radius**2 * dr
    for pair_index, (a, b) in enumerate(pair_types):
        norm = counts[a] * (counts[b] / volume) * shell
        factor = 2.0 if a == b else 1.0
        partials.append(factor * histograms[pair_index] / norm)
    return list(zip(radius.tolist(), *(partial.tolist() for partial in partials)))


def write_target(path, native_target, column):
    path.write_text(
        "# Q i(Q)\n"
        + "".join(f"{row[0]:.12g} {row[column]:.12g}\n" for row in native_target)
    )


def epsr_neutron_weights(path: Path):
    triples = []
    for line in path.read_text().splitlines():
        try:
            values = tuple(float(value) for value in line.split())
        except ValueError:
            continue
        if len(values) == 3:
            triples.append(values)
    if len(triples) != 3:
        raise ValueError(f"expected three pair weights in {path}, found {len(triples)}")
    return np.asarray([values[2] for values in triples])


def epsr_xray_form_factor_parameters(path: Path):
    parameter_rows = []
    for line in path.read_text().splitlines():
        try:
            values = tuple(float(value) for value in line.split())
        except ValueError:
            continue
        if len(values) == 11:
            parameter_rows.append(values)
    if len(parameter_rows) != 2:
        raise ValueError(f"expected Si and O X-ray parameters in {path}")
    return parameter_rows


def epsr_xray_pair_weights(q_values, parameter_path: Path):
    parameters = epsr_xray_form_factor_parameters(parameter_path)
    s_squared = (np.asarray(q_values) / (4.0 * math.pi)) ** 2

    def form_factor(values):
        constant = values[0]
        a = np.asarray(values[1::2])
        b = np.asarray(values[2::2])
        return constant + np.sum(a[:, None] * np.exp(-b[:, None] * s_squared), axis=0)

    silicon, oxygen = (form_factor(values) for values in parameters)
    c_si, c_o = 1.0 / 3.0, 2.0 / 3.0
    # EPSR DTBsilicaNX XWTS uses normtot=2: single-atom scattering <f^2>,
    # rather than the Faber-Ziman <f>^2 denominator used by rsmith outputs.
    denominator = c_si * silicon**2 + c_o * oxygen**2
    return np.column_stack(
        (
            c_si * c_si * silicon**2 / denominator,
            2.0 * c_si * c_o * silicon * oxygen / denominator,
            c_o * c_o * oxygen**2 / denominator,
        )
    )


def weighted_total(partial_iq, weights):
    q_values = np.asarray([row[0] for row in partial_iq])
    partial_values = np.asarray([row[1:4] for row in partial_iq])
    if np.ndim(weights) == 1:
        values = partial_values @ weights
    else:
        values = np.sum(partial_values * weights, axis=1)
    return list(zip(q_values.tolist(), values.tolist()))


def write_curve(path: Path, curve, label):
    path.write_text(
        f"# Q(1/A) {label}\n"
        + "".join(f"{row[0]:.12g} {row[1]:.12g}\n" for row in curve)
    )


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


parser = argparse.ArgumentParser()
parser.add_argument("--binary", type=Path)
parser.add_argument("--native-run-dir", type=Path)
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
repo_root = case_root.parents[1]
binary = args.binary or repo_root / "target/release/rsmith"
native = args.native_run_dir or case_root / "results/native-zero-move"
source = case_root / "reference/local/upstream"
if not binary.is_file():
    raise SystemExit(f"rsmith binary not found: {binary}; run cargo build --release")
for required in (source / "DTBsilica.ato", native / "DTBsilica.EPSR.f01"):
    if not required.is_file():
        raise SystemExit(f"required local file is missing: {required}")

results = case_root / "results/rsmith-zero-move"
results.mkdir(parents=True, exist_ok=True)
atom_count, box, composition, positions, type_ids = convert_ato(
    source / "DTBsilica.ato", results / "dtbsilica.data"
)
native_target = read_curve(native / "DTBsilica.EPSR.t01")
write_target(results / "native-neutron-iq.dat", native_target, 1)
write_target(results / "native-xray-iq.dat", native_target, 3)

rdf_cutoff = math.floor((0.5 * box - 1.0e-8) / 0.03) * 0.03
rdf_bins = round(rdf_cutoff / 0.03)
config = results / "forward.toml"
config.write_text(
    f'''[system]
structure = "{results / "dtbsilica.data"}"
format = "lammps"
types = {{ "1" = "Si", "2" = "O" }}

[data]
[data.neutron_sq]
file = "{results / "native-neutron-iq.dat"}"
sigma = 1.0
convention = "iq"
fit_min = 0.1
fit_max = 29.95

[data.xray_sq]
file = "{results / "native-xray-iq.dat"}"
sigma = 1.0
convention = "iq"
fit_min = 0.1
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
rdf_nbins = {rdf_bins}
'''
)
subprocess.run(
    [str(binary), str(config), "--output-dir", str(results), "--quiet"], check=True
)

native_partial_iq = read_curve(native / "DTBsilica.EPSR.f01")
native_gr = read_curve(native / "DTBsilica.EPSR.g01")
native_totals = read_curve(native / "DTBsilica.EPSR.u01")
rsmith_partial_sq = read_curve(results / "start_partial_sq.dat")
rsmith_partial_iq = [
    tuple([row[0], *(value - 1.0 for value in row[1:])]) for row in rsmith_partial_sq
]
rsmith_gr = read_curve(results / "start_gr.dat")
rsmith_neutron = read_curve(results / "start_neutron_sq.dat")
rsmith_xray = read_curve(results / "start_xray_sq.dat")
neutron_weights = epsr_neutron_weights(source / "DTBsilica.NWTStot.wts")
xray_weights = epsr_xray_pair_weights(
    [row[0] for row in rsmith_partial_iq], source / "DTBsilica.XWTS.wts"
)
rsmith_epsr_neutron = weighted_total(rsmith_partial_iq, neutron_weights)
rsmith_epsr_xray = weighted_total(rsmith_partial_iq, xray_weights)
write_curve(
    results / "rsmith-epsr-normalized-neutron-iq.dat",
    rsmith_epsr_neutron,
    "i_neutron(Q)",
)
write_curve(
    results / "rsmith-epsr-normalized-xray-iq.dat", rsmith_epsr_xray, "i_xray(Q)"
)

native_gr_rebinned = block_average(native_gr, 0.12, NATIVE_PAIR_COLUMNS)
rsmith_gr_rebinned = block_average(rsmith_gr, 0.12, (1, 2, 3))
oracle_gr = direct_partial_rdf_oracle(positions, type_ids, box, rdf_cutoff, rdf_bins)

summary = {
    "status": "diagnostic_forward_baseline_no_frozen_equivalence_limits",
    "reference": {
        "program": "EPSR26",
        "case": "DTBsilicaNX",
        "ato_sha256": sha256(source / "DTBsilica.ato"),
        "input_sha256": sha256(source / "DTBsilica.EPSR.inp"),
    },
    "system": {
        "atom_count": atom_count,
        "composition": composition,
        "box_angstrom": box,
        "number_density_atoms_per_a3": atom_count / box**3,
        "rdf_cutoff_angstrom": rdf_cutoff,
        "rdf_bins": rdf_bins,
    },
    "totals": {
        "neutron_iq": {
            "rsmith_default_faber_ziman": metrics(
                native_totals, rsmith_neutron, 0.5, 29.95, 1, 1
            ),
            "epsr_nrtype5_from_rsmith_partials": metrics(
                native_totals, rsmith_epsr_neutron, 0.5, 29.95, 1, 1
            ),
            "epsr_pair_weights": neutron_weights.tolist(),
        },
        "xray_iq": {
            "rsmith_default_faber_ziman": metrics(
                native_totals, rsmith_xray, 0.5, 29.95, 3, 1
            ),
            "epsr_single_atom_normalized_from_rsmith_partials": metrics(
                native_totals, rsmith_epsr_xray, 0.5, 29.95, 3, 1
            ),
        },
    },
    "partial_iq": {
        pair: metrics(
            native_partial_iq, rsmith_partial_iq, 0.5, 29.95, native_column, index + 1
        )
        for index, (pair, native_column) in enumerate(
            zip(PAIR_NAMES, NATIVE_PAIR_COLUMNS)
        )
    },
    "partial_gr_rebinned_0.12A": {
        pair: metrics(
            native_gr_rebinned,
            rsmith_gr_rebinned,
            1.2,
            rdf_cutoff - 0.12,
            index + 1,
            index + 1,
        )
        for index, pair in enumerate(PAIR_NAMES)
    },
    "rsmith_gr_vs_independent_oracle": {
        pair: metrics(
            oracle_gr, rsmith_gr, 0.5, rdf_cutoff - 0.03, index + 1, index + 1
        )
        for index, pair in enumerate(PAIR_NAMES)
    },
}

oracle_limit = 2.0e-6
oracle_failures = [
    f"{pair}: {values['rms_over_dynamic_range']:.6g}"
    for pair, values in summary["rsmith_gr_vs_independent_oracle"].items()
    if values["rms_over_dynamic_range"] > oracle_limit
]
summary["integrity_checks"] = {
    "independent_rdf_oracle_limit": oracle_limit,
    "independent_rdf_oracle_passed": not oracle_failures,
}
(results / "summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(summary, indent=2, sort_keys=True))
if oracle_failures:
    raise SystemExit(
        "independent RDF oracle failed:\n  " + "\n  ".join(oracle_failures)
    )
print("completed deterministic EPSR26/rsmith DTBsilicaNX forward comparison")
