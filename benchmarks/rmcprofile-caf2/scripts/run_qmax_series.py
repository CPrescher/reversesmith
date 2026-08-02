#!/usr/bin/env python3
"""Run and analyse the RMCProfile-inspired CaF2 finite-Q PDF benchmark."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import tomllib
from pathlib import Path

import numpy as np


QMAX_VALUES = (12.0, 14.0, 16.0, 18.0, 20.0, 30.0)
B_CA = 4.70
B_F = 5.654
C_CA = 1.0 / 3.0
C_F = 2.0 / 3.0


def load_columns(path: Path) -> np.ndarray:
    rows = []
    for raw in path.read_text().splitlines():
        if raw.strip() and not raw.lstrip().startswith("#"):
            rows.append([float(value) for value in raw.split()])
    if not rows:
        raise ValueError(f"no numeric rows in {path}")
    return np.asarray(rows, dtype=float)


def write_curve(path: Path, header: str, x: np.ndarray, y: np.ndarray) -> None:
    with path.open("w") as stream:
        stream.write(f"# {header}\n")
        for xi, yi in zip(x, y):
            stream.write(f"{xi:.6f} {yi:.12g}\n")


def inverse_pdf(q: np.ndarray, sq: np.ndarray, r: np.ndarray, qmax: float) -> np.ndarray:
    selected = q <= qmax + 1e-12
    q_use = q[selected]
    fq = q_use * (sq[selected] - 1.0)
    dq = q_use[1] - q_use[0]
    return (2.0 * dq / math.pi) * (np.sin(np.outer(r, q_use)) @ fq)


def error_metrics(reference: np.ndarray, candidate: np.ndarray, selected: np.ndarray) -> dict[str, float | int]:
    residual = candidate[selected] - reference[selected]
    dynamic_range = float(np.ptp(reference[selected]))
    rms = float(np.sqrt(np.mean(residual * residual)))
    maximum = float(np.max(np.abs(residual)))
    return {
        "n_points": int(np.count_nonzero(selected)),
        "rms_error": rms,
        "max_absolute_error": maximum,
        "rms_over_dynamic_range": rms / dynamic_range,
        "max_over_dynamic_range": maximum / dynamic_range,
    }


parser = argparse.ArgumentParser()
parser.add_argument("--cells", type=int, default=15)
parser.add_argument("--lattice", type=float, default=5.462)
parser.add_argument("--uiso-ca", type=float, default=0.005)
parser.add_argument("--uiso-f", type=float, default=0.007)
parser.add_argument("--seed", type=int, default=20260802)
parser.add_argument("--binary", type=Path)
parser.add_argument("--skip-run", action="store_true")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
repo_root = case_root.parents[1]
binary = args.binary or repo_root / "target/release/rsmith"
results = case_root / "results"
results.mkdir(parents=True, exist_ok=True)
structure = case_root / "rsmith" / f"caf2_{args.cells}x{args.cells}x{args.cells}.data"
config = case_root / "rsmith" / "forward.toml"

if not args.skip_run:
    subprocess.run(
        [
            "python3",
            str(case_root / "scripts/generate_caf2_supercell.py"),
            "--cells",
            str(args.cells),
            "--lattice",
            str(args.lattice),
            "--uiso-ca",
            str(args.uiso_ca),
            "--uiso-f",
            str(args.uiso_f),
            "--seed",
            str(args.seed),
            "--output",
            str(structure),
        ],
        check=True,
    )
    config.write_text(
        f'''[system]
structure = "{structure}"
format = "lammps"
types = {{ "1" = "Ca", "2" = "F" }}

[data]
[data.neutron_sq]
file = "{case_root / 'rsmith/neutron-flat.dat'}"
sigma_column = 3
fit_min = 0.01
fit_max = 29.99

[rmc]
max_moves = 0
seed = {args.seed}

[sq]
qmin = 0.01
qmax = 30.0
nq = 3000
lorch = false
rdf_cutoff = 20.0
rdf_nbins = 2000
'''
    )
    subprocess.run(
        [str(binary), str(config), "--output-dir", str(results), "--quiet"],
        check=True,
    )

gr = load_columns(results / "start_gr.dat")
sq_data = load_columns(results / "start_neutron_sq.dat")
r = gr[:, 0]
q = sq_data[:, 0]
sq = sq_data[:, 1]
if gr.shape[1] != 4:
    raise ValueError("expected Ca-Ca, Ca-F, and F-F partial RDF columns")

b_average = C_CA * B_CA + C_F * B_F
weights = np.asarray(
    [
        C_CA * C_CA * B_CA * B_CA,
        2.0 * C_CA * C_F * B_CA * B_F,
        C_F * C_F * B_F * B_F,
    ]
) / (b_average * b_average)
g_total = gr[:, 1:] @ weights
rho = 12.0 / args.lattice**3
histogram_pdf = 4.0 * math.pi * rho * r * (g_total - 1.0)
write_curve(results / "histogram_G.dat", "r [A]  neutron histogram G(r) [1/A]", r, histogram_pdf)

fit_selection = (r >= 1.5) & (r <= 15.0)
series = []
for qmax in QMAX_VALUES:
    finite_q_pdf = inverse_pdf(q, sq, r, qmax)
    write_curve(
        results / f"backtransform_qmax_{qmax:g}.dat",
        f"r [A]  back-transform G(r), Qmax={qmax:g} 1/A",
        r,
        finite_q_pdf,
    )
    metrics = error_metrics(histogram_pdf, finite_q_pdf, fit_selection)
    series.append({"qmax_inverse_angstrom": qmax, **metrics})
    print(
        f"Qmax={qmax:>4g}  RMS/range={metrics['rms_over_dynamic_range']:.6%}  "
        f"max/range={metrics['max_over_dynamic_range']:.6%}"
    )

expected_path = case_root / "expected/qmax-series.toml"
expected = tomllib.loads(expected_path.read_text())
if len(expected["case"]) != len(series):
    raise AssertionError("expected and computed Qmax series have different lengths")
for observed, limit in zip(series, expected["case"]):
    if observed["qmax_inverse_angstrom"] != limit["qmax_inverse_angstrom"]:
        raise AssertionError("expected and computed Qmax grids differ")
    if observed["rms_over_dynamic_range"] > limit["rms_regression_max"]:
        raise AssertionError(f"Qmax={observed['qmax_inverse_angstrom']:g} RMS regression limit exceeded")
    if observed["max_over_dynamic_range"] > limit["max_regression_max"]:
        raise AssertionError(f"Qmax={observed['qmax_inverse_angstrom']:g} maximum-error limit exceeded")
rms_values = [item["rms_over_dynamic_range"] for item in series]
if not all(left > right for left, right in zip(rms_values, rms_values[1:])):
    raise AssertionError("finite-Q RMS error must decrease strictly as Qmax rises")
if rms_values[-1] > expected["qmax_30_rms_over_dynamic_range_max"]:
    raise AssertionError("Qmax=30 histogram/back-transform agreement regressed")
if rms_values[0] / rms_values[-1] < expected["qmax_12_to_30_rms_ratio_min"]:
    raise AssertionError("the expected low-Q termination effect is too weak")

summary = {
    "benchmark_scope": "transparent reconstruction of RMCProfile7 Figure 11 protocol",
    "native_rmcprofile_curve_comparison": False,
    "cells": args.cells,
    "atom_count": 12 * args.cells**3,
    "box_angstrom": args.cells * args.lattice,
    "lattice_angstrom": args.lattice,
    "uiso_ca_square_angstrom": args.uiso_ca,
    "uiso_f_square_angstrom": args.uiso_f,
    "seed": args.seed,
    "neutron_scattering_lengths_fm": {"Ca": B_CA, "F": B_F},
    "fit_min_angstrom": 1.5,
    "fit_max_angstrom": 15.0,
    "series": series,
}
(results / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
with (results / "summary.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(series[0]))
    writer.writeheader()
    writer.writerows(series)
print(f"wrote curves and metrics to {results}")
print(f"passed frozen regression checks from {expected_path}")
