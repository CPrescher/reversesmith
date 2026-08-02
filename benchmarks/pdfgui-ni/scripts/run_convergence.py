#!/usr/bin/env python3
"""One-factor-at-a-time convergence study for the PDFgui Ni forward benchmark."""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
import subprocess
import tomllib
from pathlib import Path


def read_curve(path: Path):
    rows = []
    for raw in path.read_text().splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        fields = raw.split()
        try:
            rows.append((float(fields[0]), float(fields[1])))
        except (ValueError, IndexError):
            continue
    if len(rows) < 2:
        raise ValueError(f"{path} contains fewer than two numeric rows")
    return rows


def interpolate(curve, x):
    xs = [point[0] for point in curve]
    hi = bisect.bisect_left(xs, x)
    if hi < len(curve) and curve[hi][0] == x:
        return curve[hi][1]
    if hi == 0 or hi == len(curve):
        raise ValueError(f"{x} is outside the candidate grid")
    x0, y0 = curve[hi - 1]
    x1, y1 = curve[hi]
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0)


def metrics(reference, candidate, fit_min=1.7, fit_max=20.0):
    overlap_max = min(fit_max, candidate[-1][0])
    selected = [point for point in reference if fit_min <= point[0] <= overlap_max]
    residuals = [interpolate(candidate, x) - y for x, y in selected]
    dynamic_range = max(y for _, y in selected) - min(y for _, y in selected)
    rms = math.sqrt(sum(value * value for value in residuals) / len(residuals))
    maximum = max(abs(value) for value in residuals)
    return {
        "n_points": len(residuals),
        "rms_error": rms,
        "max_absolute_error": maximum,
        "rms_error_over_dynamic_range": rms / dynamic_range,
        "max_absolute_error_over_dynamic_range": maximum / dynamic_range,
    }


def average_curves(curves):
    if any(len(curve) != len(curves[0]) for curve in curves):
        raise ValueError("replica grids have different lengths")
    averaged = []
    for points in zip(*curves):
        r = points[0][0]
        if any(abs(point[0] - r) > 1e-12 for point in points[1:]):
            raise ValueError("replica r grids differ")
        averaged.append((r, sum(point[1] for point in points) / len(points)))
    return averaged


parser = argparse.ArgumentParser()
parser.add_argument("--replicas", type=int, default=3)
parser.add_argument("--first-seed", type=int, default=20260802)
parser.add_argument("--binary", type=Path)
args = parser.parse_args()
if args.replicas < 1:
    raise SystemExit("replicas must be positive")

root = Path(__file__).resolve().parents[1]
repo_root = root.parents[1]
binary = args.binary or repo_root / "target/release/rsmith"
reference_path = root / "reference/pdfgui-calculated.dat"
data_path = root / "reference/ni-xray-r-G-sigma.dat"
generator = root / "scripts/generate_ni_supercell.py"
work = root / "results/convergence"
work.mkdir(parents=True, exist_ok=True)
reference = read_curve(reference_path)

cases = [
    {"id": "base", "cells": 12, "rdf_nbins": 2000, "nq": 2700},
    {"id": "box-14", "cells": 14, "rdf_nbins": 2000, "nq": 2700},
    {"id": "box-16", "cells": 16, "rdf_nbins": 2000, "nq": 2700},
    {"id": "dr-0.020", "cells": 12, "rdf_nbins": 1000, "nq": 2700},
    {"id": "dr-0.005", "cells": 12, "rdf_nbins": 4000, "nq": 2700},
    {"id": "dq-0.020", "cells": 12, "rdf_nbins": 2000, "nq": 1350},
    {"id": "dq-0.005", "cells": 12, "rdf_nbins": 2000, "nq": 5400},
]

results = []
for case in cases:
    curves = []
    individual = []
    for replica in range(args.replicas):
        seed = args.first_seed + replica
        structure = work / f"ni-c{case['cells']}-s{seed}.data"
        if not structure.exists():
            subprocess.run(
                [
                    "python3",
                    str(generator),
                    "--cells",
                    str(case["cells"]),
                    "--lattice",
                    "3.52",
                    "--uiso",
                    "0.00126651",
                    "--seed",
                    str(seed),
                    "--output",
                    str(structure),
                ],
                check=True,
            )
        case_dir = work / case["id"] / f"seed-{seed}"
        case_dir.mkdir(parents=True, exist_ok=True)
        config = case_dir / "forward.toml"
        config.write_text(
            f'''[system]
structure = "{structure}"
format = "lammps"
types = {{ "1" = "Ni" }}

[data]
[data.xray_fr]
file = "{data_path}"
sigma_column = 3
fit_min = 1.7
fit_max = 20.0
qmax = 27.0
lorch = false
qdamp = 0.08

[rmc]
max_moves = 0
seed = {seed}

[sq]
qmin = 0.01
qmax = 27.0
nq = {case["nq"]}
lorch = false
rdf_cutoff = 20.0
rdf_nbins = {case["rdf_nbins"]}
'''
        )
        subprocess.run(
            [str(binary), str(config), "--output-dir", str(case_dir), "--quiet"],
            check=True,
        )
        curve = read_curve(case_dir / "start_total_fr.dat")
        curves.append(curve)
        individual.append({"seed": seed, **metrics(reference, curve)})

    ensemble = average_curves(curves)
    ensemble_path = work / case["id"] / "ensemble_mean_fr.dat"
    ensemble_path.write_text(
        f"# mean of {args.replicas} replicas\n"
        + "".join(f"{r:.8f} {value:.12g}\n" for r, value in ensemble)
    )
    result = {
        **case,
        "atom_count": 4 * case["cells"] ** 3,
        "rdf_dr": 20.0 / case["rdf_nbins"],
        "nominal_dq": (27.0 - 0.01) / case["nq"],
        "replicas": args.replicas,
        "individual": individual,
        "ensemble": metrics(reference, ensemble),
    }
    results.append(result)
    print(
        f"{case['id']}: RMS/range={result['ensemble']['rms_error_over_dynamic_range']:.6%}, "
        f"max/range={result['ensemble']['max_absolute_error_over_dynamic_range']:.6%}"
    )

summary = {
    "reference": str(reference_path),
    "fit_min": 1.7,
    "fit_max": 20.0,
    "first_seed": args.first_seed,
    "cases": results,
}

expected_path = root / "expected/convergence.toml"
expected = tomllib.loads(expected_path.read_text())
if len(expected["case"]) != len(results):
    raise AssertionError("expected and computed convergence series have different lengths")
for observed, recorded in zip(results, expected["case"]):
    if observed["id"] != recorded["id"]:
        raise AssertionError("expected and computed convergence case order differs")
    ensemble = observed["ensemble"]
    if ensemble["rms_error_over_dynamic_range"] > expected["rms_over_dynamic_range_max"]:
        raise AssertionError(f"{observed['id']} RMS regression limit exceeded")
    if ensemble["max_absolute_error_over_dynamic_range"] > expected["max_over_dynamic_range_max"]:
        raise AssertionError(f"{observed['id']} maximum-error regression limit exceeded")
(work / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
with (work / "summary.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(
        [
            "case",
            "atoms",
            "rdf_dr_A",
            "nominal_dq_invA",
            "replicas",
            "rms_error",
            "max_absolute_error",
            "rms_over_range",
            "max_over_range",
        ]
    )
    for result in results:
        ensemble = result["ensemble"]
        writer.writerow(
            [
                result["id"],
                result["atom_count"],
                result["rdf_dr"],
                result["nominal_dq"],
                result["replicas"],
                ensemble["rms_error"],
                ensemble["max_absolute_error"],
                ensemble["rms_error_over_dynamic_range"],
                ensemble["max_absolute_error_over_dynamic_range"],
            ]
        )
print(f"wrote {work / 'summary.json'} and {work / 'summary.csv'}")
print(f"passed frozen regression checks from {expected_path}")
