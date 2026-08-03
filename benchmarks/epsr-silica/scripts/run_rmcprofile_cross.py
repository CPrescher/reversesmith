#!/usr/bin/env python3
"""Run native RMCProfile cross-recovery zero-move gates and smoke fits."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import time
from pathlib import Path


def read_iq(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            q, value, *_ = (float(field) for field in line.split())
            rows.append((q, value))
    return rows


def write_rmc_curve(path: Path, title: str, rows, transform):
    path.write_text(
        f"{len(rows):12d}\n{title}\n"
        + "".join(f"{q:20.12g} {transform(q, value):20.12g}\n" for q, value in rows)
    )


def control_text(density: float, moves: int, include_neutron: bool, include_xray: bool):
    blocks = []
    if include_neutron:
        blocks.append(
            """NEUTRON_RECIPROCAL_SPACE_DATA :: 1
  > FILENAME :: target-neutron.sq
  > DATA_TYPE :: S(Q) normalized
  > FIT_TYPE :: S(Q) normalized
  > START_POINT :: 1
  > END_POINT :: 490
  > CONSTANT_OFFSET :: 0.0
  > WEIGHT :: 1.0
  > NO_FITTED_OFFSET
  > NO_FITTED_SCALE
"""
        )
    if include_xray:
        blocks.append(
            """XRAY_RECIPROCAL_SPACE_DATA ::
  > FILENAME :: target-xray.fq
  > DATA_TYPE :: F(Q)
  > FIT_TYPE :: F(Q)
  > NO_FITTED_OFFSET
  > NO_FITTED_SCALE
  > NORMALIZATION_TYPE :: <f>^2
  > RECIPROCAL_SPACE_FIT :: 1 490 1
  > RECIPROCAL_SPACE_PARAMETERS :: 1 490 1.0
"""
        )
    # SAVE_PERIOD=0 makes RMCProfile save after nearly every accepted move.
    # A 0.05-minute interval captures recovered coordinates without dominating
    # this short smoke run with thousands of writes.
    return f"""TITLE :: SiO2 symmetric GAP/Pedone cross recovery
MATERIAL :: SiO2
PHASE :: glass
TEMPERATURE :: 300K

NUMBER_DENSITY :: {density:.15g} Angstrom^(-3)
MINIMUM_DISTANCES :: 2.0 1.35 2.0 Angstrom
MAXIMUM_MOVES :: 0.10 0.10 Angstrom
R_SPACING :: 0.0200 Angstrom
PRINT_PERIOD :: {max(1, moves)}
SAVE_PERIOD :: 0.05 MINUTES
TIME_LIMIT :: {0.0 if moves == 0 else 1000.0} MINUTES
ITERATION_LIMIT :: {moves}

INPUT_CONFIGURATION_FORMAT :: rmc6f
SAVE_CONFIGURATION_FORMAT :: rmc6f
ATOMS :: Si O
IGNORE_HISTORY_FILE ::

FLAGS ::
  > NO_MOVEOUT
  > NO_RESOLUTION_CONVOLUTION

{"".join(blocks)}
END ::
"""


def read_rmcprofile_calculated(path: Path):
    rows = []
    with path.open(newline="") as stream:
        for row in csv.reader(stream):
            try:
                rows.append((float(row[0]), float(row[1])))
            except (ValueError, IndexError):
                continue
    return rows


def curve_rms(first, second):
    if len(first) != len(second):
        raise ValueError("RMCProfile and independent target grids differ")
    residual = [a[1] - b[1] for a, b in zip(first, second)]
    return math.sqrt(sum(value * value for value in residual) / len(residual))


def prepare_run(
    case: Path,
    name: str,
    structure: Path,
    moves: int,
    fit_mode: str,
    force: bool,
    native_targets: bool = False,
):
    run_dir = case / name
    if run_dir.exists():
        if not force:
            raise FileExistsError(f"output exists: {run_dir} (pass --force)")
        shutil.rmtree(run_dir)
    run_dir.mkdir()
    shutil.copy2(structure, run_dir / "cross.rmc6f")
    if native_targets:
        shutil.copy2(
            case / "rmcprofile-native-target-neutron.sq", run_dir / "target-neutron.sq"
        )
        shutil.copy2(
            case / "rmcprofile-native-target-xray.fq", run_dir / "target-xray.fq"
        )
    else:
        neutron = read_iq(case / "target-neutron-iq.dat")
        xray = read_iq(case / "target-xray-iq.dat")
        write_rmc_curve(
            run_dir / "target-neutron.sq",
            "provisional independent SiO2 neutron S(Q)",
            neutron,
            lambda _q, iq: iq + 1.0,
        )
        write_rmc_curve(
            run_dir / "target-xray.fq",
            "provisional independent SiO2 X-ray F(Q)=Q[S(Q)-1]",
            xray,
            lambda q, iq: q * iq,
        )
    box_line = next(
        line
        for line in structure.read_text().splitlines()
        if line.startswith("Cell (Ang/deg):")
    )
    box = [float(value) for value in box_line.split(":", 1)[1].split()[:3]]
    density = 3000.0 / (box[0] * box[1] * box[2])
    (run_dir / "cross.dat").write_text(
        control_text(
            density,
            moves,
            fit_mode in {"neutron_only", "joint"},
            fit_mode in {"xray_only", "joint"},
        )
    )
    return run_dir


def run(executable: Path, run_dir: Path):
    started = time.perf_counter()
    completed = subprocess.run(
        [str(executable), "cross"],
        cwd=run_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env={**os.environ, "PATH": f"{executable.parent}:{os.environ.get('PATH', '')}"},
        check=False,
    )
    wall = time.perf_counter() - started
    (run_dir / "driver.log").write_text(completed.stdout)
    (run_dir / "wall-seconds.txt").write_text(f"{wall:.9f}\n")
    if completed.returncode != 0:
        raise RuntimeError(f"RMCProfile failed in {run_dir}; see driver.log")
    return wall


parser = argparse.ArgumentParser()
parser.add_argument(
    "--executable",
    type=Path,
    default=Path("/Applications/RMCProfile.app/Contents/MacOS/exe/rmcprofile"),
)
parser.add_argument("--moves", type=int, default=6000)
parser.add_argument("--force", action="store_true")
parser.add_argument("--zero-only", action="store_true")
args = parser.parse_args()
if args.moves <= 0:
    raise SystemExit("--moves must be positive")
if not args.executable.is_file():
    raise SystemExit(f"RMCProfile executable not found: {args.executable}")

case_root = Path(__file__).resolve().parents[1]
fixture_root = case_root / "results/cross-recovery"
summary_path = fixture_root / "rmcprofile-run-summary.json"
if args.zero_only and summary_path.is_file():
    summary = json.loads(summary_path.read_text())
    summary.update(
        {
            "program": "RMCProfile",
            "executable": str(args.executable),
            "moves": args.moves,
        }
    )
else:
    summary = {
        "program": "RMCProfile",
        "executable": str(args.executable),
        "moves": args.moves,
        "cases": {},
    }
for case in sorted(fixture_root.glob("target-*_*")):
    case_summary = summary["cases"].get(case.name, {})
    zero = prepare_run(
        case,
        "rmcprofile-hidden-zero-move",
        case / "hidden-target.rmc6f",
        0,
        "joint",
        args.force,
        native_targets=False,
    )
    case_summary["native_target_generation_wall_seconds"] = run(args.executable, zero)
    native_neutron = read_rmcprofile_calculated(zero / "cross_SQ1.csv")
    native_xray = read_rmcprofile_calculated(zero / "cross_FQ1.csv")
    independent_neutron = [
        (q, iq + 1.0) for q, iq in read_iq(case / "target-neutron-iq.dat")
    ]
    independent_xray = [(q, q * iq) for q, iq in read_iq(case / "target-xray-iq.dat")]
    case_summary["native_vs_independent_neutron_rms"] = curve_rms(
        native_neutron, independent_neutron
    )
    case_summary["native_vs_independent_xray_fq_rms"] = curve_rms(
        native_xray, independent_xray
    )
    write_rmc_curve(
        case / "rmcprofile-native-target-neutron.sq",
        "RMCProfile-native hidden-coordinate neutron S(Q)",
        native_neutron,
        lambda _q, value: value,
    )
    write_rmc_curve(
        case / "rmcprofile-native-target-xray.fq",
        "RMCProfile-native hidden-coordinate X-ray F(Q)",
        native_xray,
        lambda _q, value: value,
    )
    zero = prepare_run(
        case,
        "rmcprofile-hidden-zero-move",
        case / "hidden-target.rmc6f",
        0,
        "joint",
        True,
        native_targets=True,
    )
    case_summary["hidden_zero_move_wall_seconds"] = run(args.executable, zero)
    checked_neutron = read_rmcprofile_calculated(zero / "cross_SQ1.csv")
    checked_xray = read_rmcprofile_calculated(zero / "cross_FQ1.csv")
    case_summary["hidden_zero_move_neutron_rms"] = curve_rms(
        checked_neutron, native_neutron
    )
    case_summary["hidden_zero_move_xray_fq_rms"] = curve_rms(checked_xray, native_xray)
    cross_zero = prepare_run(
        case,
        "rmcprofile-cross-zero-move",
        case / "cross-start.rmc6f",
        0,
        "joint",
        args.force,
        native_targets=True,
    )
    case_summary["cross_zero_move_wall_seconds"] = run(args.executable, cross_zero)
    if not args.zero_only:
        for fit_mode in ("neutron_only", "xray_only", "joint"):
            run_dir = prepare_run(
                case,
                f"rmcprofile-{fit_mode}",
                case / "cross-start.rmc6f",
                args.moves,
                fit_mode,
                args.force,
                native_targets=True,
            )
            case_summary[f"{fit_mode}_wall_seconds"] = run(args.executable, run_dir)
    summary["cases"][case.name] = case_summary

summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(json.dumps(summary, indent=2, sort_keys=True))
