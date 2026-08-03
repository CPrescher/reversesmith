#!/usr/bin/env python3
"""Run reproducible, exactly saved RMCProfile SiO2 ensemble trajectories."""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import time
import tomllib
from pathlib import Path


def numeric_curve(path: Path) -> list[tuple[float, float]]:
    rows = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) < 2:
            continue
        try:
            rows.append((float(fields[0]), float(fields[1])))
        except ValueError:
            continue
    return rows


def metadata_count(text: str, label: str) -> int:
    match = re.search(rf"Number of moves {label}:\s+(\d+)", text)
    if match is None:
        raise ValueError(f"missing RMCProfile {label} counter")
    return int(match.group(1))


def reset_metadata(path: Path) -> None:
    text = path.read_text()
    for label in ("generated", "tried", "accepted"):
        text, count = re.subn(
            rf"(?m)^(Number of moves {label}:)\s+\d+$",
            r"\1           0",
            text,
        )
        if count != 1:
            raise ValueError(f"expected one {label} counter in {path}")
    text, count = re.subn(
        r"(?m)^(Accumulated time \(s\) in running loop:)\s+\S+$",
        r"\1 0.0",
        text,
    )
    if count != 1:
        raise ValueError(f"expected one accumulated-time counter in {path}")
    path.write_text(text)


def control_text(density: float, moves: int, save_period: float, xray_q: list[float]) -> str:
    sigma_fq = 0.02 * math.sqrt(sum(q * q for q in xray_q) / len(xray_q))
    return f"""TITLE :: SiO2 symmetric cross-program multi-seed recovery
MATERIAL :: SiO2
PHASE :: glass
TEMPERATURE :: 300K

NUMBER_DENSITY :: {density:.15g} Angstrom^(-3)
MINIMUM_DISTANCES :: 2.0 1.35 2.0 Angstrom
MAXIMUM_MOVES :: 0.10 0.10 Angstrom
R_SPACING :: 0.0200 Angstrom
PRINT_PERIOD :: {max(1, moves)}
SAVE_PERIOD :: {save_period:g} MINUTES
TIME_LIMIT :: 1000.0 MINUTES
ITERATION_LIMIT :: {moves}

INPUT_CONFIGURATION_FORMAT :: rmc6f
SAVE_CONFIGURATION_FORMAT :: rmc6f
ATOMS :: Si O
IGNORE_HISTORY_FILE ::

FLAGS ::
  > NO_MOVEOUT
  > NO_RESOLUTION_CONVOLUTION

NEUTRON_RECIPROCAL_SPACE_DATA :: 1
  > FILENAME :: target-neutron.sq
  > DATA_TYPE :: S(Q) normalized
  > FIT_TYPE :: S(Q) normalized
  > START_POINT :: 1
  > END_POINT :: 490
  > CONSTANT_OFFSET :: 0.0
  > WEIGHT :: 0.02
  > NO_FITTED_OFFSET
  > NO_FITTED_SCALE

XRAY_RECIPROCAL_SPACE_DATA ::
  > FILENAME :: target-xray.fq
  > DATA_TYPE :: F(Q)
  > FIT_TYPE :: F(Q)
  > NO_FITTED_OFFSET
  > NO_FITTED_SCALE
  > NORMALIZATION_TYPE :: <f>^2
  > RECIPROCAL_SPACE_FIT :: 1 490 1
  > RECIPROCAL_SPACE_PARAMETERS :: 1 490 {sigma_fq:.12g}

END ::
"""


def run(executable: Path, run_dir: Path, seed: int) -> tuple[float, str]:
    started = time.perf_counter()
    completed = subprocess.run(
        [str(executable), "cross", str(seed)],
        cwd=run_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env={**os.environ, "OMP_NUM_THREADS": "1"},
        check=False,
    )
    wall = time.perf_counter() - started
    if completed.returncode != 0:
        (run_dir / "failed-driver.log").write_text(completed.stdout)
        raise RuntimeError(f"RMCProfile failed in {run_dir}; see failed-driver.log")
    return wall, completed.stdout


def completed_count(log: str, label: str) -> int:
    values = re.findall(rf"Moves {label}\s+::\s+(\d+)", log)
    if not values:
        raise ValueError(f"missing completed {label} count")
    return int(values[-1])


parser = argparse.ArgumentParser()
parser.add_argument(
    "--executable",
    type=Path,
    default=Path("/Applications/RMCProfile.app/Contents/MacOS/exe/rmcprofile.x"),
)
parser.add_argument("--seed", type=int, action="append")
parser.add_argument("--force", action="store_true")
parser.add_argument("--only-missing", action="store_true")
args = parser.parse_args()
if not args.executable.is_file():
    raise SystemExit(f"RMCProfile serial executable not found: {args.executable}")

case_root = Path(__file__).resolve().parents[1]
protocol = tomllib.loads(
    (case_root / "expected/multiseed-comparison.toml").read_text()
)
seeds = args.seed or [int(value) for value in protocol["budget"]["seeds"]]
moves = int(protocol["budget"]["moves"])
cross_root = case_root / "results/cross-recovery"
root = case_root / "results/multiseed-comparison"
summary_path = root / "rmcprofile-run-summary.json"
summary = (
    json.loads(summary_path.read_text())
    if summary_path.is_file()
    else {
        "program": "RMCProfile 6.7.9.5 serial",
        "threads": 1,
        "executable": str(args.executable),
        "cases": {},
    }
)

for case in sorted(cross_root.glob("target-*_*")):
    box_line = next(
        line
        for line in (case / "cross-start.rmc6f").read_text().splitlines()
        if line.startswith("Cell (Ang/deg):")
    )
    box = [float(value) for value in box_line.split(":", 1)[1].split()[:3]]
    density = 3000.0 / math.prod(box)
    xray_q = [q for q, _ in numeric_curve(case / "rmcprofile-native-target-xray.fq")]
    if len(xray_q) != 490:
        raise ValueError(f"expected 490 RMCProfile X-ray points, got {len(xray_q)}")
    method_root = root / case.name / "rmcprofile"
    method_root.mkdir(parents=True, exist_ok=True)
    case_summary = summary["cases"].setdefault(case.name, {})
    for seed in seeds:
        run_dir = method_root / f"seed-{seed}"
        if args.only_missing and (run_dir / "move-audit.json").is_file():
            continue
        if run_dir.exists():
            if not args.force:
                raise SystemExit(f"output exists: {run_dir} (pass --force to replace it)")
            shutil.rmtree(run_dir)
        run_dir.mkdir()
        shutil.copy2(case / "cross-start.rmc6f", run_dir / "cross.rmc6f")
        shutil.copy2(
            case / "rmcprofile-native-target-neutron.sq",
            run_dir / "target-neutron.sq",
        )
        shutil.copy2(
            case / "rmcprofile-native-target-xray.fq",
            run_dir / "target-xray.fq",
        )
        (run_dir / "cross.dat").write_text(control_text(density, moves, 0.05, xray_q))

        stage1_wall, stage1_log = run(args.executable, run_dir, seed)
        if completed_count(stage1_log, "generated") != moves:
            raise RuntimeError(f"RMCProfile stage 1 did not complete {moves} moves")
        checkpoint_text = (run_dir / "cross.rmc6f").read_text()
        checkpoint_moves = metadata_count(checkpoint_text, "generated")
        checkpoint_accepted = metadata_count(checkpoint_text, "accepted")
        if not 0 < checkpoint_moves <= moves:
            raise RuntimeError(f"invalid checkpoint move count {checkpoint_moves}")
        shutil.copy2(run_dir / "cross.rmc6f", run_dir / "stage1-checkpoint.rmc6f")
        (run_dir / "stage1-driver.log").write_text(stage1_log)
        (run_dir / "sampling-wall-seconds.txt").write_text(f"{stage1_wall:.9f}\n")

        tail_moves = moves - checkpoint_moves
        tail_wall = 0.0
        tail_log = ""
        tail_accepted = 0
        if tail_moves:
            reset_metadata(run_dir / "cross.rmc6f")
            history = run_dir / "cross.his6f"
            if history.exists():
                history.unlink()
            (run_dir / "cross.dat").write_text(
                control_text(density, tail_moves, 0.0, xray_q)
            )
            tail_seed = seed + 1_000_000 + checkpoint_moves
            tail_wall, tail_log = run(args.executable, run_dir, tail_seed)
            tail_generated = metadata_count(
                (run_dir / "cross.rmc6f").read_text(), "generated"
            )
            tail_accepted = metadata_count(
                (run_dir / "cross.rmc6f").read_text(), "accepted"
            )
            if tail_generated != tail_moves:
                raise RuntimeError(
                    f"RMCProfile exact tail saved {tail_generated}, expected {tail_moves}"
                )
            (run_dir / "tail-driver.log").write_text(tail_log)
        audit = {
            "seed": seed,
            "checkpoint_moves": checkpoint_moves,
            "checkpoint_accepted": checkpoint_accepted,
            "tail_moves": tail_moves,
            "tail_accepted": tail_accepted,
            "total_moves": checkpoint_moves + tail_moves,
            "total_accepted": checkpoint_accepted + tail_accepted,
            "sampling_wall_seconds": stage1_wall,
            "exact_tail_wall_seconds": tail_wall,
            "total_wall_seconds": stage1_wall + tail_wall,
            "neutron_sigma_s_q": 0.02,
            "xray_sigma_fq": "0.02*sqrt(mean(q^2)) as one stable total-information section",
        }
        if audit["total_moves"] != moves:
            raise RuntimeError(f"RMCProfile move audit failed: {audit}")
        (run_dir / "move-audit.json").write_text(
            json.dumps(audit, indent=2, sort_keys=True) + "\n"
        )
        (run_dir / "driver.log").write_text(
            "=== uninterrupted 6000-move sampling stage ===\n"
            + stage1_log
            + "\n=== exact checkpoint-completion tail ===\n"
            + tail_log
        )
        (run_dir / "wall-seconds.txt").write_text(
            f"{stage1_wall + tail_wall:.9f}\n"
        )
        case_summary[f"seed-{seed}"] = {"status": "completed", **audit}
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        print(json.dumps({case.name: {f"seed-{seed}": audit}}, indent=2))
