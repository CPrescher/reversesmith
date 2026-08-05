"""Native RMCProfile adapter for the high-pressure SiO2 series."""

from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import subprocess
import time
from pathlib import Path


def numeric_curve(path: Path):
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


def write_curve(path: Path, title: str, rows):
    path.write_text(
        f"{len(rows):12d}\n{title}\n"
        + "".join(f"{q:20.12g} {value:20.12g}\n" for q, value in rows)
    )


def read_calculated(path: Path):
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
        raise ValueError("RMCProfile and target grids differ")
    residual = [a[1] - b[1] for a, b in zip(first, second)]
    return math.sqrt(sum(value * value for value in residual) / len(residual))


def control_text(density, moves, save_period, xray_q, minimum_distances):
    sigma_fq = 0.02 * math.sqrt(sum(q * q for q in xray_q) / len(xray_q))
    minima = " ".join(f"{value:.8g}" for value in minimum_distances)
    return f"""TITLE :: SiO2 0-to-70 GPa pressure-step recovery
MATERIAL :: SiO2
PHASE :: glass
TEMPERATURE :: 300K

NUMBER_DENSITY :: {density:.15g} Angstrom^(-3)
MINIMUM_DISTANCES :: {minima} Angstrom
MAXIMUM_MOVES :: 0.10 0.10 Angstrom
R_SPACING :: 0.0200 Angstrom
PRINT_PERIOD :: {max(1, moves)}
SAVE_PERIOD :: {save_period:g} MINUTES
TIME_LIMIT :: {0.0 if moves == 0 else 1000.0} MINUTES
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


def prepare_run(
    run: Path,
    structure: Path,
    targets,
    moves: int,
    save_period: float,
    minimum_distances,
):
    if run.exists():
        shutil.rmtree(run)
    run.mkdir(parents=True)
    shutil.copy2(structure, run / "cross.rmc6f")
    write_curve(run / "target-neutron.sq", "RMCProfile neutron S(Q)", targets[0])
    write_curve(run / "target-xray.fq", "RMCProfile X-ray F(Q)", targets[1])
    box_line = next(
        line
        for line in structure.read_text().splitlines()
        if line.startswith("Cell (Ang/deg):")
    )
    box = [float(value) for value in box_line.split(":", 1)[1].split()[:3]]
    density = 3000.0 / math.prod(box)
    xray_q = [q for q, _ in targets[1]]
    (run / "cross.dat").write_text(
        control_text(density, moves, save_period, xray_q, minimum_distances)
    )


def run(executable: Path, run_dir: Path, seed: int):
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


def metadata_count(text: str, label: str):
    match = re.search(rf"Number of moves {label}:\s+(\d+)", text)
    if match is None:
        raise ValueError(f"missing RMCProfile {label} counter")
    return int(match.group(1))


def completed_count(log: str, label: str):
    values = re.findall(rf"Moves {label}\s+::\s+(\d+)", log)
    if not values:
        raise ValueError(f"missing completed {label} count")
    return int(values[-1])


def reset_metadata(path: Path):
    text = path.read_text()
    for label in ("generated", "tried", "accepted"):
        text, count = re.subn(
            rf"(?m)^(Number of moves {label}:)\s+\d+$", r"\1           0", text
        )
        if count != 1:
            raise ValueError(f"expected one {label} counter in {path}")
    text, _ = re.subn(
        r"(?m)^(Accumulated time \(s\) in running loop:)\s+\S+$", r"\1 0.0", text
    )
    path.write_text(text)


def run_exact(
    executable: Path,
    run_dir: Path,
    seed: int,
    moves: int,
    targets,
    minimum_distances,
):
    stage_wall, stage_log = run(executable, run_dir, seed)
    if completed_count(stage_log, "generated") != moves:
        raise RuntimeError(f"RMCProfile did not complete {moves} moves")
    checkpoint_text = (run_dir / "cross.rmc6f").read_text()
    checkpoint_moves = metadata_count(checkpoint_text, "generated")
    checkpoint_accepted = metadata_count(checkpoint_text, "accepted")
    if not 0 < checkpoint_moves <= moves:
        raise RuntimeError(f"invalid checkpoint move count {checkpoint_moves}")
    shutil.copy2(run_dir / "cross.rmc6f", run_dir / "stage1-checkpoint.rmc6f")
    (run_dir / "stage1-driver.log").write_text(stage_log)
    tail_moves, tail_wall, tail_accepted = moves - checkpoint_moves, 0.0, 0
    tail_log = ""
    if tail_moves:
        reset_metadata(run_dir / "cross.rmc6f")
        history = run_dir / "cross.his6f"
        if history.exists():
            history.unlink()
        box_line = next(
            line
            for line in (run_dir / "cross.rmc6f").read_text().splitlines()
            if line.startswith("Cell (Ang/deg):")
        )
        box = [float(value) for value in box_line.split(":", 1)[1].split()[:3]]
        density = 3000.0 / math.prod(box)
        (run_dir / "cross.dat").write_text(
            control_text(
                density,
                tail_moves,
                0.0,
                [q for q, _ in targets[1]],
                minimum_distances,
            )
        )
        tail_wall, tail_log = run(
            executable, run_dir, seed + 1_000_000 + checkpoint_moves
        )
        final_text = (run_dir / "cross.rmc6f").read_text()
        if metadata_count(final_text, "generated") != tail_moves:
            raise RuntimeError("RMCProfile exact completion tail failed")
        tail_accepted = metadata_count(final_text, "accepted")
        (run_dir / "tail-driver.log").write_text(tail_log)
    audit = {
        "checkpoint_accepted": checkpoint_accepted,
        "checkpoint_moves": checkpoint_moves,
        "seed": seed,
        "tail_accepted": tail_accepted,
        "tail_moves": tail_moves,
        "total_accepted": checkpoint_accepted + tail_accepted,
        "total_moves": checkpoint_moves + tail_moves,
        "wall_seconds": stage_wall + tail_wall,
    }
    (run_dir / "move-audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n"
    )
    (run_dir / "driver.log").write_text(stage_log + "\n" + tail_log)
    (run_dir / "wall-seconds.txt").write_text(f"{audit['wall_seconds']:.9f}\n")
    return audit
