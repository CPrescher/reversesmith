#!/usr/bin/env python3
"""Run the frozen 0-to-10 GPa rsmith pilot."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
import tomllib
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--binary", type=Path)
parser.add_argument("--only-missing", action="store_true")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
repo_root = case_root.parents[1]
protocol = tomllib.loads((case_root / "expected/ten-gpa-pilot.toml").read_text())
binary = (args.binary or repo_root / "target/release/rsmith").resolve()
if not binary.is_file():
    raise SystemExit(f"rsmith binary not found: {binary}")

root = case_root / "results/ten-gpa-pilot"
if not (root / "inputs/fixture-summary.json").is_file():
    raise SystemExit("pilot is not prepared; run prepare_ten_gpa_pilot.py first")
summary_path = root / "run-summary.json"
summary = (
    json.loads(summary_path.read_text())
    if summary_path.is_file()
    else {"status": "hp_sio2_10gpa_pilot_runs", "runs": {}}
)
environment = {
    **os.environ,
    "RAYON_NUM_THREADS": str(protocol["pilot"]["threads_per_run"]),
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}

for method in protocol["pilot"]["methods"]:
    for moves in protocol["pilot"]["checkpoints_moves"]:
        name = f"{method}/moves-{int(moves):06d}"
        run = root / "runs" / name
        refined = run / "refined.xyz"
        if args.only_missing and refined.is_file():
            continue
        started = time.perf_counter()
        completed = subprocess.run(
            [str(binary), str(run / "run.toml"), "--output-dir", str(run), "--quiet"],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env=environment,
            check=False,
        )
        wall = time.perf_counter() - started
        (run / "driver.log").write_text(completed.stdout)
        (run / "wall-seconds.txt").write_text(f"{wall:.9f}\n")
        result = {
            "binary": str(binary),
            "returncode": completed.returncode,
            "status": "completed" if completed.returncode == 0 else "failed",
            "wall_seconds": wall,
        }
        summary["runs"][name] = result
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        print(json.dumps({name: result}, indent=2))
        if completed.returncode != 0:
            raise RuntimeError(f"rsmith failed in {run}; see driver.log")
