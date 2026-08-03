#!/usr/bin/env python3
"""Run one rsmith family from the frozen SiO2 multi-seed comparison."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from pathlib import Path


METHODS = ("rsmith-rmc", "rsmith-pedone-w3", "rsmith-pace-w3", "rsmith-epsr")
parser = argparse.ArgumentParser()
parser.add_argument("--method", choices=METHODS, required=True)
parser.add_argument("--binary", type=Path)
parser.add_argument("--only-missing", action="store_true")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
repo_root = case_root.parents[1]
binary = args.binary or repo_root / "target/release/rsmith"
if not binary.is_file():
    raise SystemExit(f"rsmith binary not found: {binary}")
root = case_root / "results/multiseed-comparison"
summary_path = root / "rsmith-run-summary.json"
summary = (
    json.loads(summary_path.read_text())
    if summary_path.is_file()
    else {"program": "rsmith", "threads": 1, "cases": {}}
)
summary.setdefault("invocations", []).append(
    {"method": args.method, "binary": str(binary)}
)
environment = {
    **os.environ,
    "RAYON_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}
for case in sorted(root.glob("target-*_*")):
    case_summary = summary["cases"].setdefault(case.name, {})
    method_summary = case_summary.setdefault(args.method, {})
    for run in sorted((case / args.method).glob("seed-*")):
        if args.only_missing and (run / "refined.xyz").is_file():
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
        if completed.returncode != 0:
            raise RuntimeError(f"rsmith failed in {run}; see driver.log")
        method_summary[run.name] = {
            "status": "completed",
            "wall_seconds": wall,
            "binary": str(binary),
        }
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        print(json.dumps({case.name: {args.method: method_summary[run.name]}}, indent=2))
