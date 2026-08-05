#!/usr/bin/env python3
"""Run the frozen S(Q)-only, g(r)-only, and joint PACE arms."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--binary", type=Path)
parser.add_argument("--only-missing", action="store_true")
args = parser.parse_args()
case_root = Path(__file__).resolve().parents[1]
binary = (args.binary or case_root.parents[1] / "target/release/rsmith").resolve()
root = case_root / "results/sq-gr-domain"
if not binary.is_file() or not (root / "input-manifest.json").is_file():
    raise SystemExit("missing release binary or prepared domain comparison")
summary_path = root / "run-summary.json"
summary = (
    json.loads(summary_path.read_text())
    if summary_path.is_file()
    else {"status": "hp_sio2_sq_gr_domain_runs", "runs": {}}
)
environment = {
    **os.environ,
    "RAYON_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}
for run in sorted((root / "runs").iterdir()):
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
    summary["runs"][run.name] = {
        "returncode": completed.returncode,
        "status": "completed" if completed.returncode == 0 else "failed",
        "wall_seconds": wall,
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({run.name: summary["runs"][run.name]}, indent=2), flush=True)
    if completed.returncode != 0:
        raise RuntimeError(f"rsmith failed in {run}; see driver.log")
