#!/usr/bin/env python3
"""Dispatch independent sensitivity arm/seed jobs with bounded concurrency."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tomllib
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--workers", type=int, default=8)
parser.add_argument(
    "--method", choices=("both", "native", "rsmith"), default="both"
)
parser.add_argument("--force", action="store_true")
args = parser.parse_args()
if args.workers < 1:
    raise SystemExit("--workers must be positive")

case_root = Path(__file__).resolve().parents[1]
protocol = tomllib.loads(
    (case_root / "expected/epsr-control-start-sensitivity.toml").read_text()
)
script_root = case_root / "scripts"
result_root = case_root / "results/epsr-control-start-sensitivity"
log_root = result_root / "dispatch-logs"
log_root.mkdir(parents=True, exist_ok=True)

methods = (
    ("native", "rsmith") if args.method == "both" else (args.method,)
)
jobs = [
    (method, arm["name"], int(seed))
    for arm in protocol["arms"]
    for seed in protocol["design"]["refinement_seeds"]
    for method in methods
]


def command(method: str, arm: str, seed: int):
    runner = (
        script_root / "run_native_epsr_cross.py"
        if method == "native"
        else script_root / "run_rsmith_epsr_convergence.py"
    )
    result = [
        sys.executable,
        str(runner),
        "--control-start-sensitivity",
        "--sensitivity-arm",
        arm,
        "--convergence-seed",
        str(seed),
    ]
    result.append("--force" if args.force else "--only-missing")
    return result


def run_job(job):
    method, arm, seed = job
    completed = subprocess.run(
        command(method, arm, seed),
        cwd=case_root.parents[1],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log = log_root / f"{arm}-{method}-seed-{seed}.log"
    log.write_text(completed.stdout)
    return {
        "method": method,
        "arm": arm,
        "seed": seed,
        "returncode": completed.returncode,
        "log": str(log),
    }


failures = []
with ThreadPoolExecutor(max_workers=args.workers) as executor:
    futures = {executor.submit(run_job, job): job for job in jobs}
    for future in as_completed(futures):
        record = future.result()
        print(json.dumps(record), flush=True)
        if record["returncode"] != 0:
            failures.append(record)

if failures:
    raise SystemExit(f"{len(failures)} sensitivity jobs failed; inspect dispatch logs")
print(f"completed {len(jobs)} sensitivity jobs with {args.workers} workers")
