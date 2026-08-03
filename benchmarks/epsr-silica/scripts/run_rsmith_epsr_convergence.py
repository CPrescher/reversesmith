#!/usr/bin/env python3
"""Run independent-prefix rsmith EPSR convergence pilots for SiO2."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import time
import tomllib
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--binary", type=Path)
parser.add_argument("--force", action="store_true")
parser.add_argument("--only-missing", action="store_true")
parser.add_argument("--no-hard-constraints", action="store_true")
parser.add_argument("--checkpoint", type=int, action="append")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
repo_root = case_root.parents[1]
binary = args.binary or repo_root / "target/release/rsmith"
if not binary.is_file():
    raise SystemExit(f"rsmith binary not found: {binary}")
protocol = tomllib.loads(
    (case_root / "expected/epsr-convergence-pilot.toml").read_text()
)
fixture_root = case_root / "results/cross-recovery"
pilot_root = case_root / "results/epsr-convergence-pilot"
method = "rsmith-epsr-unconstrained" if args.no_hard_constraints else "rsmith-epsr"
summary_path = pilot_root / f"{method}-run-summary.json"
summary = {
    "program": "rsmith",
    "binary": str(binary),
    "threads": 1,
    "sampling": "independent deterministic prefixes",
    "cases": {},
}
seed = int(protocol["design"]["seed"])
moves_per_iteration = int(protocol["sampling"]["moves_per_refinement"])
environment = {
    **os.environ,
    "RAYON_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}


def replace_one(text: str, pattern: str, replacement: str) -> str:
    rewritten, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise ValueError(f"expected one configuration match for {pattern!r}")
    return rewritten


for case in sorted(fixture_root.glob("target-*_*")):
    source = case / "rsmith-epsr-joint/run.toml"
    if not source.is_file():
        raise SystemExit(f"missing prerequisite: {source}")
    prefix_root = pilot_root / case.name / method / f"seed-{seed}"
    prefix_root.mkdir(parents=True, exist_ok=True)
    case_summary = summary["cases"].setdefault(case.name, {})
    checkpoints = args.checkpoint or protocol["design"]["checkpoints"]
    for refinements in checkpoints:
        refinements = int(refinements)
        run = prefix_root / f"iter-{refinements:03d}"
        if args.only_missing and (run / "refined.xyz").is_file():
            continue
        if run.exists():
            if not args.force:
                raise SystemExit(f"output exists: {run} (pass --force)")
            shutil.rmtree(run)
        run.mkdir(parents=True)
        config = source.read_text()
        if args.no_hard_constraints:
            config, count = re.subn(
                r"\n\[constraints\]\n.*?(?=\n\[potential\])",
                "\n",
                config,
                count=1,
                flags=re.DOTALL,
            )
            if count != 1:
                raise ValueError("expected one [constraints] section")
        config = replace_one(
            config,
            r"^max_moves = \d+$",
            f"max_moves = {refinements * moves_per_iteration}",
        )
        config = replace_one(config, r"^seed = \d+$", f"seed = {seed}")
        config = replace_one(
            config,
            r"^print_every = \d+$",
            f"print_every = {moves_per_iteration}",
        )
        config = replace_one(
            config,
            r"^iterations = \d+$",
            f"iterations = {refinements}",
        )
        config = replace_one(
            config,
            r"^moves_per_iteration = \d+$",
            f"moves_per_iteration = {moves_per_iteration}",
        )
        config_path = run / "run.toml"
        config_path.write_text(config)
        started = time.perf_counter()
        completed = subprocess.run(
            [str(binary), str(config_path), "--output-dir", str(run), "--quiet"],
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
        case_summary[run.name] = {
            "status": "completed",
            "wall_seconds": wall,
            "seed": seed,
            "refinements": refinements,
            "attempted_moves": refinements * moves_per_iteration,
        }
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        print(json.dumps({case.name: {run.name: case_summary[run.name]}}, indent=2))
