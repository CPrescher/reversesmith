#!/usr/bin/env python3
"""Run all prepared rsmith cross-recovery zero gates and smoke refinements."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--binary", type=Path)
parser.add_argument("--zero-only", action="store_true")
parser.add_argument("--gap-only", action="store_true")
parser.add_argument(
    "--no-gap", action="store_true", help="run every prepared arm except GAP/QUIP"
)
args = parser.parse_args()
if args.gap_only and args.no_gap:
    raise SystemExit("--gap-only and --no-gap are mutually exclusive")

case_root = Path(__file__).resolve().parents[1]
repo_root = case_root.parents[1]
binary = args.binary or repo_root / "target/release/rsmith"
if not binary.is_file():
    raise SystemExit(f"rsmith release binary not found: {binary}")
fixture_root = case_root / "results/cross-recovery"
names = (
    ["rsmith-gap-joint"]
    if args.gap_only
    else ["rsmith-hidden-zero-move", "rsmith-cross-zero-move"]
)
if not args.zero_only and not args.gap_only:
    names.extend(
        (
            "rsmith-neutron_only",
            "rsmith-xray_only",
            "rsmith-joint",
            "rsmith-pedone-joint",
            "rsmith-gap-joint",
            "rsmith-epsr-joint",
        )
    )
if args.no_gap:
    names = [name for name in names if name != "rsmith-gap-joint"]

summary_path = fixture_root / "rsmith-run-summary.json"
if summary_path.is_file():
    summary = json.loads(summary_path.read_text())
    summary.setdefault("cases", {})
    summary.setdefault("invocations", [])
else:
    summary = {"program": "rsmith", "threads": 1, "cases": {}, "invocations": []}
summary["invocations"].append({"binary": str(binary), "methods": names})
for case in sorted(fixture_root.glob("target-*_*")):
    case_summary = summary["cases"].setdefault(case.name, {})
    for name in names:
        run_dir = case / name
        config = run_dir / "run.toml"
        if not config.is_file():
            case_summary[name] = {"status": "not_prepared"}
            continue
        started = time.perf_counter()
        completed = subprocess.run(
            [str(binary), str(config), "--output-dir", str(run_dir), "--quiet"],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env={
                **os.environ,
                "RAYON_NUM_THREADS": "1",
                "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
                "VECLIB_MAXIMUM_THREADS": "1",
            },
            check=False,
        )
        wall = time.perf_counter() - started
        (run_dir / "driver.log").write_text(completed.stdout)
        (run_dir / "wall-seconds.txt").write_text(f"{wall:.9f}\n")
        if completed.returncode != 0:
            raise RuntimeError(f"rsmith failed in {run_dir}; see driver.log")
        case_summary[name] = {
            "status": "completed",
            "wall_seconds": wall,
            "binary": str(binary),
        }
    summary["cases"][case.name] = case_summary

summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(json.dumps(summary, indent=2, sort_keys=True))
