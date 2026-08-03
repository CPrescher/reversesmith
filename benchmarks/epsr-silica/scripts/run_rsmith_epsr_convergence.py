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
parser.add_argument("--convergence-ensemble", action="store_true")
parser.add_argument("--control-start-sensitivity", action="store_true")
parser.add_argument("--sensitivity-arm", action="append")
parser.add_argument("--convergence-seed", type=int, action="append")
parser.add_argument("--checkpoint", type=int, action="append")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
repo_root = case_root.parents[1]
binary = args.binary or repo_root / "target/release/rsmith"
if not binary.is_file():
    raise SystemExit(f"rsmith binary not found: {binary}")
if args.convergence_ensemble and args.control_start_sensitivity:
    raise SystemExit("choose only one convergence protocol")
protocol_name = (
    "epsr-control-start-sensitivity.toml"
    if args.control_start_sensitivity
    else (
        "epsr-convergence-ensemble.toml"
        if args.convergence_ensemble
        else "epsr-convergence-pilot.toml"
    )
)
protocol = tomllib.loads((case_root / "expected" / protocol_name).read_text())
fixture_root = case_root / "results/cross-recovery"
convergence_root = case_root / "results" / (
    "epsr-control-start-sensitivity"
    if args.control_start_sensitivity
    else (
        "epsr-convergence-ensemble"
        if args.convergence_ensemble
        else "epsr-convergence-pilot"
    )
)
unconstrained = (
    args.no_hard_constraints
    or args.convergence_ensemble
    or args.control_start_sensitivity
)
method = "rsmith-epsr-unconstrained" if unconstrained else "rsmith-epsr"
configured_seeds = (
    [int(seed) for seed in protocol["design"]["refinement_seeds"]]
    if args.control_start_sensitivity
    else (
        [int(seed) for seed in protocol["design"]["seeds"]]
        if args.convergence_ensemble
        else [int(protocol["design"]["seed"])]
    )
)
selected_seeds = args.convergence_seed or configured_seeds
unknown_seeds = set(selected_seeds) - set(configured_seeds)
if unknown_seeds:
    raise SystemExit(f"seeds not present in {protocol_name}: {sorted(unknown_seeds)}")
moves_per_iteration = int(
    protocol["design"]["moves_per_refinement"]
    if args.convergence_ensemble or args.control_start_sensitivity
    else protocol["sampling"]["moves_per_refinement"]
)
checkpoints = args.checkpoint or (
    protocol["methods"]["rsmith_epsr_unconstrained"]["checkpoints"]
    if args.convergence_ensemble or args.control_start_sensitivity
    else protocol["design"]["checkpoints"]
)
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


if args.control_start_sensitivity:
    inputs_root = convergence_root / "inputs"
    if not (inputs_root / "input-manifest.json").is_file():
        raise SystemExit("run prepare_epsr_control_start_sensitivity.py first")
    arms = {arm["name"]: arm for arm in protocol["arms"]}
    selected_arms = args.sensitivity_arm or list(arms)
    unknown_arms = set(selected_arms) - set(arms)
    if unknown_arms:
        raise SystemExit(f"unknown sensitivity arms: {sorted(unknown_arms)}")
    for arm_name in selected_arms:
        arm = arms[arm_name]
        reference_label = str(float(arm["reference_scale"])).replace(".", "p")
        reference_root = (
            inputs_root / "reference-potentials" / f"scale-{reference_label}"
        )
        for seed in selected_seeds:
            summary_path = (
                convergence_root
                / "run-summaries"
                / f"{arm_name}-rsmith-epsr-unconstrained-seed-{seed}.json"
            )
            summary = {
                "program": "rsmith",
                "binary": str(binary),
                "arm": arm_name,
                "feedback": float(arm["feedback"]),
                "reference_scale": float(arm["reference_scale"]),
                "start": arm["start"],
                "threads": 1,
                "seed": seed,
                "cases": {},
            }
            for case in sorted(fixture_root.glob("target-*_*")):
                source = case / "rsmith-epsr-joint/run.toml"
                if not source.is_file():
                    raise SystemExit(f"missing prerequisite: {source}")
                structure = (
                    case / "cross-start.data"
                    if arm["start"] == "original"
                    else inputs_root
                    / "starts"
                    / case.name
                    / arm["start"]
                    / "start-charge.data"
                )
                if not structure.is_file():
                    raise SystemExit(f"missing sensitivity start: {structure}")
                prefix_root = (
                    convergence_root
                    / case.name
                    / arm_name
                    / "rsmith-epsr-unconstrained"
                    / f"seed-{seed}"
                )
                prefix_root.mkdir(parents=True, exist_ok=True)
                case_summary = summary["cases"].setdefault(case.name, {})
                for refinements in checkpoints:
                    refinements = int(refinements)
                    run = prefix_root / f"iter-{refinements:03d}"
                    if args.only_missing:
                        if (run / "refined.xyz").is_file():
                            continue
                        if run.exists():
                            shutil.rmtree(run)
                    elif run.exists():
                        if not args.force:
                            raise SystemExit(f"output exists: {run} (pass --force)")
                        shutil.rmtree(run)
                    run.mkdir(parents=True)
                    config = source.read_text()
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
                        r'^structure = ".*"$',
                        f'structure = "{structure}"',
                    )
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
                        r"^feedback = [0-9.eE+-]+$",
                        f"feedback = {float(arm['feedback']):.8f}",
                    )
                    config = replace_one(
                        config,
                        r"^moves_per_iteration = \d+$",
                        f"moves_per_iteration = {moves_per_iteration}",
                    )
                    for pair in ("Si-Si", "Si-O", "O-O"):
                        path = reference_root / f"epsr-reference-{pair}.dat"
                        config = replace_one(
                            config,
                            rf'^file = ".*epsr-reference-{pair}\.dat"$',
                            f'file = "{path}"',
                        )
                    config_path = run / "run.toml"
                    config_path.write_text(config)
                    started = time.perf_counter()
                    completed = subprocess.run(
                        [
                            str(binary),
                            str(config_path),
                            "--output-dir",
                            str(run),
                            "--quiet",
                        ],
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
                    summary_path.parent.mkdir(parents=True, exist_ok=True)
                    summary_path.write_text(
                        json.dumps(summary, indent=2, sort_keys=True) + "\n"
                    )
                    print(
                        json.dumps(
                            {case.name: {run.name: case_summary[run.name]}},
                            indent=2,
                        )
                    )
    raise SystemExit(0)


for seed in selected_seeds:
    summary_path = convergence_root / (
        f"{method}-seed-{seed}-run-summary.json"
        if args.convergence_ensemble
        else f"{method}-run-summary.json"
    )
    summary = {
        "program": "rsmith",
        "binary": str(binary),
        "threads": 1,
        "sampling": "independent deterministic prefixes",
        "seed": seed,
        "cases": {},
    }
    for case in sorted(fixture_root.glob("target-*_*")):
        source = case / "rsmith-epsr-joint/run.toml"
        if not source.is_file():
            raise SystemExit(f"missing prerequisite: {source}")
        prefix_root = convergence_root / case.name / method / f"seed-{seed}"
        prefix_root.mkdir(parents=True, exist_ok=True)
        case_summary = summary["cases"].setdefault(case.name, {})
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
            if unconstrained:
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
            summary_path.write_text(
                json.dumps(summary, indent=2, sort_keys=True) + "\n"
            )
            print(json.dumps({case.name: {run.name: case_summary[run.name]}}, indent=2))
