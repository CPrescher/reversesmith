#!/usr/bin/env python3
"""Prepare rsmith inputs for the frozen five-seed 10 GPa comparison."""

from __future__ import annotations

import argparse
import re
import shutil
import tomllib
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--force", action="store_true")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
protocol = tomllib.loads((case_root / "expected/ten-gpa-comparison.toml").read_text())
pilot_root = case_root / "results/ten-gpa-pilot"
root = case_root / "results/ten-gpa-comparison"
if not (pilot_root / "inputs/fixture-summary.json").is_file():
    raise SystemExit("prepare the 10 GPa pilot first")
if root.exists():
    if not args.force:
        raise SystemExit(f"output exists: {root} (pass --force)")
    shutil.rmtree(root)
(root / "runs").mkdir(parents=True)


def replace_one(text: str, pattern: str, replacement: str):
    result, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise ValueError(f"expected one match for {pattern!r}")
    return result


base = (pilot_root / "runs/rsmith-pace-w3/moves-030000/run.toml").read_text()
base = replace_one(base, r"^weight = [0-9.eE+-]+$", "weight = 30.0")
primary = int(protocol["design"]["primary_seed"])
endpoint = int(protocol["design"]["endpoint_moves"])
for method in ("rsmith-rmc", "rsmith-pace-w30"):
    for seed in protocol["design"]["seeds"]:
        checkpoints = (
            protocol["design"]["primary_seed_checkpoints_moves"]
            if method == "rsmith-pace-w30" and int(seed) == primary
            else [endpoint]
        )
        if method == "rsmith-rmc" and int(seed) == primary:
            continue
        for moves in checkpoints:
            run = root / "runs" / method / f"seed-{int(seed)}" / f"moves-{int(moves):06d}"
            run.mkdir(parents=True)
            config = base
            if method == "rsmith-rmc":
                config = re.sub(r"\n\[ml_potential\].*\Z", "\n", config, flags=re.DOTALL)
            config = replace_one(config, r"^max_moves = \d+$", f"max_moves = {int(moves)}")
            config = replace_one(config, r"^seed = \d+$", f"seed = {int(seed)}")
            config = replace_one(config, r"^print_every = \d+$", f"print_every = {int(moves)}")
            (run / "run.toml").write_text(config)

print(root)
