#!/usr/bin/env python3
"""Prepare the frozen rsmith arms of the SiO2 multi-seed comparison."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import tomllib
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--force", action="store_true")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
protocol = tomllib.loads(
    (case_root / "expected/multiseed-comparison.toml").read_text()
)
cross_root = case_root / "results/cross-recovery"
production_root = case_root / "results/hrmc-production-unscreened"
output = case_root / "results/multiseed-comparison"
if output.exists():
    if not args.force:
        raise SystemExit(f"output exists: {output} (pass --force to replace it)")
    shutil.rmtree(output)
output.mkdir(parents=True)


def rewrite(text: str, seed: int, moves: int) -> str:
    for pattern, replacement in (
        (r"(?m)^max_moves = \d+$", f"max_moves = {moves}"),
        (r"(?m)^seed = \d+$", f"seed = {seed}"),
        (r"(?m)^print_every = \d+$", f"print_every = {moves}"),
    ):
        text, count = re.subn(pattern, replacement, text)
        if count != 1:
            raise ValueError(f"expected one match for {pattern}")
    text = re.sub(r"(?m)^moves_per_iteration = \d+$", f"moves_per_iteration = {moves}", text)
    return text


manifest = {"status": "prepared_rsmith_arms", "cases": {}}
for case in sorted(cross_root.glob("target-*_*")):
    case_output = output / case.name
    case_output.mkdir()
    sources = {
        "rsmith-rmc": case / "rsmith-joint/run.toml",
        "rsmith-pedone-w3": production_root / case.name / "pedone-w3/run.toml",
        "rsmith-pace-w3": production_root / case.name / "pace-w3/run.toml",
        "rsmith-epsr": case / "rsmith-epsr-joint/run.toml",
    }
    manifest["cases"][case.name] = {}
    for method, source in sources.items():
        if not source.is_file():
            raise SystemExit(f"missing prerequisite {source}")
        prepared = []
        for seed in protocol["budget"]["seeds"]:
            run = case_output / method / f"seed-{seed}"
            run.mkdir(parents=True)
            (run / "run.toml").write_text(
                rewrite(
                    source.read_text(),
                    int(seed),
                    int(protocol["budget"]["moves"]),
                )
            )
            prepared.append(run.name)
        manifest["cases"][case.name][method] = prepared

(output / "manifest.json").write_text(
    json.dumps(manifest, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(manifest, indent=2, sort_keys=True))
