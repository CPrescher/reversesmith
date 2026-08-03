#!/usr/bin/env python3
"""Prepare the frozen Pedone/GAP HRMC weight-pilot configurations."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import tomllib
from pathlib import Path


def weight_label(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def rewrite_budget(text: str, moves: int) -> str:
    text, count = re.subn(r"(?m)^max_moves = \d+$", f"max_moves = {moves}", text)
    if count != 1:
        raise ValueError("expected one max_moves setting")
    text, count = re.subn(r"(?m)^print_every = \d+$", f"print_every = {moves}", text)
    if count != 1:
        raise ValueError("expected one print_every setting")
    return text


def rewrite_weight(text: str, value: float) -> str:
    text, count = re.subn(r"(?m)^weight = [0-9.eE+-]+$", f"weight = {value:g}", text)
    if count != 1:
        raise ValueError("expected one energy weight setting")
    return text


parser = argparse.ArgumentParser()
parser.add_argument("--force", action="store_true")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
protocol = tomllib.loads((case_root / "expected/hrmc-weight-sweep.toml").read_text())
cross_root = case_root / "results/cross-recovery"
output = case_root / "results/hrmc-weight-sweep"
if output.exists():
    if not args.force:
        raise SystemExit(f"output exists: {output} (pass --force to replace it)")
    shutil.rmtree(output)
output.mkdir(parents=True)

moves = int(protocol["budget"]["moves"])
manifest = {"status": "prepared", "moves": moves, "cases": {}}
for case in sorted(cross_root.glob("target-*_*")):
    case_output = output / case.name
    case_output.mkdir()
    prepared = []
    sources = {
        "pure": case / "rsmith-joint/run.toml",
        "pedone": case / "rsmith-pedone-joint/run.toml",
        "gap": case / "rsmith-gap-joint/run.toml",
    }
    for source in sources.values():
        if not source.is_file():
            raise SystemExit(f"missing prerequisite {source}")

    pure = case_output / "pure-rmc"
    pure.mkdir()
    (pure / "run.toml").write_text(rewrite_budget(sources["pure"].read_text(), moves))
    prepared.append(pure.name)

    for model in ("pedone", "gap"):
        for weight in protocol["weights"][model]:
            run = case_output / f"{model}-w{weight_label(weight)}"
            run.mkdir()
            text = rewrite_budget(sources[model].read_text(), moves)
            (run / "run.toml").write_text(rewrite_weight(text, float(weight)))
            prepared.append(run.name)
    manifest["cases"][case.name] = prepared

(output / "manifest.json").write_text(
    json.dumps(manifest, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(manifest, indent=2, sort_keys=True))
print(f"HRMC weight-pilot configurations: {output}")
