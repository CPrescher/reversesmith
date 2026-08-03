#!/usr/bin/env python3
"""Prepare the frozen 6,000-move Pedone/GAP HRMC production bracket."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import tomllib
from pathlib import Path


def weight_label(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def rewrite_config(text: str, weight: float, moves: int) -> str:
    replacements = (
        (r"(?m)^max_moves = \d+$", f"max_moves = {moves}"),
        (r"(?m)^print_every = \d+$", f"print_every = {moves}"),
        (r"(?m)^weight = [0-9.eE+-]+$", f"weight = {weight:g}"),
    )
    for pattern, replacement in replacements:
        text, count = re.subn(pattern, replacement, text)
        if count != 1:
            raise ValueError(f"expected one match for {pattern}")
    text, count = re.subn(
        r"(?m)^adjust_step_every = \d+$",
        "adjust_step_every = 1000\ndelayed_acceptance = true\nenergy_calibration_moves = 0",
        text,
    )
    if count != 1:
        raise ValueError("expected one adjust_step_every setting")
    return text


parser = argparse.ArgumentParser()
parser.add_argument("--force", action="store_true")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
protocol = tomllib.loads(
    (case_root / "expected/hrmc-production-bracket.toml").read_text()
)
cross_root = case_root / "results/cross-recovery"
output = case_root / "results/hrmc-production-bracket"
if output.exists():
    if not args.force:
        raise SystemExit(f"output exists: {output} (pass --force to replace it)")
    shutil.rmtree(output)
output.mkdir(parents=True)

manifest = {"status": "prepared", "cases": {}}
for case in sorted(cross_root.glob("target-*_*")):
    case_output = output / case.name
    case_output.mkdir()
    prepared = []
    for model, source_name in (
        ("pedone", "rsmith-pedone-joint"),
        ("gap", "rsmith-gap-joint"),
    ):
        source = case / source_name / "run.toml"
        if not source.is_file():
            raise SystemExit(f"missing prerequisite {source}")
        for weight in protocol["weights"][model]:
            run = case_output / f"{model}-w{weight_label(weight)}"
            run.mkdir()
            (run / "run.toml").write_text(
                rewrite_config(
                    source.read_text(), float(weight), int(protocol["budget"]["moves"])
                )
            )
            prepared.append(run.name)
    manifest["cases"][case.name] = prepared

(output / "manifest.json").write_text(
    json.dumps(manifest, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(manifest, indent=2, sort_keys=True))
