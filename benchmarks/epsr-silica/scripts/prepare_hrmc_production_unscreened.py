#!/usr/bin/env python3
"""Prepare the frozen joint-acceptance Pedone/GAP/PACE production bracket."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import tomllib
from pathlib import Path


def weight_label(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def rewrite_common(text: str, weight: float, moves: int) -> str:
    for pattern, replacement in (
        (r"(?m)^max_moves = \d+$", f"max_moves = {moves}"),
        (r"(?m)^print_every = \d+$", f"print_every = {moves}"),
        (r"(?m)^weight = [0-9.eE+-]+$", f"weight = {weight:g}"),
    ):
        text, count = re.subn(pattern, replacement, text)
        if count != 1:
            raise ValueError(f"expected one match for {pattern}")
    return text


def rewrite_pace(text: str, weight: float, moves: int, model: Path) -> str:
    text = rewrite_common(text, weight, moves)
    text, count = re.subn(
        r"(?ms)^\[ml_potential\]\n.*\Z",
        "[ml_potential]\n"
        'backend = "pace_native"\n'
        f'model = "{model}"\n'
        f"weight = {weight:g}\n",
        text,
    )
    if count != 1:
        raise ValueError("expected one terminal ml_potential table")
    return text


parser = argparse.ArgumentParser()
parser.add_argument("--force", action="store_true")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
protocol = tomllib.loads(
    (case_root / "expected/hrmc-production-unscreened.toml").read_text()
)
pace_model = case_root / "reference/local/public-ace2024/SiOx_potential.yace"
if not pace_model.is_file():
    raise SystemExit(f"missing verified Erhard-2024 model {pace_model}")
cross_root = case_root / "results/cross-recovery"
output = case_root / "results/hrmc-production-unscreened"
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
    for regularizer, source_name in (
        ("pedone", "rsmith-pedone-joint"),
        ("gap", "rsmith-gap-joint"),
        ("pace", "rsmith-gap-joint"),
    ):
        source = case / source_name / "run.toml"
        if not source.is_file():
            raise SystemExit(f"missing prerequisite {source}")
        for weight in protocol["weights"][regularizer]:
            run = case_output / f"{regularizer}-w{weight_label(weight)}"
            run.mkdir()
            if regularizer == "pace":
                text = rewrite_pace(
                    source.read_text(),
                    float(weight),
                    int(protocol["budget"]["moves"]),
                    pace_model,
                )
            else:
                text = rewrite_common(
                    source.read_text(),
                    float(weight),
                    int(protocol["budget"]["moves"]),
                )
            (run / "run.toml").write_text(text)
            prepared.append(run.name)
    manifest["cases"][case.name] = prepared

(output / "manifest.json").write_text(
    json.dumps(manifest, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(manifest, indent=2, sort_keys=True))
