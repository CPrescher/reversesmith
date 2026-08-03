#!/usr/bin/env python3
"""Prepare the frozen Erhard-2024 PACE weight pilot."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import tomllib
from pathlib import Path


def weight_label(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def rewrite_config(text: str, weight: float, moves: int, model: Path) -> str:
    for pattern, replacement in (
        (r"(?m)^max_moves = \d+$", f"max_moves = {moves}"),
        (r"(?m)^print_every = \d+$", f"print_every = {moves}"),
    ):
        text, count = re.subn(pattern, replacement, text)
        if count != 1:
            raise ValueError(f"expected one match for {pattern}")
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
    (case_root / "expected/hrmc-pace2024-weight-pilot.toml").read_text()
)
model = case_root / "reference/local/public-ace2024/SiOx_potential.yace"
if not model.is_file():
    raise SystemExit(f"missing verified Erhard-2024 model {model}")
cross_root = case_root / "results/cross-recovery"
output = case_root / "results/hrmc-pace2024-weight-pilot"
if output.exists():
    if not args.force:
        raise SystemExit(f"output exists: {output} (pass --force to replace it)")
    shutil.rmtree(output)
output.mkdir(parents=True)

manifest = {"status": "prepared", "cases": {}}
for case in sorted(cross_root.glob("target-*_*")):
    source = case / "rsmith-gap-joint/run.toml"
    if not source.is_file():
        raise SystemExit(f"missing prerequisite {source}")
    case_output = output / case.name
    case_output.mkdir()
    prepared = []
    for weight in protocol["weights"]["pace"]:
        run = case_output / f"pace-w{weight_label(weight)}"
        run.mkdir()
        (run / "run.toml").write_text(
            rewrite_config(
                source.read_text(),
                float(weight),
                int(protocol["budget"]["moves"]),
                model,
            )
        )
        prepared.append(run.name)
    manifest["cases"][case.name] = prepared

(output / "manifest.json").write_text(
    json.dumps(manifest, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(manifest, indent=2, sort_keys=True))
