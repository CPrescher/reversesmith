#!/usr/bin/env python3
"""Add the frozen Erhard-2024 PACE weight follow-up to the pilot directory."""

from __future__ import annotations

import argparse
import re
import tomllib
from pathlib import Path


def weight_label(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def rewrite_weight(text: str, value: float) -> str:
    text, count = re.subn(r"(?m)^weight = [0-9.eE+-]+$", f"weight = {value:g}", text)
    if count != 1:
        raise ValueError("expected one energy weight setting")
    return text


parser = argparse.ArgumentParser()
parser.add_argument("--force", action="store_true")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
protocol = tomllib.loads(
    (case_root / "expected/hrmc-pace2024-weight-followup.toml").read_text()
)
root = case_root / "results/hrmc-pace2024-weight-pilot"
if not root.is_dir():
    raise SystemExit("prepare and run the parent PACE weight pilot first")

for case in sorted(root.glob("target-*_*")):
    source = case / "pace-w3/run.toml"
    if not source.is_file():
        raise SystemExit(f"missing parent configuration {source}")
    for weight in protocol["weights"]["pace"]:
        run = case / f"pace-w{weight_label(weight)}"
        run.mkdir(exist_ok=args.force)
        config = run / "run.toml"
        if config.exists() and not args.force:
            raise SystemExit(
                f"configuration exists: {config} (pass --force to replace it)"
            )
        config.write_text(rewrite_weight(source.read_text(), float(weight)))
        print(config)
