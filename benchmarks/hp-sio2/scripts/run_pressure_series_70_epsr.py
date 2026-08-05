#!/usr/bin/env python3
"""Run stock EPSR26 over the frozen 0-to-70 GPa pressure series."""

from __future__ import annotations

import argparse
import json
import tomllib
from pathlib import Path

from epsr26_common import (
    calculated_targets,
    prepare_run,
    provisional_targets,
    read_iq,
    rms,
    run_epsr,
    write_epsr_data,
)


parser = argparse.ArgumentParser()
parser.add_argument("--force", action="store_true")
parser.add_argument("--only-missing", action="store_true")
args = parser.parse_args()
case_root = Path(__file__).resolve().parents[1]
root = case_root / "results/pressure-series-70"
protocol = tomllib.loads((case_root / "expected/pressure-series-70.toml").read_text())
source = case_root.parent / "epsr-silica/reference/local/upstream"
record_path = case_root.parent / "epsr-silica/reference/local/IMPORT.txt"
if not source.is_dir() or not record_path.is_file():
    raise SystemExit("import the accepted EPSR26 tutorial first")
record = dict(
    line.split("=", 1) for line in record_path.read_text().splitlines() if "=" in line
)
binary = Path(record["epsr_root"]) / "bin/epsr"
summary_path = root / "native-epsr-run-summary.json"
summary = (
    json.loads(summary_path.read_text())
    if summary_path.is_file()
    else {"status": "hp_sio2_pressure_series_70_native_epsr_runs", "steps": {}}
)
feedback = float(protocol["native_epsr26"]["feedback"])
seed = int(protocol["design"]["seed"])
for step in protocol["steps"]:
    name, step_root = step["name"], root / step["name"]
    target_root = step_root / "native-targets"
    target_root.mkdir(exist_ok=True)
    target_files = (
        target_root / "target-neutron.dat",
        target_root / "target-xray.dat",
    )
    forward_summary_path = step_root / "native-epsr-forward-summary.json"
    if args.force or not all(path.is_file() for path in target_files):
        provisional = provisional_targets(step_root / "inputs", source)
        first = step_root / "native-epsr-forward/provisional"
        prepare_run(
            source,
            first,
            step_root / "inputs/hidden-target.data",
            0,
            1,
            seed,
            provisional,
            feedback,
        )
        first_wall = run_epsr(binary, first)
        native_targets = calculated_targets(first / "DTBsilica.EPSR.v01", provisional)
        write_epsr_data(target_files[0], "EPSR-native hidden neutron i(Q)", native_targets[0])
        write_epsr_data(target_files[1], "EPSR-native hidden X-ray i(Q)", native_targets[1])
        replay = step_root / "native-epsr-forward/replay"
        prepare_run(
            source,
            replay,
            step_root / "inputs/hidden-target.data",
            0,
            1,
            seed,
            native_targets,
            feedback,
        )
        replay_wall = run_epsr(binary, replay)
        checked = calculated_targets(replay / "DTBsilica.EPSR.v01", native_targets)
        forward = {
            "first_wall_seconds": first_wall,
            "hidden_replay_neutron_rms": rms(checked[0], native_targets[0]),
            "hidden_replay_wall_seconds": replay_wall,
            "hidden_replay_xray_rms": rms(checked[1], native_targets[1]),
        }
        forward_summary_path.write_text(json.dumps(forward, indent=2, sort_keys=True) + "\n")
    else:
        native_targets = tuple(
            [row for row in read_iq(path) if 0.5 <= row[0] < 25.0]
            for path in target_files
        )
    run = step_root / "runs/native-epsr26"
    if args.only_missing and (run / "DTBsilica.EPSR.v01").is_file():
        continue
    if run.exists() and not args.force:
        raise SystemExit(f"output exists: {run} (pass --force or --only-missing)")
    prepare_run(
        source,
        run,
        step_root / "inputs/cross-start.data",
        int(protocol["native_epsr26"]["moves_per_refinement"]),
        int(protocol["native_epsr26"]["refinements"]),
        seed,
        native_targets,
        feedback,
    )
    wall = run_epsr(binary, run)
    summary["steps"][name] = {"status": "completed", "wall_seconds": wall}
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({name: summary["steps"][name]}, indent=2), flush=True)
