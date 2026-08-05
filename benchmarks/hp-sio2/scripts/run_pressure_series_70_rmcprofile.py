#!/usr/bin/env python3
"""Run native RMCProfile over the frozen 0-to-70 GPa pressure series."""

from __future__ import annotations

import argparse
import hashlib
import json
import tomllib
from pathlib import Path

from rmcprofile_common import (
    curve_rms,
    numeric_curve,
    prepare_run,
    read_calculated,
    run,
    run_exact,
    write_curve,
)


parser = argparse.ArgumentParser()
parser.add_argument(
    "--executable",
    type=Path,
    default=Path("/Applications/RMCProfile.app/Contents/MacOS/exe/rmcprofile.x"),
)
parser.add_argument("--force", action="store_true")
parser.add_argument("--only-missing", action="store_true")
parser.add_argument(
    "--forward-only",
    action="store_true",
    help="validate native hidden-coordinate forward/replay without refining",
)
args = parser.parse_args()
case_root = Path(__file__).resolve().parents[1]
root = case_root / "results/pressure-series-70"
protocol = tomllib.loads((case_root / "expected/pressure-series-70.toml").read_text())
if not args.executable.is_file():
    raise SystemExit(f"RMCProfile executable not found: {args.executable}")
digest = hashlib.sha256(args.executable.read_bytes()).hexdigest()
if digest != protocol["native_rmcprofile"]["executable_sha256"]:
    raise SystemExit("RMCProfile executable differs from frozen hash")
seed, moves = int(protocol["design"]["seed"]), int(protocol["design"]["moves"])
minima = protocol["native_rmcprofile"]["minimum_distances_a"]
summary_path = root / "native-rmcprofile-run-summary.json"
summary = (
    json.loads(summary_path.read_text())
    if summary_path.is_file()
    else {"status": "hp_sio2_pressure_series_70_native_rmcprofile_runs", "steps": {}}
)
for step in protocol["steps"]:
    name, step_root = step["name"], root / step["name"]
    input_root = step_root / "inputs"
    target_root = step_root / "rmcprofile-native-targets"
    target_root.mkdir(exist_ok=True)
    target_files = target_root / "target-neutron.sq", target_root / "target-xray.fq"
    forward_path = step_root / "native-rmcprofile-forward-summary.json"
    if args.force or not all(path.is_file() for path in target_files):
        neutron_iq = numeric_curve(input_root / "target-neutron-iq.dat")
        xray_iq = numeric_curve(input_root / "target-xray-iq.dat")
        provisional = (
            [(q, value + 1.0) for q, value in neutron_iq],
            [(q, q * value) for q, value in xray_iq],
        )
        first = step_root / "native-rmcprofile-forward/provisional"
        prepare_run(first, input_root / "hidden-target.rmc6f", provisional, 0, 0.0, minima)
        first_wall, first_log = run(args.executable, first, seed)
        (first / "driver.log").write_text(first_log)
        native_targets = (
            read_calculated(first / "cross_SQ1.csv"),
            read_calculated(first / "cross_FQ1.csv"),
        )
        write_curve(target_files[0], "RMCProfile-native hidden neutron S(Q)", native_targets[0])
        write_curve(target_files[1], "RMCProfile-native hidden X-ray F(Q)", native_targets[1])
        replay = step_root / "native-rmcprofile-forward/replay"
        prepare_run(replay, input_root / "hidden-target.rmc6f", native_targets, 0, 0.0, minima)
        replay_wall, replay_log = run(args.executable, replay, seed)
        (replay / "driver.log").write_text(replay_log)
        checked = (
            read_calculated(replay / "cross_SQ1.csv"),
            read_calculated(replay / "cross_FQ1.csv"),
        )
        forward = {
            "first_wall_seconds": first_wall,
            "hidden_replay_neutron_rms": curve_rms(checked[0], native_targets[0]),
            "hidden_replay_wall_seconds": replay_wall,
            "hidden_replay_xray_fq_rms": curve_rms(checked[1], native_targets[1]),
        }
        forward_path.write_text(json.dumps(forward, indent=2, sort_keys=True) + "\n")
    else:
        native_targets = numeric_curve(target_files[0]), numeric_curve(target_files[1])
    if args.forward_only:
        summary["steps"][name] = {
            "status": "forward_replay_validated",
            **json.loads(forward_path.read_text()),
        }
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        print(json.dumps({name: summary["steps"][name]}, indent=2), flush=True)
        continue
    run_dir = step_root / "runs/native-rmcprofile"
    if args.only_missing and (run_dir / "move-audit.json").is_file():
        continue
    if run_dir.exists() and not args.force:
        raise SystemExit(f"output exists: {run_dir} (pass --force or --only-missing)")
    prepare_run(
        run_dir,
        input_root / "cross-start.rmc6f",
        native_targets,
        moves,
        0.05,
        minima,
    )
    audit = run_exact(args.executable, run_dir, seed, moves, native_targets, minima)
    summary["steps"][name] = {"status": "completed", **audit}
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({name: summary["steps"][name]}, indent=2), flush=True)
