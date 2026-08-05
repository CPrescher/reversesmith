#!/usr/bin/env python3
"""Verify the frozen 0-to-70 GPa observations against local raw outputs."""

from __future__ import annotations

import json
import math
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
root = case_root / "results/pressure-series-70"
observed = tomllib.loads(
    (case_root / "expected/pressure-series-70-observed.toml").read_text()
)
scores = json.loads((root / "score-summary.json").read_text())
step_names = list(scores["steps"])
method_names = {
    "rsmith_rmc": "rsmith-rmc",
    "rsmith_pace_w30": "rsmith-pace-w30",
    "native_epsr26": "native-epsr26",
    "native_rmcprofile": "native-rmcprofile",
}
score_names = {
    "mean_partial_rdf_rms": "mean_partial_rdf_rms",
    "joint_iq_rms": "joint_iq_rms",
    "si_coordination_tv": "si_coordination_total_variation",
    "mean_lower_tail_error_a": "mean_lower_tail_quantile_error_a",
    "wall_s": "wall_seconds",
}


def close(first, second):
    if not math.isclose(float(first), float(second), rel_tol=1e-11, abs_tol=1e-12):
        raise AssertionError(f"{first} != {second}")


for recorded_name, score_name in score_names.items():
    if recorded_name == "wall_s":
        continue
    for expected, step_name in zip(observed["start"][recorded_name], step_names):
        close(expected, scores["steps"][step_name]["start"][score_name])

for recorded_method, score_method in method_names.items():
    record = observed["methods"][recorded_method]
    for recorded_name, score_name in score_names.items():
        for expected, step_name in zip(record[recorded_name], step_names):
            close(expected, scores["steps"][step_name]["methods"][score_method][score_name])
    for pair, recorded_name in (
        ("Si-Si", "minimum_si_si_a"),
        ("Si-O", "minimum_si_o_a"),
        ("O-O", "minimum_o_o_a"),
    ):
        for expected, step_name in zip(record[recorded_name], step_names):
            close(
                expected,
                scores["steps"][step_name]["methods"][score_method][
                    "minimum_distance_a"
                ][pair],
            )
    for expected, step_name in zip(
        record["hidden_rdf_gap_recovered_fraction"], step_names
    ):
        step = scores["steps"][step_name]
        start = step["start"]["mean_partial_rdf_rms"]
        final = step["methods"][score_method]["mean_partial_rdf_rms"]
        close(expected, (start - final) / start)

for step_name in step_names:
    audit = json.loads(
        (root / step_name / "runs/native-rmcprofile/move-audit.json").read_text()
    )
    if audit["total_moves"] != observed["moves"]:
        raise AssertionError(f"RMCProfile move audit failed for {step_name}")
    epsr_forward = json.loads(
        (root / step_name / "native-epsr-forward-summary.json").read_text()
    )
    rmcprofile_forward = json.loads(
        (root / step_name / "native-rmcprofile-forward-summary.json").read_text()
    )
    if epsr_forward["hidden_replay_neutron_rms"] > 2e-8:
        raise AssertionError(f"EPSR neutron replay failed for {step_name}")
    if epsr_forward["hidden_replay_xray_rms"] > 5e-8:
        raise AssertionError(f"EPSR X-ray replay failed for {step_name}")
    if rmcprofile_forward["hidden_replay_neutron_rms"] != 0.0:
        raise AssertionError(f"RMCProfile neutron replay failed for {step_name}")
    if rmcprofile_forward["hidden_replay_xray_fq_rms"] != 0.0:
        raise AssertionError(f"RMCProfile X-ray replay failed for {step_name}")

decision = observed["decision"]
if decision["pace_beats_native_epsr26_hidden_rdf"] != 7:
    raise AssertionError("frozen EPSR comparison count changed")
if decision["pace_beats_native_rmcprofile_hidden_rdf"] != 7:
    raise AssertionError("frozen RMCProfile comparison count changed")
print("verified 0-to-70 GPa four-program pressure series")
