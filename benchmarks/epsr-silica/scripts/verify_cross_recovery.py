#!/usr/bin/env python3
"""Verify cross-recovery adapter smoke outputs against committed guards."""

from __future__ import annotations

import json
import math
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
expected = tomllib.loads((case_root / "expected/cross-recovery-smoke.toml").read_text())
results = case_root / "results/cross-recovery"
score = json.loads((results / "score-summary.json").read_text())
epsr = json.loads((results / "native-epsr-run-summary.json").read_text())
rmcprofile = json.loads((results / "rmcprofile-run-summary.json").read_text())
pdfgui = json.loads((results / "pdfgui-forward-summary.json").read_text())
guards = expected["guards"]

failures = []
for case_name, case in score["cases"].items():
    common = case["common_reciprocal_distance_from_hidden_target"]
    self_check = common["hidden_target_self_check"]
    if (
        max(
            self_check["neutron_iq_rms"],
            self_check["xray_iq_rms"],
            self_check["mean_partial_sq_rms"],
        )
        > guards["hidden_self_common_rms_max"]
    ):
        failures.append(f"{case_name}: independent hidden-target self-check failed")
    missing = set(guards["required_joint_coordinate_methods"]) - set(common)
    if missing:
        failures.append(f"{case_name}: missing coordinate outputs {sorted(missing)}")
    start_mean = 0.5 * (
        common["cross_start"]["neutron_iq_rms"] + common["cross_start"]["xray_iq_rms"]
    )
    for method in guards["required_joint_coordinate_methods"]:
        if method not in common:
            continue
        value = common[method]
        method_mean = 0.5 * (value["neutron_iq_rms"] + value["xray_iq_rms"])
        if (
            guards["require_each_joint_method_to_improve_mean_common_total_rms"]
            and not method_mean < start_mean
        ):
            failures.append(
                f"{case_name}: {method} did not improve the common total-scattering RMS"
            )
    for method_name, fit in case["fits"].items():
        wall = fit.get("wall_seconds")
        if wall is not None and (not math.isfinite(wall) or wall <= 0.0):
            failures.append(f"{case_name}: invalid wall time for {method_name}")

for case_name, values in epsr["cases"].items():
    if (
        max(values["hidden_zero_move_neutron_rms"], values["hidden_zero_move_xray_rms"])
        > guards["native_epsr_hidden_replay_rms_max"]
    ):
        failures.append(f"{case_name}: native EPSR hidden replay failed")
for case_name, values in rmcprofile["cases"].items():
    if (
        max(
            values["hidden_zero_move_neutron_rms"],
            values["hidden_zero_move_xray_fq_rms"],
        )
        > guards["rmcprofile_hidden_replay_rms_max"]
    ):
        failures.append(f"{case_name}: RMCProfile hidden replay failed")
for case_name, values in pdfgui["cases"].items():
    for radiation in ("N", "X"):
        if (
            values[radiation]["rms_over_target_range"]
            < guards["pdfgui_cross_start_rms_over_target_range_min"]
        ):
            failures.append(
                f"{case_name}: PDFgui {radiation} failed to distinguish target and cross start"
            )

if failures:
    raise SystemExit(
        "cross-recovery smoke verification failed:\n- " + "\n- ".join(failures)
    )
print("passed symmetric GAP/Pedone cross-program adapter smoke checks")
