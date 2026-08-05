#!/usr/bin/env python3
"""Verify frozen incremental pressure-step observations."""

from __future__ import annotations

import json
import math
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
observed = tomllib.loads((case_root / "expected/pressure-steps-observed.toml").read_text())
score_path = case_root / "results/pressure-steps/score-summary.json"
if not score_path.is_file():
    raise SystemExit("missing local pressure-step score summary")
scores = json.loads(score_path.read_text())
for key in ("all_steps_informative", "pace_recovery_advantage", "lower_tail_guard"):
    assert scores["decision"][key] == observed["decision"][key]
method_names = {
    "rsmith_pace_w30": "rsmith-pace-w30",
    "native_epsr26": "native-epsr26",
    "rsmith_rmc": "rsmith-rmc",
}
for step_name, expected_step in observed["steps"].items():
    actual_step = scores["steps"][step_name]
    for expected_name, actual_name in method_names.items():
        expected = expected_step[expected_name]
        actual = actual_step["methods"][actual_name]
        for metric in (
            "joint_iq_rms",
            "mean_partial_rdf_rms",
            "mean_lower_tail_quantile_error_a",
        ):
            assert math.isclose(
                actual[metric], expected[metric], rel_tol=0.0, abs_tol=5.0e-10
            )
        for pair, value in expected["minimum_distance_a"].items():
            assert math.isclose(
                actual["minimum_distance_a"][pair],
                value,
                rel_tol=0.0,
                abs_tol=5.0e-10,
            )
    for kind in ("neutron", "xray"):
        actual = actual_step["native_forward"][f"hidden_replay_{kind}_rms"]
        expected = expected_step["native_forward"][f"{kind}_rms"]
        assert math.isclose(actual, expected, rel_tol=0.0, abs_tol=1.0e-15)
        assert actual < 1.0e-6
print("incremental pressure-step verification passed")
