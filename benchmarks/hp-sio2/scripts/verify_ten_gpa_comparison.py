#!/usr/bin/env python3
"""Verify frozen 10 GPa comparison observations against local raw scores."""

from __future__ import annotations

import json
import math
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
observed = tomllib.loads(
    (case_root / "expected/ten-gpa-comparison-observed.toml").read_text()
)
score_path = case_root / "results/ten-gpa-pilot/score-summary.json"
if not score_path.is_file():
    raise SystemExit("missing local score summary; run score_ten_gpa_pilot.py")
comparison = json.loads(score_path.read_text())["comparison"]

assert comparison["decision"]["paired_hidden_rdf_wins"] == observed["decision"][
    "paired_hidden_rdf_wins"
]
assert (
    comparison["decision"]["rsmith_pace_superior_to_epsr_for_frozen_task"]
    == observed["decision"]["rsmith_pace_superior_to_epsr_for_frozen_task"]
)
for method in ("rsmith-pace-w30", "native-epsr26", "rsmith-rmc"):
    expected = observed["aggregate"][method.replace("-", "_")]
    actual = comparison["aggregate"][method]
    for observed_name, actual_name in (
        ("joint_iq_rms", "joint_iq_rms"),
        ("mean_partial_rdf_rms", "mean_partial_rdf_rms"),
        ("si_coordination_tv", "si_coordination_total_variation"),
        ("wall_s", "wall_seconds"),
    ):
        for statistic in ("median", "q1", "q3"):
            assert math.isclose(
                actual[actual_name][statistic],
                expected[observed_name][statistic],
                rel_tol=0.0,
                abs_tol=5.0e-10,
            )
    assert actual["safety_passes"] == expected["safety_passes"]
for kind in ("neutron", "xray"):
    value = comparison["native_forward"][f"hidden_replay_{kind}_rms"]
    assert math.isclose(
        value,
        observed["native_forward"][f"hidden_replay_{kind}_rms"],
        rel_tol=0.0,
        abs_tol=1.0e-15,
    )
    assert value <= observed["native_forward"]["gate_max"]
print("10 GPa SiO2 comparison verification passed")
