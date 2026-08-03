#!/usr/bin/env python3
"""Verify the local DTBsilicaNX forward result against regression guards."""

from __future__ import annotations

import json
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
expected = tomllib.loads((case_root / "expected/native-forward.toml").read_text())
summary_path = case_root / "results/rsmith-zero-move/summary.json"
if not summary_path.is_file():
    raise SystemExit("run scripts/compare_forward.py first")
summary = json.loads(summary_path.read_text())
limits = expected["regression_limits"]

checks = []
for pair, values in summary["partial_iq"].items():
    checks.append(
        (
            f"partial_iq.{pair}",
            values["rms_over_dynamic_range"],
            limits["partial_iq_rms_over_dynamic_range_max"],
        )
    )
checks.extend(
    (
        name,
        summary["totals"][kind][field]["rms_over_dynamic_range"],
        limits["epsr_normalized_total_iq_rms_over_dynamic_range_max"],
    )
    for name, kind, field in (
        ("total.neutron", "neutron_iq", "epsr_nrtype5_from_rsmith_partials"),
        (
            "total.xray",
            "xray_iq",
            "epsr_single_atom_normalized_from_rsmith_partials",
        ),
    )
)
for pair, values in summary["partial_gr_rebinned_0.12A"].items():
    checks.append(
        (
            f"partial_gr.{pair}",
            values["rms_over_dynamic_range"],
            limits["partial_gr_rebinned_rms_over_dynamic_range_max"],
        )
    )
for pair, values in summary["rsmith_gr_vs_independent_oracle"].items():
    checks.append(
        (
            f"oracle_gr.{pair}",
            values["rms_over_dynamic_range"],
            limits["oracle_gr_rms_over_dynamic_range_max"],
        )
    )

failures = [
    f"{name}: observed {observed:.8g} > limit {limit:.8g}"
    for name, observed, limit in checks
    if observed > limit
]
for name, observed, limit in checks:
    print(f"{name:24s} {observed:.8g} <= {limit:.8g}")
if failures:
    raise SystemExit("forward regression checks failed:\n  " + "\n  ".join(failures))
print("passed DTBsilicaNX deterministic forward regression checks")
