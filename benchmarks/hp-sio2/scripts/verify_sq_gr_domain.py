#!/usr/bin/env python3
"""Verify the frozen S(Q)/g(r) domain observations against raw outputs."""

from __future__ import annotations

import json
import math
import re
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
root = case_root / "results/sq-gr-domain"
observed = tomllib.loads(
    (case_root / "expected/sq-gr-domain-observed.toml").read_text()
)
scores = json.loads((root / "score-summary.json").read_text())
method_names = {
    "sq_only": "sq-only",
    "gr_only": "gr-only",
    "sq_plus_gr": "sq-plus-gr",
}
score_names = {
    "full_partial_rdf_rms": "joint_hidden_partial_rdf_rms",
    "local_partial_rdf_rms": "local_hidden_partial_rdf_rms",
    "xray_iq_rms": "xray_iq_rms",
    "xray_gr_local_rms": "xray_gr_local_rms",
    "si_coordination_tv": "si_coordination_total_variation",
    "mean_lower_tail_error_a": "mean_lower_tail_quantile_error_a",
    "energy_change_ev_per_atom": "energy_change_ev_per_atom",
    "accepted_moves": "accepted_moves",
    "wall_s": "wall_seconds",
}


def close(first, second):
    if not math.isclose(float(first), float(second), rel_tol=1e-11, abs_tol=1e-12):
        raise AssertionError(f"{first} != {second}")


for recorded_name, score_name in score_names.items():
    if recorded_name in {"energy_change_ev_per_atom", "accepted_moves", "wall_s"}:
        continue
    close(observed["start"][recorded_name], scores["start"][score_name])

for recorded_method, score_method in method_names.items():
    record = observed["methods"][recorded_method]
    score = scores["methods"][score_method]
    for recorded_name, score_name in score_names.items():
        close(record[recorded_name], score[score_name])
    log = (root / "runs" / score_method / "rsmith.log").read_text()
    initial = re.search(r"Initial chi2 =\s*([+\-0-9.eE]+)", log)
    final = re.search(r"Final chi2 =\s*([+\-0-9.eE]+)", log)
    if initial is None or final is None:
        raise AssertionError(f"missing chi2 audit for {score_method}")
    close(record["initial_chi2"], initial.group(1))
    close(record["final_chi2"], final.group(1))

joint_log = (root / "runs/sq-plus-gr/rsmith.log").read_text()
components = re.search(
    r"Initial chi2 =\s*[+\-0-9.eE]+ \(sq:\s*([+\-0-9.eE]+), gr:\s*([+\-0-9.eE]+)\)",
    joint_log,
)
if components is None:
    raise AssertionError("missing joint initial-cost decomposition")
close(observed["methods"]["sq_plus_gr"]["initial_sq_chi2"], components.group(1))
close(observed["methods"]["sq_plus_gr"]["initial_gr_chi2"], components.group(2))
if abs(float(components.group(1)) - float(components.group(2))) > 0.1:
    raise AssertionError("joint domain costs are not balanced")
print("verified S(Q)/g(r) domain comparison")
