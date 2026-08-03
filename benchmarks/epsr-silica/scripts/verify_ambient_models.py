#!/usr/bin/env python3
"""Verify the locally regenerated GAP/Pedone endpoint observations."""

from __future__ import annotations

import json
import math
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
with (case_root / "expected/ambient-model-endpoints.toml").open("rb") as stream:
    expected = tomllib.load(stream)
summary_path = case_root / "results/ambient-models/summary.json"
if not summary_path.is_file():
    raise SystemExit("run analyze_ambient_models.py first")
observed = json.loads(summary_path.read_text())

if observed["source_commit"] != expected["source_commit"]:
    raise SystemExit("ambient-model source commit does not match the frozen record")
if observed["analysis"]["upstream_rdf_sha256"] != expected["analysis"]["upstream_rdf_sha256"]:
    raise SystemExit("ambient-model upstream RDF hash does not match the frozen record")


def check_close(label, actual, wanted, tolerance=1.0e-10):
    if not math.isclose(actual, wanted, rel_tol=tolerance, abs_tol=tolerance):
        raise SystemExit(f"{label}: observed {actual}, expected {wanted}")


if (
    observed["analysis"]["rdf_max_absolute_difference_vs_upstream"]
    > expected["analysis"]["rdf_max_absolute_difference_vs_upstream_max"]
):
    raise SystemExit("independent RDF analysis no longer matches the pinned upstream curves")


for model in ("gap", "pedone"):
    actual = observed["models"][model]
    wanted = expected[model]
    if actual["sha256"] != wanted["sha256"]:
        raise SystemExit(f"{model}: structure hash does not match the frozen record")
    if actual["atom_count"] != wanted["atom_count"]:
        raise SystemExit(f"{model}: atom count does not match the frozen record")
    for key in (
        "number_density_atoms_a3",
        "mass_density_g_cm3",
    ):
        check_close(f"{model}.{key}", actual[key], wanted[key])
    for actual_key, expected_key in (
        ("Si_defect_fraction", "si_coordination_defect_fraction"),
        ("O_defect_fraction", "o_coordination_defect_fraction"),
    ):
        check_close(
            f"{model}.coordination.{actual_key}",
            actual["coordination"][actual_key],
            wanted[expected_key],
        )
    for angle, expected_prefix in (("O-Si-O", "o_si_o"), ("Si-O-Si", "si_o_si")):
        check_close(
            f"{model}.{angle}.mean",
            actual["angles_degrees"][angle]["mean"],
            wanted[f"{expected_prefix}_mean_degrees"],
        )
        check_close(
            f"{model}.{angle}.standard_deviation",
            actual["angles_degrees"][angle]["standard_deviation"],
            wanted[f"{expected_prefix}_standard_deviation_degrees"],
        )
    for pair, expected_prefix in (("Si-Si", "si_si"), ("Si-O", "si_o"), ("O-O", "o_o")):
        check_close(
            f"{model}.{pair}.peak",
            actual["rdf_first_peak"][pair]["r_a"],
            wanted[f"{expected_prefix}_peak_a"],
        )
    ring = actual["shortest_ring_diagnostic"]
    if (
        ring["network_edges"] != wanted["si_network_edges"]
        or ring["acyclic_edges"] != wanted["si_network_acyclic_edges"]
    ):
        raise SystemExit(f"{model}: shortest-ring graph size changed")
    expected_cycles = {
        key.removeprefix("edge_cycle_size_"): value
        for key, value in wanted.items()
        if key.startswith("edge_cycle_size_")
    }
    if ring["cycle_size_by_edge"] != expected_cycles:
        raise SystemExit(f"{model}: shortest-ring distribution changed")

print("Ambient GAP/Pedone endpoint observations match the frozen record.")
