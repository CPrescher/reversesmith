#!/usr/bin/env python3
"""Verify a local independent-start summary against the frozen observations."""

from __future__ import annotations

import json
import math
import tomllib
from pathlib import Path


ABS_TOLERANCE = 5.0e-12


def check_close(label: str, observed: float, expected: float, tolerance=ABS_TOLERANCE):
    if not math.isclose(observed, expected, rel_tol=0.0, abs_tol=tolerance):
        raise AssertionError(
            f"{label}: observed {observed:.15g}, expected {expected:.15g}"
        )


def main():
    case_root = Path(__file__).resolve().parents[1]
    summary_path = case_root / "results/independent-starts/summary.json"
    expected_path = case_root / "expected/independent-starts.toml"
    if not summary_path.is_file():
        raise SystemExit(
            "run run_independent_starts.py before verifying the local result"
        )

    summary = json.loads(summary_path.read_text())
    with expected_path.open("rb") as stream:
        expected = tomllib.load(stream)

    protocol = expected["protocol"]
    if summary["protocol"]["seeds"] != protocol["seeds"]:
        raise AssertionError("seed list differs from frozen protocol")
    if summary["protocol"]["checkpoints"] != protocol["checkpoints"]:
        raise AssertionError("checkpoint list differs from frozen protocol")
    if summary["protocol"]["replicates"] != protocol["replicates"]:
        raise AssertionError("replicate count differs from frozen protocol")

    replicate_expected = expected["observed"]["replicates"]
    if summary["start_sha256"] != {
        str(seed): digest
        for seed, digest in zip(
            replicate_expected["seeds"], replicate_expected["start_sha256"], strict=True
        )
    }:
        raise AssertionError("generated-start hashes differ from the frozen record")

    for index, trace in enumerate(summary["traces"]):
        if trace["seed"] != replicate_expected["seeds"][index]:
            raise AssertionError(f"replicate {index} seed differs")
        for method, field in (
            ("native", "native_total_fit_by_checkpoint"),
            ("rsmith", "rsmith_total_fit_by_checkpoint"),
        ):
            observed_trace = [
                row[method]["total_fit"]["rms_over_dynamic_range"]
                for row in trace["checkpoints"]
            ]
            for epoch, observed_value, expected_value in zip(
                protocol["checkpoints"],
                observed_trace,
                replicate_expected[field][index],
                strict=True,
            ):
                check_close(
                    f"seed {trace['seed']} {method} epoch {epoch} residual",
                    observed_value,
                    expected_value,
                )

    observed_record = expected["observed"]
    for method in ("native", "rsmith"):
        frozen_final = observed_record["final"][method]
        local_final = summary["checkpoint_distributions"]["40"][method]
        for metric in ("total_fit_rms_over_dynamic_range", "wall_seconds"):
            for statistic in ("min", "median", "max", "mean", "sample_stddev"):
                check_close(
                    f"{method} final {metric} {statistic}",
                    local_final[metric][statistic],
                    frozen_final[metric][statistic],
                    tolerance=5.0e-9 if metric == "wall_seconds" else ABS_TOLERANCE,
                )

    for suffix in ("0_040", "0_030", "0_025", "0_020", "0_015"):
        frozen_target = observed_record["time_to_target"][f"target_{suffix}"]
        target_key = f"{float(suffix.replace('_', '.')):.6f}"
        local_target = summary["time_to_target"][target_key]
        for method in ("native", "rsmith"):
            if (
                local_target[method]["successes"]
                != frozen_target[f"{method}_successes"]
            ):
                raise AssertionError(f"target {suffix} {method} success count differs")
            wall_key = f"{method}_median_wall_seconds"
            if wall_key in frozen_target:
                check_close(
                    f"target {suffix} {method} median wall time",
                    local_target[method]["wall_seconds"]["median"],
                    frozen_target[wall_key],
                    tolerance=5.0e-9,
                )
        if "rsmith_speedup_over_native" in frozen_target:
            check_close(
                f"target {suffix} speedup",
                local_target["rsmith_speedup_over_native_median"],
                frozen_target["rsmith_speedup_over_native"],
            )

    instability = observed_record["native_instability"]
    failed_trace = next(
        trace for trace in summary["traces"] if trace["seed"] == instability["seed"]
    )["checkpoints"][-1]["native"]
    check_close(
        "native instability minimum distance",
        failed_trace["configuration"]["nearest_neighbor_angstrom"]["min"],
        instability["minimum_distance_angstrom"],
    )
    check_close(
        "native instability energy",
        failed_trace["combined_energy_kj_mol_per_atom"],
        instability["combined_energy_kj_mol_per_atom"],
    )

    print(
        "Independent-start benchmark verified: 10 paired traces, five targets, "
        "final distributions, start hashes, and native-instability diagnostics match."
    )


if __name__ == "__main__":
    main()
