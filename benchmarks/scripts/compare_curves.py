#!/usr/bin/env python3
"""Compare two tabulated curves on the reference grid using linear interpolation."""

from __future__ import annotations

import argparse
import bisect
import json
import math
from pathlib import Path


def read_columns(path: Path, x_column: int, y_column: int, sigma_column: int | None):
    required = max(x_column, y_column, sigma_column or 0)
    rows = []
    for line_number, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < required:
            raise ValueError(f"{path}:{line_number}: expected at least {required} columns")
        try:
            x = float(fields[x_column - 1])
            y = float(fields[y_column - 1])
            sigma = float(fields[sigma_column - 1]) if sigma_column else None
        except ValueError as error:
            raise ValueError(f"{path}:{line_number}: non-numeric comparison row") from error
        if not math.isfinite(x) or not math.isfinite(y):
            raise ValueError(f"{path}:{line_number}: non-finite coordinate or value")
        if sigma is not None and (not math.isfinite(sigma) or sigma <= 0.0):
            raise ValueError(f"{path}:{line_number}: uncertainty must be finite and positive")
        rows.append((x, y, sigma))
    if len(rows) < 2:
        raise ValueError(f"{path}: expected at least two numeric rows")
    if any(rows[index][0] <= rows[index - 1][0] for index in range(1, len(rows))):
        raise ValueError(f"{path}: coordinates must be strictly increasing")
    return rows


def interpolate(rows, x):
    xs = [row[0] for row in rows]
    upper = bisect.bisect_left(xs, x)
    if upper < len(rows) and rows[upper][0] == x:
        return rows[upper][1]
    if upper == 0 or upper == len(rows):
        raise ValueError(f"coordinate {x} lies outside candidate range")
    x0, y0, _ = rows[upper - 1]
    x1, y1, _ = rows[upper]
    fraction = (x - x0) / (x1 - x0)
    return y0 + fraction * (y1 - y0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--reference-x-column", type=int, default=1)
    parser.add_argument("--reference-y-column", type=int, default=2)
    parser.add_argument("--reference-sigma-column", type=int)
    parser.add_argument("--candidate-x-column", type=int, default=1)
    parser.add_argument("--candidate-y-column", type=int, default=2)
    parser.add_argument("--fit-min", type=float, default=-math.inf)
    parser.add_argument("--fit-max", type=float, default=math.inf)
    parser.add_argument("--max-abs-tolerance", type=float)
    parser.add_argument("--rms-tolerance", type=float)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    reference = read_columns(
        args.reference,
        args.reference_x_column,
        args.reference_y_column,
        args.reference_sigma_column,
    )
    candidate = read_columns(args.candidate, args.candidate_x_column, args.candidate_y_column, None)
    overlap_min = max(args.fit_min, candidate[0][0])
    overlap_max = min(args.fit_max, candidate[-1][0])
    selected = [row for row in reference if overlap_min <= row[0] <= overlap_max]
    if len(selected) < 2:
        raise ValueError("fewer than two reference points lie in the comparison range")

    residuals = [interpolate(candidate, x) - y for x, y, _ in selected]
    reference_values = [row[1] for row in selected]
    dynamic_range = max(reference_values) - min(reference_values)
    result = {
        "reference": str(args.reference),
        "candidate": str(args.candidate),
        "n_points": len(residuals),
        "fit_min": selected[0][0],
        "fit_max": selected[-1][0],
        "mean_error": sum(residuals) / len(residuals),
        "mean_absolute_error": sum(abs(value) for value in residuals) / len(residuals),
        "rms_error": math.sqrt(sum(value * value for value in residuals) / len(residuals)),
        "max_absolute_error": max(abs(value) for value in residuals),
        "reference_dynamic_range": dynamic_range,
    }
    if dynamic_range > 0.0:
        result["rms_error_over_dynamic_range"] = result["rms_error"] / dynamic_range
        result["max_absolute_error_over_dynamic_range"] = (
            result["max_absolute_error"] / dynamic_range
        )
    if args.reference_sigma_column:
        normalized = [residual / row[2] for residual, row in zip(residuals, selected)]
        result["reduced_chi_squared"] = sum(value * value for value in normalized) / len(normalized)

    failures = []
    if args.max_abs_tolerance is not None and result["max_absolute_error"] > args.max_abs_tolerance:
        failures.append("max_absolute_error")
    if args.rms_tolerance is not None and result["rms_error"] > args.rms_tolerance:
        failures.append("rms_error")
    result["passed"] = not failures
    result["failed_tolerances"] = failures
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
