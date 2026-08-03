#!/usr/bin/env python3
"""Verify frozen target-blind inputs for the EPSR sensitivity matrix."""

from __future__ import annotations

import hashlib
import json
import math
import tomllib
from pathlib import Path

import numpy as np


case_root = Path(__file__).resolve().parents[1]
result_root = case_root / "results/epsr-control-start-sensitivity"
input_root = result_root / "inputs"
manifest_path = input_root / "input-manifest.json"
manifest = json.loads(manifest_path.read_text())
expected = tomllib.loads(
    (case_root / "expected/epsr-control-start-sensitivity-inputs.toml").read_text()
)
protocol_path = case_root / "expected/epsr-control-start-sensitivity.toml"
protocol = tomllib.loads(protocol_path.read_text())


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(label: str, actual: float, wanted: float, tolerance: float = 5e-9):
    if not math.isclose(actual, wanted, rel_tol=0.0, abs_tol=tolerance):
        raise SystemExit(f"{label}: {actual} != {wanted}")


def read_lammps(path: Path):
    lines = path.read_text().splitlines()
    count = int(next(line.split()[0] for line in lines if line.strip().endswith("atoms")))
    bounds = {}
    for line in lines:
        fields = line.split()
        if len(fields) >= 4 and fields[-2:] in (
            ["xlo", "xhi"],
            ["ylo", "yhi"],
            ["zlo", "zhi"],
        ):
            bounds[fields[-2][0]] = float(fields[1]) - float(fields[0])
    header = next(line for line in lines if line.startswith("Atoms"))
    offset = 3 if "charge" in header else 2
    start = lines.index(header) + 2
    rows = [line.split() for line in lines[start : start + count]]
    types = np.asarray([int(row[1]) for row in rows], dtype=int)
    positions = np.asarray(
        [[float(value) for value in row[offset : offset + 3]] for row in rows]
    )
    box = np.asarray([bounds[axis] for axis in "xyz"])
    return types, positions, box


def periodic_rms(first, second, box) -> float:
    difference = first - second
    difference -= box * np.rint(difference / box)
    return float(np.sqrt(np.mean(np.sum(difference * difference, axis=1))))


def minimum_distance(positions, box) -> float:
    minimum = math.inf
    for first in range(0, len(positions), 250):
        block = positions[first : first + 250]
        differences = block[:, None, :] - positions[None, :, :]
        differences -= box * np.rint(differences / box)
        distances = np.sqrt(np.sum(differences * differences, axis=2))
        distances[np.arange(len(block)), np.arange(first, first + len(block))] = np.inf
        minimum = min(minimum, float(np.min(distances)))
    return minimum


def potential_rows(path: Path):
    return np.asarray(
        [
            [float(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
    )


if sha256(protocol_path) != expected["protocol_sha256"]:
    raise SystemExit("sensitivity protocol changed after input generation")
if manifest["protocol_sha256"] != expected["protocol_sha256"]:
    raise SystemExit("input manifest carries the wrong protocol hash")
if manifest["binary_sha256"] != expected["binary_sha256"]:
    raise SystemExit("input manifest carries the wrong generator-binary hash")

baseline_tables = {}
for scale_label, scale in (("0p5", 0.5), ("1p0", 1.0), ("1p5", 1.5)):
    scale_root = input_root / "reference-potentials" / f"scale-{scale_label}"
    for pair in ("Si-Si", "Si-O", "O-O"):
        path = scale_root / f"epsr-reference-{pair}.dat"
        wanted_hash = expected["reference_potentials"][f"scale_{scale_label}"][pair]
        if sha256(path) != wanted_hash:
            raise SystemExit(f"reference-potential hash changed: {scale_label}/{pair}")
        rows = potential_rows(path)
        if scale_label == "1p0":
            baseline_tables[pair] = rows
        else:
            baseline = potential_rows(
                input_root
                / "reference-potentials/scale-1p0"
                / f"epsr-reference-{pair}.dat"
            )
            if not np.array_equal(rows[:, 0], baseline[:, 0]):
                raise SystemExit(f"reference grid changed: {scale_label}/{pair}")
            if not np.allclose(rows[:, 1], scale * baseline[:, 1], rtol=1e-14, atol=1e-14):
                raise SystemExit(f"reference energies not scaled exactly: {scale_label}/{pair}")

generated_names = protocol["starts"]["names"][1:]
for case_name, wanted_case in expected["cases"].items():
    original_path = result_root.parent / "cross-recovery" / case_name / "cross-start.data"
    if sha256(original_path) != wanted_case["original_sha256"]:
        raise SystemExit(f"{case_name}: original start changed")
    original_types, original_positions, original_box = read_lammps(original_path)
    generated = {}
    for start_name in generated_names:
        wanted = wanted_case[start_name]
        run = input_root / "starts" / case_name / start_name
        start_path = run / "start.data"
        if sha256(start_path) != wanted["sha256"]:
            raise SystemExit(f"{case_name}/{start_name}: generated-start hash changed")
        types, positions, box = read_lammps(start_path)
        if len(types) != expected["atom_count"]:
            raise SystemExit(f"{case_name}/{start_name}: wrong atom count")
        if int(np.sum(types == 1)) != expected["si_count"] or int(
            np.sum(types == 2)
        ) != expected["o_count"]:
            raise SystemExit(f"{case_name}/{start_name}: wrong composition")
        if not np.array_equal(types, original_types):
            raise SystemExit(f"{case_name}/{start_name}: atom identities reordered")
        if not np.allclose(box, original_box, rtol=0.0, atol=5e-7):
            raise SystemExit(f"{case_name}/{start_name}: box changed")
        close(
            f"{case_name}/{start_name}/minimum distance",
            minimum_distance(positions, box),
            wanted["minimum_distance_a"],
        )
        close(
            f"{case_name}/{start_name}/displacement from original",
            periodic_rms(positions, original_positions, box),
            wanted["from_original_rms_a"],
        )
        config = tomllib.loads((run / "run.toml").read_text())
        if "constraints" in config:
            raise SystemExit(f"{case_name}/{start_name}: constraints unexpectedly active")
        if config["epsr"]["mode"] != "pure" or config["epsr"]["feedback"] != 0.0:
            raise SystemExit(f"{case_name}/{start_name}: generator was not target-blind")
        if config["epsr"]["iterations"] != 1:
            raise SystemExit(f"{case_name}/{start_name}: wrong generator iterations")
        if config["rmc"]["max_moves"] != protocol["starts"]["generator_moves"]:
            raise SystemExit(f"{case_name}/{start_name}: wrong generator move count")
        generated[start_name] = positions
    close(
        f"{case_name}/difference between generated starts",
        periodic_rms(
            generated[generated_names[0]], generated[generated_names[1]], original_box
        ),
        wanted_case["between_generated_starts_rms_a"],
    )

print(
    "EPSR control/start inputs: PASS "
    "(hashes, exact scaling, target blindness, composition, boxes, and geometry)"
)
