#!/usr/bin/env python3
"""Prepare scaled reference tables and target-blind common sensitivity starts."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import tomllib
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--binary", type=Path)
parser.add_argument("--force", action="store_true")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
repo_root = case_root.parents[1]
binary = args.binary or repo_root / "target/release/rsmith"
if not binary.is_file():
    raise SystemExit(f"rsmith binary not found: {binary}")
protocol_path = case_root / "expected/epsr-control-start-sensitivity.toml"
protocol = tomllib.loads(protocol_path.read_text())
fixture_root = case_root / "results/cross-recovery"
result_root = case_root / "results/epsr-control-start-sensitivity"
input_root = result_root / "inputs"
environment = {
    **os.environ,
    "RAYON_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_one(text: str, pattern: str, replacement: str) -> str:
    rewritten, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise ValueError(f"expected one configuration match for {pattern!r}")
    return rewritten


def scale_table(source: Path, target: Path, scale: float) -> None:
    if scale == 1.0:
        shutil.copy2(source, target)
        return
    lines = []
    for line in source.read_text().splitlines():
        fields = line.split()
        if not fields or line.lstrip().startswith("#"):
            lines.append(line)
            continue
        if len(fields) != 2:
            raise ValueError(f"unexpected potential row in {source}: {line!r}")
        lines.append(f"{float(fields[0]):.8f} {float(fields[1]) * scale:.17g}")
    target.write_text("\n".join(lines) + "\n")


def write_lammps_from_xyz(source: Path, target: Path, *, charge_style: bool) -> None:
    lines = source.read_text().splitlines()
    count = int(lines[0])
    match = re.search(r'Lattice="([^"]+)"', lines[1])
    if match is None:
        raise ValueError(f"missing lattice in {source}")
    lattice = [float(value) for value in match.group(1).split()]
    box = [lattice[0], lattice[4], lattice[8]]
    atoms = []
    types = {"Si": 1, "O": 2}
    for atom_id, line in enumerate(lines[2 : 2 + count], start=1):
        fields = line.split()
        atoms.append(
            (atom_id, types[fields[0]], *(float(value) for value in fields[1:4]))
        )
    output = [
        "Target-blind reference-potential-equilibrated SiO2 start",
        "",
        f"{count} atoms",
        "2 atom types",
        "",
        f"0.0 {box[0]:.12f} xlo xhi",
        f"0.0 {box[1]:.12f} ylo yhi",
        f"0.0 {box[2]:.12f} zlo zhi",
        "",
        "Masses",
        "",
        "1 28.0855",
        "2 15.999",
        "",
        "Atoms # charge" if charge_style else "Atoms # atomic",
        "",
    ]
    if charge_style:
        output.extend(
            f"{atom_id} {atom_type} 0.0 {x:.12f} {y:.12f} {z:.12f}"
            for atom_id, atom_type, x, y, z in atoms
        )
    else:
        output.extend(
            f"{atom_id} {atom_type} {x:.12f} {y:.12f} {z:.12f}"
            for atom_id, atom_type, x, y, z in atoms
        )
    target.write_text("\n".join(output) + "\n")


reference_sources = {
    pair: case_root
    / "results/reference-potential-smoke"
    / f"epsr-reference-{pair}.dat"
    for pair in ("Si-Si", "Si-O", "O-O")
}
reference_scales = sorted({float(arm["reference_scale"]) for arm in protocol["arms"]})
reference_manifest = {}
for scale in reference_scales:
    label = str(scale).replace(".", "p")
    scale_root = input_root / "reference-potentials" / f"scale-{label}"
    scale_root.mkdir(parents=True, exist_ok=True)
    reference_manifest[label] = {}
    for pair, source in reference_sources.items():
        target = scale_root / source.name
        scale_table(source, target, scale)
        reference_manifest[label][pair] = {
            "path": str(target),
            "sha256": sha256(target),
        }

start_names = protocol["starts"]["names"][1:]
generator_seeds = [int(seed) for seed in protocol["starts"]["generator_seeds"]]
if len(start_names) != len(generator_seeds):
    raise SystemExit("start names and generator seeds differ in length")
generator_moves = int(protocol["starts"]["generator_moves"])
start_manifest = {}
scale_one_root = input_root / "reference-potentials/scale-1p0"
for case in sorted(fixture_root.glob("target-*_*")):
    source_config = case / "rsmith-epsr-joint/run.toml"
    start_manifest[case.name] = {}
    for start_name, seed in zip(start_names, generator_seeds):
        run = input_root / "starts" / case.name / start_name
        if run.exists():
            if not args.force:
                raise SystemExit(f"output exists: {run} (pass --force)")
            shutil.rmtree(run)
        run.mkdir(parents=True)
        config = source_config.read_text()
        config, count = re.subn(
            r"\n\[constraints\]\n.*?(?=\n\[potential\])",
            "\n",
            config,
            count=1,
            flags=re.DOTALL,
        )
        if count != 1:
            raise ValueError("expected one [constraints] section")
        config = replace_one(config, r"^max_moves = \d+$", f"max_moves = {generator_moves}")
        config = replace_one(config, r"^seed = \d+$", f"seed = {seed}")
        config = replace_one(
            config, r"^print_every = \d+$", f"print_every = {generator_moves}"
        )
        config = replace_one(config, r"^iterations = \d+$", "iterations = 1")
        config = replace_one(config, r"^feedback = [0-9.eE+-]+$", "feedback = 0.0")
        config = replace_one(
            config,
            r"^moves_per_iteration = \d+$",
            f"moves_per_iteration = {generator_moves}",
        )
        for pair in ("Si-Si", "Si-O", "O-O"):
            config = replace_one(
                config,
                rf'^file = ".*epsr-reference-{pair}\.dat"$',
                f'file = "{scale_one_root / f"epsr-reference-{pair}.dat"}"',
            )
        config_path = run / "run.toml"
        config_path.write_text(config)
        completed = subprocess.run(
            [str(binary), str(config_path), "--output-dir", str(run), "--quiet"],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env=environment,
            check=False,
        )
        (run / "driver.log").write_text(completed.stdout)
        if completed.returncode != 0:
            raise RuntimeError(f"start generation failed in {run}; see driver.log")
        start_path = run / "start.data"
        rsmith_start_path = run / "start-charge.data"
        write_lammps_from_xyz(
            run / "refined.xyz", start_path, charge_style=False
        )
        write_lammps_from_xyz(
            run / "refined.xyz", rsmith_start_path, charge_style=True
        )
        start_manifest[case.name][start_name] = {
            "generator_seed": seed,
            "generator_moves": generator_moves,
            "path": str(start_path),
            "sha256": sha256(start_path),
            "rsmith_path": str(rsmith_start_path),
            "rsmith_sha256": sha256(rsmith_start_path),
            "coordinate_equivalence": "exact atom ids, types, box, and coordinates; only the LAMMPS atom-style representation differs",
        }

manifest = {
    "status": "epsr_control_start_sensitivity_inputs_prepared",
    "protocol": str(protocol_path),
    "protocol_sha256": sha256(protocol_path),
    "binary": str(binary),
    "binary_sha256": sha256(binary),
    "reference_potentials": reference_manifest,
    "starts": start_manifest,
}
manifest_path = input_root / "input-manifest.json"
manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
print(manifest_path)
