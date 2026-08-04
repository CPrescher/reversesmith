#!/usr/bin/env python3
"""Run stock EPSR26 for the frozen five-seed 10 GPa comparison."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import time
import tomllib
from pathlib import Path


def numeric_rows(path: Path, minimum_columns: int = 2):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        try:
            values = tuple(float(value) for value in line.split())
        except ValueError:
            continue
        if len(values) >= minimum_columns:
            rows.append(values)
    return rows


def replace_setting(text: str, name: str, value: str):
    pattern = re.compile(rf"^(\s*{re.escape(name)}\s+)([^\r\n]*)(\r?$)", re.MULTILINE)
    match = pattern.search(text)
    if match is None:
        raise ValueError(f"setting {name!r} not found in EPSR input")
    old = match.group(2)
    separator = "               "
    comment = separator + old.split(separator, 1)[1] if separator in old else ""
    replacement = match.group(1) + value + comment + match.group(3)
    return text[: match.start()] + replacement + text[match.end() :]


def read_lammps(path: Path):
    bounds, atoms = {}, []
    atom_style, in_atoms = None, False
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 4 and fields[-2:] in (
            ["xlo", "xhi"],
            ["ylo", "yhi"],
            ["zlo", "zhi"],
        ):
            bounds[fields[-2][0]] = (float(fields[0]), float(fields[1]))
            continue
        if line.strip().startswith("Atoms"):
            atom_style = "charge" if "charge" in line else "atomic"
            in_atoms = True
            continue
        if not in_atoms or not line.strip():
            continue
        if line.lstrip()[0].isalpha():
            in_atoms = False
            continue
        fields = line.split()
        offset = 3 if atom_style == "charge" else 2
        atoms.append(
            (
                int(fields[0]),
                int(fields[1]),
                tuple(float(value) for value in fields[offset : offset + 3]),
            )
        )
    atoms.sort()
    return atoms, tuple(bounds[axis][1] - bounds[axis][0] for axis in "xyz")


def write_ato(path: Path, structure: Path, template: Path, seed: int):
    atoms, box = read_lammps(structure)
    if len(atoms) != 3000 or max(box) - min(box) > 1.0e-9:
        raise ValueError("EPSR comparison requires 3000 atoms in a cubic box")
    template_lines = template.read_text().splitlines()
    template_count = int(template_lines[0].split()[0])
    tail = template_lines[2 + 5 * template_count :]
    for index, line in enumerate(tail):
        fields = line.split()
        if len(fields) != 35:
            continue
        try:
            [int(value) for value in fields]
        except ValueError:
            continue
        tail[index] = " " + " ".join(str(value) for value in [0, -seed, 0, *([0] * 32)])
        break
    with path.open("w") as stream:
        stream.write(f"{len(atoms):8d} {box[0]:15.8E} {300.0:15.8E}\n")
        stream.write(template_lines[1] + "\n")
        for atom_id, atom_type, position in atoms:
            centered = [position[axis] - 0.5 * box[axis] for axis in range(3)]
            stream.write(
                " 1 "
                + " ".join(f"{value:15.8E}" for value in centered)
                + "  0.00000000E+00  0.00000000E+00  0.00000000E+00"
                + f"  F      0    1.00000 {atom_id:5d}\n"
            )
            stream.write(f" {('Si' if atom_type == 1 else 'O'):<8s} 1      0\n")
            stream.write("     0.00000000     0.00000000     0.00000000\n    0\n    0\n")
        stream.write("\n".join(tail) + "\n")


def read_iq(path: Path):
    return [(row[0], row[1], row[2]) for row in numeric_rows(path, 3)]


def xray_form_factor(values, q):
    s2 = (q / (4.0 * math.pi)) ** 2
    return values[0] + sum(
        amplitude * math.exp(-decay * s2)
        for amplitude, decay in zip(values[1::2], values[2::2])
    )


def provisional_targets(input_root: Path, source: Path):
    neutron = read_iq(input_root / "target-neutron-iq.dat")
    xray = read_iq(input_root / "target-xray-iq.dat")
    neutron_weights = [
        row[2]
        for row in numeric_rows(source / "DTBsilica.NWTStot.wts", 3)
        if len(row) == 3
    ]
    neutron_factor = sum(neutron_weights)
    xparams = [
        row for row in numeric_rows(source / "DTBsilica.XWTS.wts", 11) if len(row) == 11
    ]
    if len(xparams) != 2:
        raise ValueError("expected Si and O EPSR X-ray parameter rows")
    neutron_native = [
        (q, value * neutron_factor, sigma * neutron_factor)
        for q, value, sigma in neutron
    ]
    xray_native = []
    for q, value, sigma in xray:
        f_si = xray_form_factor(xparams[0], q)
        f_o = xray_form_factor(xparams[1], q)
        mean_f = (f_si + 2.0 * f_o) / 3.0
        mean_f2 = (f_si * f_si + 2.0 * f_o * f_o) / 3.0
        factor = mean_f2 / (mean_f * mean_f)
        xray_native.append((q, value / factor, sigma / factor))
    return neutron_native, xray_native


def write_epsr_data(path: Path, title: str, rows):
    step = rows[-1][0] - rows[-2][0]
    extended = [*rows, (rows[-1][0] + step, rows[-1][1], rows[-1][2])]
    path.write_text(
        f"# {title}\n"
        + "".join(
            f"{q:.12g} {value:.12g} {sigma:.12g}\n"
            for q, value, sigma in extended
        )
    )


def configure_input(path: Path, ato_name: str, moves: int, refinements: int):
    text = path.read_text()
    settings = (
        ("feedback", "0.90000000"),
        ("potfac", "0.0 0.0" if moves == 0 else "1.0 1.0"),
        ("num_threds", "1"),
        ("nq", "500"),
        ("qstep", "0.05"),
        ("ireset", "1" if moves == 0 else "2"),
        ("iinit", "1"),
        ("ntimes", "0" if moves == 0 else str(max(1, moves // 3000))),
        ("niter", str(refinements)),
        ("nsumt", "0"),
        ("rho", "0.0"),
        ("cellst", "0.02"),
        ("rmaxgr", "16.0"),
        ("ngrsamples", "3000"),
        ("fwhm", "0.0"),
        ("fwhmq", "0.0"),
        ("nsmoop", "0"),
        ("fnameato", ato_name),
        ("qmin", "0.5 0.0"),
    )
    for name, value in settings:
        text = replace_setting(text, name, value)
    text = text.replace("SilicaGlassRT.mint01", "target-neutron.dat")
    text = text.replace("SiO2XRD.int01", "target-xray.dat")
    text = re.sub(r"(\bnrtype\s+)5", r"\g<1>6", text)
    path.write_text(text)


def prepare_run(
    source: Path,
    run: Path,
    structure: Path,
    moves: int,
    refinements: int,
    seed: int,
    targets,
):
    if run.exists():
        shutil.rmtree(run)
    shutil.copytree(source, run)
    for output in run.glob("DTBsilica.EPSR.*"):
        if output.suffix not in {".inp", ".inpa"}:
            output.unlink()
    write_ato(run / "Cross.ato", structure, source / "DTBsilica.ato", seed)
    write_epsr_data(run / "target-neutron.dat", "EPSR-native neutron i(Q)", targets[0])
    write_epsr_data(run / "target-xray.dat", "EPSR-native X-ray i(Q)", targets[1])
    configure_input(run / "DTBsilica.EPSR.inp", "Cross.ato", moves, refinements)


def run_epsr(binary: Path, run: Path):
    started = time.perf_counter()
    completed = subprocess.run(
        [str(binary), str(run) + "/", "epsr", "DTBsilica"],
        cwd=run,
        text=True,
        input="",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    wall = time.perf_counter() - started
    (run / "driver.log").write_text(completed.stdout)
    (run / "wall-seconds.txt").write_text(f"{wall:.9f}\n")
    if completed.returncode != 0 or not (run / "DTBsilica.EPSR.v01").is_file():
        raise RuntimeError(f"native EPSR failed in {run}; see driver.log")
    return wall


def calculated_targets(v01: Path, input_targets, sigma: float = 0.02):
    rows = [row for row in numeric_rows(v01, 5) if 0.5 <= row[0] < 25.0]
    if len(rows) != len(input_targets[0]) or len(rows) != len(input_targets[1]):
        raise ValueError("EPSR residual and target grids differ")
    return (
        [
            (row[0], input_targets[0][index][1] - row[1], sigma)
            for index, row in enumerate(rows)
        ],
        [
            (row[0], input_targets[1][index][1] - row[3], sigma)
            for index, row in enumerate(rows)
        ],
    )


def rms(first, second):
    residual = [a[1] - b[1] for a, b in zip(first, second)]
    return math.sqrt(sum(value * value for value in residual) / len(residual))


parser = argparse.ArgumentParser()
parser.add_argument("--only-missing", action="store_true")
parser.add_argument("--force", action="store_true")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
comparison = tomllib.loads((case_root / "expected/ten-gpa-comparison.toml").read_text())
pilot_root = case_root / "results/ten-gpa-pilot"
root = case_root / "results/ten-gpa-comparison"
source = case_root.parent / "epsr-silica/reference/local/upstream"
record_path = case_root.parent / "epsr-silica/reference/local/IMPORT.txt"
if not source.is_dir() or not record_path.is_file():
    raise SystemExit("import the accepted EPSR26 local-testing tutorial first")
record = dict(
    line.split("=", 1) for line in record_path.read_text().splitlines() if "=" in line
)
binary = Path(record["epsr_root"]) / "bin/epsr"
if not binary.is_file():
    raise SystemExit(f"EPSR executable not found: {binary}")
root.mkdir(parents=True, exist_ok=True)

target_root = root / "native-targets"
target_root.mkdir(exist_ok=True)
target_files = (target_root / "target-neutron.dat", target_root / "target-xray.dat")
forward_summary_path = root / "native-forward-summary.json"
if args.force or not all(path.is_file() for path in target_files):
    provisional = provisional_targets(pilot_root / "inputs", source)
    first = root / "native-forward/provisional"
    prepare_run(
        source,
        first,
        pilot_root / "inputs/hidden-target.data",
        0,
        1,
        int(comparison["design"]["primary_seed"]),
        provisional,
    )
    first_wall = run_epsr(binary, first)
    native_targets = calculated_targets(first / "DTBsilica.EPSR.v01", provisional)
    write_epsr_data(target_files[0], "EPSR-native hidden neutron i(Q)", native_targets[0])
    write_epsr_data(target_files[1], "EPSR-native hidden X-ray i(Q)", native_targets[1])
    replay = root / "native-forward/replay"
    prepare_run(
        source,
        replay,
        pilot_root / "inputs/hidden-target.data",
        0,
        1,
        int(comparison["design"]["primary_seed"]),
        native_targets,
    )
    replay_wall = run_epsr(binary, replay)
    checked = calculated_targets(replay / "DTBsilica.EPSR.v01", native_targets)
    forward_summary = {
        "first_wall_seconds": first_wall,
        "hidden_replay_neutron_rms": rms(checked[0], native_targets[0]),
        "hidden_replay_wall_seconds": replay_wall,
        "hidden_replay_xray_rms": rms(checked[1], native_targets[1]),
        "status": "epsr_native_hidden_forward_generated_and_replayed",
    }
    forward_summary_path.write_text(
        json.dumps(forward_summary, indent=2, sort_keys=True) + "\n"
    )
else:
    native_targets = tuple(
        [row for row in read_iq(path) if 0.5 <= row[0] < 25.0]
        for path in target_files
    )

summary_path = root / "native-epsr-run-summary.json"
summary = (
    json.loads(summary_path.read_text())
    if summary_path.is_file()
    else {"status": "hp_sio2_10gpa_native_epsr_runs", "runs": {}}
)
primary = int(comparison["design"]["primary_seed"])
moves_per_refinement = int(comparison["native_epsr26"]["moves_per_refinement"])
for seed in comparison["design"]["seeds"]:
    refinements_list = (
        comparison["native_epsr26"]["checkpoints_refinements"]
        if int(seed) == primary
        else [int(comparison["design"]["endpoint_moves"]) // moves_per_refinement]
    )
    for refinements in refinements_list:
        moves = int(refinements) * moves_per_refinement
        run = root / "runs/native-epsr26" / f"seed-{int(seed)}" / f"moves-{moves:06d}"
        if args.only_missing and (run / "DTBsilica.EPSR.v01").is_file():
            continue
        if run.exists() and not args.force:
            raise SystemExit(f"output exists: {run} (pass --force or --only-missing)")
        prepare_run(
            source,
            run,
            pilot_root / "inputs/cross-start.data",
            moves_per_refinement,
            int(refinements),
            int(seed),
            native_targets,
        )
        wall = run_epsr(binary, run)
        name = str(run.relative_to(root / "runs"))
        summary["runs"][name] = {
            "attempted_moves": moves,
            "refinements": int(refinements),
            "status": "completed",
            "wall_seconds": wall,
        }
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        print(json.dumps({name: summary["runs"][name]}, indent=2), flush=True)
