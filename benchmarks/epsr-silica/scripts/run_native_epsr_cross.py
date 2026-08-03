#!/usr/bin/env python3
"""Run native EPSR26 hidden-forward gates and symmetric cross-start smokes."""

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


def replace_setting(text: str, name: str, value: str) -> str:
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
    lines = path.read_text().splitlines()
    bounds = {}
    atoms = []
    atom_style = None
    in_atoms = False
    for line in lines:
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
    box = tuple(bounds[axis][1] - bounds[axis][0] for axis in "xyz")
    return atoms, box


def write_ato(path: Path, structure: Path, template: Path, seed: int):
    atoms, box = read_lammps(structure)
    if len(atoms) != 3000 or max(box) - min(box) > 1.0e-9:
        raise ValueError("EPSR cross benchmark requires 3000 atoms in a cubic box")
    template_lines = template.read_text().splitlines()
    template_count = int(template_lines[0].split()[0])
    tail = template_lines[2 + 5 * template_count :]
    # Replace the copied RNG state deterministically if the 35-integer record is found.
    for index, line in enumerate(tail):
        fields = line.split()
        if len(fields) == 35:
            try:
                [int(value) for value in fields]
            except ValueError:
                continue
            tail[index] = " " + " ".join(
                str(value) for value in [0, -seed, 0, *([0] * 32)]
            )
            break
    with path.open("w") as stream:
        stream.write(f"{len(atoms):8d} {box[0]:15.8E} {300.0:15.8E}\n")
        stream.write(template_lines[1] + "\n")
        for atom_id, atom_type, position in atoms:
            centered = [position[axis] - 0.5 * box[axis] for axis in range(3)]
            stream.write(
                " 1 "
                + " ".join(f"{value:15.8E}" for value in centered)
                + "  0.00000000E+00  0.00000000E+00  0.00000000E+00  F      0    1.00000 "
                + f"{atom_id:5d}\n"
            )
            stream.write(f" {('Si' if atom_type == 1 else 'O'):<8s} 1      0\n")
            stream.write(
                "     0.00000000     0.00000000     0.00000000\n    0\n    0\n"
            )
        stream.write("\n".join(tail) + "\n")


def read_iq(path: Path):
    return [(row[0], row[1], row[2]) for row in numeric_rows(path, 3)]


def xray_parameters(path: Path):
    rows = [row for row in numeric_rows(path, 11) if len(row) == 11]
    if len(rows) != 2:
        raise ValueError("expected Si and O EPSR X-ray parameter rows")
    return rows


def xray_form_factor(values, q):
    s2 = (q / (4.0 * math.pi)) ** 2
    return values[0] + sum(
        amplitude * math.exp(-decay * s2)
        for amplitude, decay in zip(values[1::2], values[2::2])
    )


def write_epsr_data(path: Path, title: str, rows):
    # EPSR nrtype 6 discards its final record; append one sentinel grid point.
    if len(rows) < 2:
        raise ValueError("EPSR target requires at least two points")
    step = rows[-1][0] - rows[-2][0]
    extended = [*rows, (rows[-1][0] + step, rows[-1][1], rows[-1][2])]
    path.write_text(
        f"# {title}\n"
        + "".join(
            f"{q:.12g} {value:.12g} {sigma:.12g}\n" for q, value, sigma in extended
        )
    )


def set_pcof_rminex(path: Path, values: dict[tuple[str, str], float]) -> None:
    lines = path.read_text().splitlines()
    found = set()
    for index, line in enumerate(lines[:-1]):
        fields = line.split()
        if len(fields) < 2:
            continue
        pair = (fields[0], fields[1])
        if pair not in values:
            continue
        controls = lines[index + 1].split()
        if len(controls) < 3:
            raise ValueError(f"invalid pcof pair-control record after {line!r}")
        controls[0] = f"{values[pair]:.8f}"
        lines[index + 1] = "   " + "   ".join(controls)
        found.add(pair)
    if found != set(values):
        raise ValueError(f"missing pcof pairs: {set(values) - found}")
    path.write_text("\n".join(lines) + "\n")


def provisional_targets(case: Path, source: Path):
    neutron = read_iq(case / "target-neutron-iq.dat")
    xray = read_iq(case / "target-xray-iq.dat")
    neutron_weights = [
        row[2]
        for row in numeric_rows(source / "DTBsilica.NWTStot.wts", 3)
        if len(row) == 3
    ]
    neutron_factor = sum(neutron_weights)
    xparams = xray_parameters(source / "DTBsilica.XWTS.wts")
    neutron_native = [
        (q, value * neutron_factor, sigma * neutron_factor)
        for q, value, sigma in neutron
    ]
    xray_native = []
    for q, value, sigma in xray:
        f_si, f_o = xray_form_factor(xparams[0], q), xray_form_factor(xparams[1], q)
        mean_f = (f_si + 2.0 * f_o) / 3.0
        mean_f2 = (f_si * f_si + 2.0 * f_o * f_o) / 3.0
        factor = mean_f2 / (mean_f * mean_f)
        xray_native.append((q, value / factor, sigma / factor))
    return neutron_native, xray_native


def configure_input(
    path: Path,
    ato_name: str,
    moves: int,
    native_target_names=("target-neutron.dat", "target-xray.dat"),
    refinements: int = 1,
    feedback: float = 0.9,
):
    text = path.read_text()
    settings = (
        ("feedback", f"{feedback:.8f}"),
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
    text = text.replace("SilicaGlassRT.mint01", native_target_names[0])
    text = text.replace("SiO2XRD.int01", native_target_names[1])
    text = re.sub(r"(\bnrtype\s+)5", r"\g<1>6", text)
    path.write_text(text)


def prepare_run(
    source: Path,
    case: Path,
    name: str,
    structure: Path,
    moves: int,
    seed: int,
    targets,
    force: bool,
    refinements: int = 1,
    rminex: dict[tuple[str, str], float] | None = None,
    feedback: float = 0.9,
    reference_scale: float = 1.0,
):
    run_dir = case / name
    if run_dir.exists():
        if not force:
            raise FileExistsError(f"output exists: {run_dir} (pass --force)")
        shutil.rmtree(run_dir)
    shutil.copytree(source, run_dir)
    if rminex is not None:
        set_pcof_rminex(run_dir / "DTBsilica.pcof", rminex)
    pcof_path = run_dir / "DTBsilica.pcof"
    pcof_path.write_text(
        replace_setting(
            pcof_path.read_text(), "refpotfac", f"{reference_scale:.8f}"
        )
    )
    for output in run_dir.glob("DTBsilica.EPSR.*"):
        if output.suffix not in {".inp", ".inpa"}:
            output.unlink()
    write_ato(run_dir / "Cross.ato", structure, source / "DTBsilica.ato", seed)
    write_epsr_data(
        run_dir / "target-neutron.dat",
        "native-format synthetic neutron i(Q)",
        targets[0],
    )
    write_epsr_data(
        run_dir / "target-xray.dat", "native-format synthetic X-ray i(Q)", targets[1]
    )
    configure_input(
        run_dir / "DTBsilica.EPSR.inp",
        "Cross.ato",
        moves,
        refinements=refinements,
        feedback=feedback,
    )
    return run_dir


def run_epsr(binary: Path, run_dir: Path):
    started = time.perf_counter()
    completed = subprocess.run(
        [str(binary), str(run_dir) + "/", "epsr", "DTBsilica"],
        cwd=run_dir,
        text=True,
        input="",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    wall = time.perf_counter() - started
    (run_dir / "driver.log").write_text(completed.stdout)
    (run_dir / "wall-seconds.txt").write_text(f"{wall:.9f}\n")
    if completed.returncode != 0:
        raise RuntimeError(f"native EPSR failed in {run_dir}; see driver.log")
    required = run_dir / "DTBsilica.EPSR.v01"
    if not required.is_file():
        raise RuntimeError(f"native EPSR did not create {required}")
    return wall


def calculated_targets(v01: Path, input_targets, sigma=0.02):
    rows = [row for row in numeric_rows(v01, 5) if 0.5 <= row[0] < 25.0]
    if len(rows) != len(input_targets[0]) or len(rows) != len(input_targets[1]):
        raise ValueError("EPSR residual and input target grids differ")
    # EPSR .v01 columns are data-minus-model residuals, despite the dataset
    # names in the header.  Recover the forward calculation as input-residual.
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
    if len(first) != len(second):
        raise ValueError(f"target lengths differ: {len(first)} != {len(second)}")
    residual = [a[1] - b[1] for a, b in zip(first, second)]
    return math.sqrt(sum(value * value for value in residual) / len(residual))


parser = argparse.ArgumentParser()
parser.add_argument("--epsr-root", type=Path)
parser.add_argument("--moves", type=int, default=6000)
parser.add_argument("--seed", type=int, default=20260802)
parser.add_argument("--force", action="store_true")
parser.add_argument("--zero-only", action="store_true")
parser.add_argument("--ensemble", action="store_true")
parser.add_argument("--convergence-pilot", action="store_true")
parser.add_argument("--convergence-ensemble", action="store_true")
parser.add_argument("--control-start-sensitivity", action="store_true")
parser.add_argument("--sensitivity-arm", action="append")
parser.add_argument("--convergence-seed", type=int, action="append")
parser.add_argument("--checkpoint", type=int, action="append")
parser.add_argument(
    "--native-rminex-control",
    action="store_true",
    help="test EPSR pcof rminex values; these are not assumed to be hard constraints",
)
parser.add_argument("--only-missing", action="store_true")
args = parser.parse_args()
if args.moves <= 0:
    raise SystemExit("--moves must be positive")

case_root = Path(__file__).resolve().parents[1]
source = case_root / "reference/local/upstream"
record_path = case_root / "reference/local/IMPORT.txt"
if not source.is_dir() or not record_path.is_file():
    raise SystemExit("run import_local_reference.sh --accept-local-testing-terms first")
record = dict(
    line.split("=", 1) for line in record_path.read_text().splitlines() if "=" in line
)
epsr_root = args.epsr_root or Path(record["epsr_root"])
binary = epsr_root / "bin/epsr"
if not binary.is_file():
    raise SystemExit(f"EPSR executable not found: {binary}")

fixture_root = case_root / "results/cross-recovery"
if args.control_start_sensitivity:
    if args.native_rminex_control:
        raise SystemExit("the rminex control is not part of the sensitivity matrix")
    protocol_name = "epsr-control-start-sensitivity.toml"
    protocol = tomllib.loads((case_root / "expected" / protocol_name).read_text())
    result_root = case_root / "results/epsr-control-start-sensitivity"
    inputs_root = result_root / "inputs"
    if not (inputs_root / "input-manifest.json").is_file():
        raise SystemExit("run prepare_epsr_control_start_sensitivity.py first")
    arms = {arm["name"]: arm for arm in protocol["arms"]}
    selected_arms = args.sensitivity_arm or list(arms)
    unknown_arms = set(selected_arms) - set(arms)
    if unknown_arms:
        raise SystemExit(f"unknown sensitivity arms: {sorted(unknown_arms)}")
    configured_seeds = [int(seed) for seed in protocol["design"]["refinement_seeds"]]
    selected_seeds = args.convergence_seed or configured_seeds
    unknown_seeds = set(selected_seeds) - set(configured_seeds)
    if unknown_seeds:
        raise SystemExit(f"seeds not present in {protocol_name}: {sorted(unknown_seeds)}")
    moves_per_refinement = int(protocol["design"]["moves_per_refinement"])
    checkpoints = args.checkpoint or protocol["methods"]["native_epsr26"][
        "checkpoints"
    ]
    for arm_name in selected_arms:
        arm = arms[arm_name]
        for seed in selected_seeds:
            summary_path = (
                result_root
                / "run-summaries"
                / f"{arm_name}-native-epsr26-seed-{seed}.json"
            )
            summary = {
                "program": "EPSR26",
                "binary": str(binary),
                "arm": arm_name,
                "feedback": float(arm["feedback"]),
                "reference_scale": float(arm["reference_scale"]),
                "start": arm["start"],
                "threads": 1,
                "seed": seed,
                "cases": {},
            }
            for case in sorted(fixture_root.glob("target-*_*")):
                targets = []
                for filename in (
                    "epsr-native-target-neutron.dat",
                    "epsr-native-target-xray.dat",
                ):
                    rows = read_iq(case / filename)
                    targets.append([row for row in rows if 0.5 <= row[0] < 25.0])
                structure = (
                    case / "cross-start.data"
                    if arm["start"] == "original"
                    else inputs_root
                    / "starts"
                    / case.name
                    / arm["start"]
                    / "start.data"
                )
                if not structure.is_file():
                    raise SystemExit(f"missing sensitivity start: {structure}")
                prefix_root = (
                    result_root
                    / case.name
                    / arm_name
                    / "native-epsr26"
                    / f"seed-{seed}"
                )
                prefix_root.mkdir(parents=True, exist_ok=True)
                case_summary = summary["cases"].setdefault(case.name, {})
                for refinements in checkpoints:
                    name = f"iter-{int(refinements):03d}"
                    run_dir = prefix_root / name
                    if args.only_missing and (
                        run_dir / "DTBsilica.EPSR.v01"
                    ).is_file():
                        continue
                    run_dir = prepare_run(
                        source,
                        prefix_root,
                        name,
                        structure,
                        moves_per_refinement,
                        seed,
                        tuple(targets),
                        args.force,
                        refinements=int(refinements),
                        feedback=float(arm["feedback"]),
                        reference_scale=float(arm["reference_scale"]),
                    )
                    wall = run_epsr(binary, run_dir)
                    case_summary[name] = {
                        "status": "completed",
                        "wall_seconds": wall,
                        "seed": seed,
                        "refinements": int(refinements),
                        "attempted_moves": int(refinements)
                        * moves_per_refinement,
                    }
                    summary_path.parent.mkdir(parents=True, exist_ok=True)
                    summary_path.write_text(
                        json.dumps(summary, indent=2, sort_keys=True) + "\n"
                    )
                    print(
                        json.dumps(
                            {case.name: {name: case_summary[name]}}, indent=2
                        )
                    )
    raise SystemExit(0)

if args.convergence_pilot or args.convergence_ensemble:
    if args.convergence_pilot and args.convergence_ensemble:
        raise SystemExit("choose only one convergence protocol")
    if args.convergence_ensemble and args.native_rminex_control:
        raise SystemExit("the rminex sensitivity control is not part of the ensemble")
    protocol_name = (
        "epsr-convergence-ensemble.toml"
        if args.convergence_ensemble
        else "epsr-convergence-pilot.toml"
    )
    protocol = tomllib.loads(
        (case_root / "expected" / protocol_name).read_text()
    )
    convergence_root = case_root / "results" / (
        "epsr-convergence-ensemble"
        if args.convergence_ensemble
        else "epsr-convergence-pilot"
    )
    method = (
        "native-epsr26-rminex-control"
        if args.native_rminex_control
        else "native-epsr26"
    )
    configured_seeds = (
        [int(seed) for seed in protocol["design"]["seeds"]]
        if args.convergence_ensemble
        else [int(protocol["design"]["seed"])]
    )
    selected_seeds = args.convergence_seed or configured_seeds
    unknown_seeds = set(selected_seeds) - set(configured_seeds)
    if unknown_seeds:
        raise SystemExit(f"seeds not present in {protocol_name}: {sorted(unknown_seeds)}")
    moves_per_refinement = int(
        protocol["design"]["moves_per_refinement"]
        if args.convergence_ensemble
        else protocol["sampling"]["moves_per_refinement"]
    )
    checkpoints = args.checkpoint or (
        protocol["methods"]["native_epsr26"]["checkpoints"]
        if args.convergence_ensemble
        else protocol["design"]["checkpoints"]
    )
    for seed in selected_seeds:
        summary_path = convergence_root / (
            f"{method}-seed-{seed}-run-summary.json"
            if args.convergence_ensemble
            else (
                "native-epsr26-rminex-control-run-summary.json"
                if args.native_rminex_control
                else "native-epsr-run-summary.json"
            )
        )
        summary = {
            "program": "EPSR26",
            "binary": str(binary),
            "threads": 1,
            "sampling": "independent deterministic prefixes",
            "seed": seed,
            "cases": {},
        }
        for case in sorted(fixture_root.glob("target-*_*")):
            targets = []
            for filename in (
                "epsr-native-target-neutron.dat",
                "epsr-native-target-xray.dat",
            ):
                rows = read_iq(case / filename)
                targets.append([row for row in rows if 0.5 <= row[0] < 25.0])
            prefix_root = convergence_root / case.name / method / f"seed-{seed}"
            prefix_root.mkdir(parents=True, exist_ok=True)
            case_summary = summary["cases"].setdefault(case.name, {})
            for refinements in checkpoints:
                name = f"iter-{int(refinements):03d}"
                run_dir = prefix_root / name
                if args.only_missing and (run_dir / "DTBsilica.EPSR.v01").is_file():
                    continue
                run_dir = prepare_run(
                    source,
                    prefix_root,
                    name,
                    case / "cross-start.data",
                    moves_per_refinement,
                    seed,
                    tuple(targets),
                    args.force,
                    refinements=int(refinements),
                    rminex=(
                        {
                            ("Si", "Si"): 2.0,
                            ("Si", "O"): 1.35,
                            ("O", "O"): 2.0,
                        }
                        if args.native_rminex_control
                        else None
                    ),
                )
                wall = run_epsr(binary, run_dir)
                case_summary[name] = {
                    "status": "completed",
                    "wall_seconds": wall,
                    "seed": seed,
                    "refinements": int(refinements),
                    "attempted_moves": int(refinements) * moves_per_refinement,
                }
                summary_path.parent.mkdir(parents=True, exist_ok=True)
                summary_path.write_text(
                    json.dumps(summary, indent=2, sort_keys=True) + "\n"
                )
                print(json.dumps({case.name: {name: case_summary[name]}}, indent=2))
    raise SystemExit(0)

if args.ensemble:
    protocol = tomllib.loads(
        (case_root / "expected/multiseed-comparison.toml").read_text()
    )
    ensemble_root = case_root / "results/multiseed-comparison"
    summary_path = ensemble_root / "native-epsr-run-summary.json"
    summary = (
        json.loads(summary_path.read_text())
        if summary_path.is_file()
        else {"program": "EPSR26", "binary": str(binary), "threads": 1, "cases": {}}
    )
    for case in sorted(fixture_root.glob("target-*_*")):
        targets = []
        for filename in ("epsr-native-target-neutron.dat", "epsr-native-target-xray.dat"):
            rows = read_iq(case / filename)
            targets.append([row for row in rows if 0.5 <= row[0] < 25.0])
        method_root = ensemble_root / case.name / "native-epsr26"
        method_root.mkdir(parents=True, exist_ok=True)
        case_summary = summary["cases"].setdefault(case.name, {})
        for seed in protocol["budget"]["seeds"]:
            name = f"seed-{seed}"
            run_dir = method_root / name
            if args.only_missing and (run_dir / "DTBsilica.EPSR.v01").is_file():
                continue
            run_dir = prepare_run(
                source,
                method_root,
                name,
                case / "cross-start.data",
                int(protocol["budget"]["moves"]),
                int(seed),
                tuple(targets),
                args.force,
            )
            wall = run_epsr(binary, run_dir)
            case_summary[name] = {
                "status": "completed",
                "wall_seconds": wall,
                "seed": int(seed),
            }
            summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
            print(json.dumps({case.name: {name: case_summary[name]}}, indent=2))
    raise SystemExit(0)

summary_path = fixture_root / "native-epsr-run-summary.json"
if args.zero_only and summary_path.is_file():
    summary = json.loads(summary_path.read_text())
    summary.update(
        {
            "program": "EPSR26",
            "binary": str(binary),
            "moves": args.moves,
            "seed": args.seed,
        }
    )
else:
    summary = {
        "program": "EPSR26",
        "binary": str(binary),
        "moves": args.moves,
        "seed": args.seed,
        "cases": {},
    }
for case in sorted(fixture_root.glob("target-*_*")):
    case_summary = summary["cases"].get(case.name, {})
    provisional = provisional_targets(case, source)
    zero = prepare_run(
        source,
        case,
        "native-epsr-hidden-zero-move",
        case / "hidden-target.data",
        0,
        args.seed,
        provisional,
        args.force,
    )
    case_summary["native_target_generation_wall_seconds"] = run_epsr(binary, zero)
    native_targets = calculated_targets(zero / "DTBsilica.EPSR.v01", provisional)
    case_summary["native_vs_provisional_neutron_rms"] = rms(
        native_targets[0], provisional[0]
    )
    case_summary["native_vs_provisional_xray_rms"] = rms(
        native_targets[1], provisional[1]
    )
    write_epsr_data(
        case / "epsr-native-target-neutron.dat",
        "EPSR-native hidden-coordinate neutron i(Q)",
        native_targets[0],
    )
    write_epsr_data(
        case / "epsr-native-target-xray.dat",
        "EPSR-native hidden-coordinate X-ray i(Q)",
        native_targets[1],
    )
    zero = prepare_run(
        source,
        case,
        "native-epsr-hidden-zero-move",
        case / "hidden-target.data",
        0,
        args.seed,
        native_targets,
        True,
    )
    case_summary["hidden_zero_move_wall_seconds"] = run_epsr(binary, zero)
    checked = calculated_targets(zero / "DTBsilica.EPSR.v01", native_targets)
    case_summary["hidden_zero_move_neutron_rms"] = rms(checked[0], native_targets[0])
    case_summary["hidden_zero_move_xray_rms"] = rms(checked[1], native_targets[1])
    cross_zero = prepare_run(
        source,
        case,
        "native-epsr-cross-zero-move",
        case / "cross-start.data",
        0,
        args.seed,
        native_targets,
        args.force,
    )
    case_summary["cross_zero_move_wall_seconds"] = run_epsr(binary, cross_zero)
    if not args.zero_only:
        smoke = prepare_run(
            source,
            case,
            "native-epsr-joint",
            case / "cross-start.data",
            args.moves,
            args.seed,
            native_targets,
            args.force,
        )
        case_summary["joint_wall_seconds"] = run_epsr(binary, smoke)
    summary["cases"][case.name] = case_summary

summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(json.dumps(summary, indent=2, sort_keys=True))
