#!/usr/bin/env python3
"""Run matched-start native EPSR26 and rsmith liquid-Ga replicate ensembles."""

from __future__ import annotations

import argparse
import bisect
import concurrent.futures
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import time
from pathlib import Path

import numpy as np


DEFAULT_SEEDS = tuple(range(20260802, 20260812))
NEUTRON_TOTAL_WEIGHT = 0.531149447
KJ_PER_MOL_PER_EV = 96.48533212


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


def seed_epsr_ato(path: Path, seed: int) -> None:
    lines = path.read_text().splitlines()
    candidates = []
    for index, line in enumerate(lines):
        fields = line.split()
        if len(fields) == 35:
            try:
                [int(value) for value in fields]
            except ValueError:
                continue
            candidates.append(index)
    if len(candidates) != 1:
        raise ValueError(f"expected one EPSR RNG-state record, found {len(candidates)}")
    lines[candidates[0]] = " " + " ".join(str(value) for value in [0, -seed, 0, *([0] * 32)])
    path.write_text("\n".join(lines) + "\n")


def read_curve(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        try:
            row = tuple(float(value) for value in line.split())
        except ValueError:
            continue
        if len(row) >= 2:
            rows.append(row)
    return rows


def interpolate(curve, x, column=1):
    xs = [row[0] for row in curve]
    hi = bisect.bisect_left(xs, x)
    if hi < len(curve) and abs(xs[hi] - x) < 1.0e-8:
        return curve[hi][column]
    if hi == 0 or hi == len(curve):
        raise ValueError(f"{x} outside candidate grid")
    lo = hi - 1
    fraction = (x - xs[lo]) / (xs[hi] - xs[lo])
    return curve[lo][column] + fraction * (curve[hi][column] - curve[lo][column])


def curve_metrics(reference, candidate, fit_min, fit_max, scale=1.0):
    selected = [row for row in reference if fit_min <= row[0] <= fit_max]
    expected = [row[1] for row in selected]
    residual = [scale * interpolate(candidate, row[0]) - row[1] for row in selected]
    dynamic_range = max(expected) - min(expected)
    rms = math.sqrt(sum(value * value for value in residual) / len(residual))
    maximum = max(abs(value) for value in residual)
    return {
        "n_points": len(selected),
        "rms": rms,
        "max_absolute": maximum,
        "rms_over_dynamic_range": rms / dynamic_range,
        "max_over_dynamic_range": maximum / dynamic_range,
    }


def block_average(curve, width):
    start = math.ceil(curve[0][0] / width) * width
    stop = math.floor(curve[-1][0] / width) * width
    result = []
    index = 0
    lower = start
    while lower <= stop + 1.0e-10:
        upper = lower + width
        values = []
        while index < len(curve) and curve[index][0] < upper - 1.0e-10:
            if curve[index][0] >= lower - 1.0e-10:
                values.append(curve[index][1])
            index += 1
        if values:
            result.append((lower + 0.5 * width, sum(values) / len(values)))
        lower = upper
    return result


def convert_ato(path: Path, output: Path):
    lines = path.read_text().splitlines()
    atom_count, box = int(lines[0].split()[0]), float(lines[0].split()[1])
    records = lines[2:]
    atoms = []
    for index in range(atom_count):
        fields = records[5 * index].split()
        species = records[5 * index + 1].split()[0]
        xyz = tuple((float(fields[axis]) + 0.5 * box) % box for axis in (1, 2, 3))
        atoms.append((species, xyz))
    with output.open("w") as stream:
        stream.write("EPSR26 LiquidGa50C matched stochastic start\n\n")
        stream.write(f"{atom_count} atoms\n1 atom types\n\n")
        for axis in "xyz":
            stream.write(f"0.0 {box:.12f} {axis}lo {axis}hi\n")
        stream.write("\nAtoms # charge\n\n")
        for atom_id, (species, xyz) in enumerate(atoms, 1):
            if species != "Ga":
                raise ValueError(f"unexpected species {species}")
            stream.write(f"{atom_id} 1 0.0 {xyz[0]:.12f} {xyz[1]:.12f} {xyz[2]:.12f}\n")
    return atom_count, box


def rsmith_epoch_ensemble(run_dir: Path, box: float, iterations: int):
    """Average all saved rsmith epoch configurations on the production grid."""
    gr_cache = run_dir / "epoch_mean_gr.dat"
    iq_cache = run_dir / "epoch_mean_iq.dat"
    if gr_cache.is_file() and iq_cache.is_file():
        return read_curve(gr_cache), read_curve(iq_cache)

    cutoff = 22.86
    nbins = 762
    dr = cutoff / nbins
    radius = (np.arange(nbins) + 0.5) * dr
    number_density = 5000 / box**3
    normalisation = 5000 * number_density * 4.0 * math.pi * radius**2 * dr
    gr_sum = np.zeros(nbins)
    for epoch in range(1, iterations + 1):
        path = run_dir / f"refined_iter_{epoch}.xyz"
        lines = path.read_text().splitlines()
        positions = np.asarray(
            [[float(value) for value in line.split()[1:4]] for line in lines[2:5002]]
        )
        histogram = np.zeros(nbins, dtype=np.int64)
        for index in range(len(positions) - 1):
            displacement = positions[index + 1 :] - positions[index]
            displacement -= box * np.rint(displacement / box)
            distance = np.sqrt(np.einsum("ij,ij->i", displacement, displacement))
            bins = np.floor(distance[distance < cutoff] / dr).astype(np.int64)
            histogram += np.bincount(bins, minlength=nbins)
        gr_sum += 2.0 * histogram / normalisation
    mean_gr = gr_sum / iterations

    q = np.arange(600) * 0.05
    iq = np.zeros_like(q)
    integrand = radius * (mean_gr - 1.0)
    for index in range(1, len(q)):
        iq[index] = (
            4.0
            * math.pi
            * number_density
            * dr
            * np.sum(integrand * np.sin(q[index] * radius))
            / q[index]
        )
    gr_cache.write_text(
        "# r mean_g_GaGa\n"
        + "".join(f"{r:.12g} {value:.12g}\n" for r, value in zip(radius, mean_gr))
    )
    iq_cache.write_text(
        "# Q mean_S(Q)-1\n"
        + "".join(f"{q_value:.12g} {value:.12g}\n" for q_value, value in zip(q, iq))
    )
    return list(zip(radius.tolist(), mean_gr.tolist())), list(zip(q.tolist(), iq.tolist()))


def native_acceptance(report: str) -> float:
    block = report.split("Number of rejects in MOLFIT:-", 1)[1].splitlines()
    attempted = [int(value) for value in block[1].split()][-1]
    rejected = [int(value) for value in block[2].split()][-1]
    return (attempted - rejected) / attempted


def rsmith_acceptance(log: str, iterations: int) -> float:
    match = re.findall(
        rf"EPSR iter {iterations}:.*?acceptance = ([0-9.]+)%", log
    )
    if len(match) != 1:
        raise ValueError("could not identify final rsmith acceptance")
    return float(match[0]) / 100.0


def run_native(source: Path, binary: Path, output: Path, seed: int, iterations: int, force: bool):
    if output.exists():
        if not force:
            return
        shutil.rmtree(output)
    shutil.copytree(source, output)
    input_path = output / "LiquidGa50C.EPSR.inp"
    text = input_path.read_text()
    for name, value in (
        ("potfac", "1.0 1.0"),
        ("num_threds", "1"),
        ("ireset", "2"),
        ("iinit", "1"),
        ("ntimes", "5"),
        ("niter", str(iterations)),
        ("nsumt", "0"),
        ("ngrsamples", "5000"),
    ):
        text = replace_setting(text, name, value)
    input_path.write_text(text)
    seed_epsr_ato(output / "LiquidGa50C.ato", seed)
    for path in output.glob("LiquidGa50C.EPSR.*"):
        if path.suffix not in {".inp", ".inpa"}:
            path.unlink()
    started = time.perf_counter()
    completed = subprocess.run(
        [str(binary), str(output) + "/", "epsr", "LiquidGa50C"],
        cwd=output,
        text=True,
        input="",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    elapsed = time.perf_counter() - started
    (output / "native.log").write_text(completed.stdout)
    (output / "wall_seconds.txt").write_text(f"{elapsed:.9f}\n")
    if completed.returncode != 0:
        raise RuntimeError(f"native EPSR seed {seed} failed; see {output / 'native.log'}")
    report = (output / "LiquidGa50C.EPSR.out").read_text()
    match = re.search(r"No\. of configurations in sum\s+(\d+)", report)
    if match is None or int(match.group(1)) != iterations:
        raise RuntimeError(f"native EPSR seed {seed} accumulated the wrong number of configurations")


def run_rsmith(
    source_ato: Path,
    target_iq: Path,
    binary: Path,
    output: Path,
    seed: int,
    iterations: int,
    moves_per_iteration: int,
    force: bool,
):
    if output.exists():
        if not force:
            return
        shutil.rmtree(output)
    output.mkdir(parents=True)
    atom_count, box = convert_ato(source_ato, output / "liquid-ga.data")
    if atom_count != 5000:
        raise ValueError("unexpected LiquidGa50C atom count")
    target = read_curve(target_iq)
    target_path = output / "target-partial-iq.dat"
    target_path.write_text(
        "# Q S(Q)-1\n" + "".join(f"{row[0]:.12g} {row[1]:.12g}\n" for row in target)
    )
    epsilon_ev = 0.8 / KJ_PER_MOL_PER_EV
    config = output / "run.toml"
    config.write_text(
        f'''[system]
structure = "{output / 'liquid-ga.data'}"
format = "lammps"
types = {{ "1" = "Ga" }}

[data]
[data.neutron_sq]
file = "{target_path}"
sigma = 1.0
convention = "iq"
fit_min = 1.3
fit_max = 29.95

[rmc]
max_moves = {moves_per_iteration}
max_step = 0.3
checkpoint_every = 1000000000
seed = {seed}
print_every = {moves_per_iteration}
target_acceptance = 0.25
adjust_step_every = 5000

[sq]
qmin = 0.0
qmax = 30.0
nq = 600
lorch = false
rdf_cutoff = 22.86
rdf_nbins = 762

[potential]
cutoff = 12.0
weight = 1.0

[[potential.lennard_jones]]
pair = "Ga-Ga"
epsilon = {epsilon_ev:.15g}
sigma = 3.0

[epsr]
mode = "pure"
iterations = {iterations}
feedback = 0.9
smooth_sigma = 0.03
moves_per_iteration = {moves_per_iteration}
temperature = 0.0278339856
min_r = 1.0
'''
    )
    started = time.perf_counter()
    completed = subprocess.run(
        [str(binary), str(config), "--output-dir", str(output), "--quiet"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env={**os.environ, "RAYON_NUM_THREADS": "1"},
        check=False,
    )
    elapsed = time.perf_counter() - started
    (output / "driver.log").write_text(completed.stdout)
    (output / "wall_seconds.txt").write_text(f"{elapsed:.9f}\n")
    if completed.returncode != 0:
        raise RuntimeError(f"rsmith seed {seed} failed; see {output / 'rsmith.log'}")


def distribution(values):
    ordered = sorted(values)
    if len(ordered) == 1:
        q1 = q3 = ordered[0]
    else:
        q1 = statistics.median(ordered[: len(ordered) // 2])
        q3 = statistics.median(ordered[(len(ordered) + 1) // 2 :])
    return {
        "min": min(values),
        "median": statistics.median(values),
        "max": max(values),
        "mean": statistics.fmean(values),
        "sample_stddev": statistics.stdev(values) if len(values) > 1 else 0.0,
        "q1": q1,
        "q3": q3,
    }


def summarize(case_root: Path, seeds, iterations, moves_per_iteration):
    local = case_root / "reference/local/upstream"
    target_total = read_curve(local / "LiquidGa50C.EPSR.t01")
    target_gr = block_average(read_curve(local / "LiquidGa50C.EPSR.g01"), 0.12)
    runs = []
    for seed in seeds:
        native = case_root / f"results/stochastic/native/seed-{seed}"
        rsmith = case_root / f"results/stochastic/rsmith/seed-{seed}"
        native_iq = read_curve(native / "LiquidGa50C.EPSR.u01")
        _, box = convert_ato(local / "LiquidGa50C.ato", rsmith / "summary-start.data")
        rsmith_gr_raw, rsmith_iq = rsmith_epoch_ensemble(rsmith, box, iterations)
        native_gr = block_average(read_curve(native / "LiquidGa50C.EPSR.g01"), 0.12)
        rsmith_gr = block_average(rsmith_gr_raw, 0.12)
        native_report = (native / "LiquidGa50C.EPSR.out").read_text()
        rsmith_log = (rsmith / "rsmith.log").read_text()
        runs.append(
            {
                "seed": seed,
                "native": {
                    "total_fit": curve_metrics(target_total, native_iq, 1.3, 29.95),
                    "gr_vs_worked_ensemble": curve_metrics(target_gr, native_gr, 1.5, 10.0),
                    "final_cycle_acceptance": native_acceptance(native_report),
                    "wall_seconds": float((native / "wall_seconds.txt").read_text()),
                },
                "rsmith": {
                    "total_fit": curve_metrics(
                        target_total, rsmith_iq, 1.3, 29.95, scale=NEUTRON_TOTAL_WEIGHT
                    ),
                    "gr_vs_worked_ensemble": curve_metrics(target_gr, rsmith_gr, 1.5, 10.0),
                    "final_epoch_acceptance": rsmith_acceptance(rsmith_log, iterations),
                    "wall_seconds": float((rsmith / "wall_seconds.txt").read_text()),
                },
                "matched": {
                    "iq": curve_metrics(
                        read_curve(native / "LiquidGa50C.EPSR.f01"), rsmith_iq, 1.3, 29.95
                    ),
                    "gr": curve_metrics(native_gr, rsmith_gr, 1.5, 10.0),
                },
            }
        )
    summary = {
        "protocol": {
            "seeds": list(seeds),
            "replicates": len(seeds),
            "iterations": iterations,
            "moves_per_iteration": moves_per_iteration,
            "moves_per_replicate": iterations * moves_per_iteration,
            "starting_coordinates": "identical EPSR26 worked-example endpoint for every method and seed",
        },
        "runs": runs,
        "distributions": {},
    }
    for method in ("native", "rsmith"):
        summary["distributions"][method] = {
            "total_fit_rms_over_dynamic_range": distribution(
                [run[method]["total_fit"]["rms_over_dynamic_range"] for run in runs]
            ),
            "gr_vs_worked_rms_over_dynamic_range": distribution(
                [run[method]["gr_vs_worked_ensemble"]["rms_over_dynamic_range"] for run in runs]
            ),
            "wall_seconds": distribution([run[method]["wall_seconds"] for run in runs]),
            "final_acceptance": distribution(
                [
                    run[method][
                        "final_cycle_acceptance" if method == "native" else "final_epoch_acceptance"
                    ]
                    for run in runs
                ]
            ),
        }
    summary["distributions"]["matched"] = {
        "iq_rms_over_dynamic_range": distribution(
            [run["matched"]["iq"]["rms_over_dynamic_range"] for run in runs]
        ),
        "gr_rms_over_dynamic_range": distribution(
            [run["matched"]["gr"]["rms_over_dynamic_range"] for run in runs]
        ),
    }
    output = case_root / "results/stochastic/summary.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary["distributions"], indent=2, sort_keys=True))
    print(f"Full stochastic summary: {output}")


parser = argparse.ArgumentParser()
parser.add_argument("--seeds", default=",".join(str(seed) for seed in DEFAULT_SEEDS))
parser.add_argument("--iterations", type=int, default=40)
parser.add_argument("--moves-per-iteration", type=int, default=25000)
parser.add_argument("--force", action="store_true")
parser.add_argument("--jobs", type=int, default=1)
parser.add_argument("--native-only", action="store_true")
parser.add_argument("--rsmith-only", action="store_true")
parser.add_argument("--summarize-only", action="store_true")
parser.add_argument("--binary", type=Path)
args = parser.parse_args()
if args.native_only and args.rsmith_only:
    raise SystemExit("--native-only and --rsmith-only are mutually exclusive")
seeds = tuple(int(value) for value in args.seeds.split(",") if value)
if not seeds or args.iterations <= 0 or args.moves_per_iteration <= 0 or args.jobs <= 0:
    raise SystemExit("seeds, iterations, and moves per iteration must be positive")

case_root = Path(__file__).resolve().parents[1]
repo_root = case_root.parents[1]
local = case_root / "reference/local"
source = local / "upstream"
record_file = local / "IMPORT.txt"
if not source.is_dir() or not record_file.is_file():
    raise SystemExit("run import_local_reference.sh --accept-local-testing-terms first")
record = dict(line.split("=", 1) for line in record_file.read_text().splitlines() if "=" in line)
epsr_binary = Path(record["binary"])
rsmith_binary = args.binary or repo_root / "target/release/rsmith"
if not epsr_binary.is_file() or not rsmith_binary.is_file():
    raise SystemExit("native EPSR or release rsmith executable is missing")

def run_seed(seed):
    print(f"seed {seed}: native EPSR", flush=True)
    if not args.rsmith_only:
        run_native(
            source,
            epsr_binary,
            case_root / f"results/stochastic/native/seed-{seed}",
            seed,
            args.iterations,
            args.force,
        )
    print(f"seed {seed}: rsmith", flush=True)
    if not args.native_only:
        run_rsmith(
            source / "LiquidGa50C.ato",
            source / "LiquidGa50C.EPSR.q01",
            rsmith_binary,
            case_root / f"results/stochastic/rsmith/seed-{seed}",
            seed,
            args.iterations,
            args.moves_per_iteration,
            args.force,
        )


if not args.summarize_only:
    if args.jobs == 1:
        for seed in seeds:
            run_seed(seed)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = [executor.submit(run_seed, seed) for seed in seeds]
            for future in concurrent.futures.as_completed(futures):
                future.result()

if not args.native_only and not args.rsmith_only:
    summarize(case_root, seeds, args.iterations, args.moves_per_iteration)
