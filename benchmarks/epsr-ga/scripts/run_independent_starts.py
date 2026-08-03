#!/usr/bin/env python3
"""Run paired EPSR26/rsmith convergence tests from independent liquid starts.

The original fixed-budget benchmark deliberately starts from the supplied EPSR
endpoint.  That configuration already fits the measured total well, so it
cannot answer a time-to-target question.  This script generates neutral,
periodic hard-sphere starts and evaluates both programs at frozen iteration
checkpoints from the identical coordinates for each seed.

Raw upstream inputs, generated structures, and executable outputs stay below
``results/`` and are ignored.  Only aggregate, non-redistributive observations
belong in ``expected/``.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np

from run_stochastic import (
    KJ_PER_MOL_PER_EV,
    NEUTRON_TOTAL_WEIGHT,
    block_average,
    curve_metrics,
    distribution,
    read_curve,
    rsmith_epoch_ensemble,
    run_native,
    run_rsmith,
)


DEFAULT_SEEDS = tuple(range(20260812, 20260822))
DEFAULT_CHECKPOINTS = (0, 1, 2, 4, 8, 16, 24, 32, 40)
TARGET_RESIDUALS = (0.05, 0.04, 0.03, 0.025, 0.02, 0.015)
ATOM_COUNT = 5000
BOX_ANGSTROM = 45.725258
HARD_SPHERE_MIN_ANGSTROM = 2.25
FIRST_SHELL_CUTOFF_ANGSTROM = 3.93
RDF_CUTOFF_ANGSTROM = 22.86
RDF_BINS = 762
MOVES_PER_ITERATION = 25000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_ato_positions(path: Path):
    lines = path.read_text().splitlines()
    atom_count, box = int(lines[0].split()[0]), float(lines[0].split()[1])
    records = lines[2:]
    positions = np.empty((atom_count, 3), dtype=np.float64)
    for index in range(atom_count):
        fields = records[5 * index].split()
        positions[index] = [
            (float(fields[axis]) + 0.5 * box) % box for axis in (1, 2, 3)
        ]
    return positions, box


def read_xyz_positions(path: Path):
    lines = path.read_text().splitlines()
    atom_count = int(lines[0])
    positions = np.asarray(
        [
            [float(value) for value in line.split()[1:4]]
            for line in lines[2 : atom_count + 2]
        ],
        dtype=np.float64,
    )
    return positions


def write_ato_with_positions(template: Path, output: Path, positions, seed: int):
    lines = template.read_text().splitlines()
    atom_count, box = int(lines[0].split()[0]), float(lines[0].split()[1])
    if len(positions) != atom_count:
        raise ValueError("position count does not match EPSR ATO template")
    for index, position in enumerate(positions):
        line_index = 2 + 5 * index
        fields = lines[line_index].split()
        centred = ((np.asarray(position) + 0.5 * box) % box) - 0.5 * box
        fields[1:4] = [f"{value:.10E}" for value in centred]
        lines[line_index] = " " + " ".join(fields)

    rng_records = []
    for index, line in enumerate(lines):
        fields = line.split()
        if len(fields) != 35:
            continue
        try:
            [int(value) for value in fields]
        except ValueError:
            continue
        rng_records.append(index)
    if len(rng_records) != 1:
        raise ValueError(f"expected one EPSR RNG record, found {len(rng_records)}")
    lines[rng_records[0]] = " " + " ".join(
        str(value) for value in [0, -seed, 0, *([0] * 32)]
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n")


def write_xyz(path: Path, positions, box: float):
    path.write_text(
        f"{len(positions)}\n"
        f"independent periodic hard-sphere start; box={box:.9f} A\n"
        + "".join(f"Ga {x:.12f} {y:.12f} {z:.12f}\n" for x, y, z in positions)
    )


def generate_hard_sphere_positions(
    seed: int, atom_count: int, box: float, min_distance: float
):
    """Generate a periodic random-sequential hard-sphere configuration."""
    rng = np.random.default_rng(seed)
    cells_per_axis = max(1, int(math.floor(box / min_distance)))
    cell_width = box / cells_per_axis
    cells: dict[tuple[int, int, int], list[int]] = {}
    positions = np.empty((atom_count, 3), dtype=np.float64)
    min_squared = min_distance * min_distance
    attempts = 0
    placed = 0
    max_attempts = atom_count * 20000

    while placed < atom_count:
        if attempts >= max_attempts:
            raise RuntimeError(
                f"hard-sphere insertion stalled after {attempts} attempts at {placed}/{atom_count} atoms"
            )
        attempts += 1
        candidate = rng.random(3) * box
        cell = tuple((candidate / cell_width).astype(int) % cells_per_axis)
        valid = True
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neighbour = (
                        (cell[0] + dx) % cells_per_axis,
                        (cell[1] + dy) % cells_per_axis,
                        (cell[2] + dz) % cells_per_axis,
                    )
                    for other in cells.get(neighbour, ()):
                        displacement = candidate - positions[other]
                        displacement -= box * np.rint(displacement / box)
                        if float(np.dot(displacement, displacement)) < min_squared:
                            valid = False
                            break
                    if not valid:
                        break
                if not valid:
                    break
            if not valid:
                break
        if valid:
            positions[placed] = candidate
            cells.setdefault(cell, []).append(placed)
            placed += 1
    return positions, attempts


def pair_histogram(positions, box: float, cutoff: float, bins: int):
    dr = cutoff / bins
    histogram = np.zeros(bins, dtype=np.int64)
    for index in range(len(positions) - 1):
        displacement = positions[index + 1 :] - positions[index]
        displacement -= box * np.rint(displacement / box)
        distances = np.sqrt(np.einsum("ij,ij->i", displacement, displacement))
        selected = distances[distances < cutoff]
        histogram += np.bincount(
            np.floor(selected / dr).astype(np.int64), minlength=bins
        )
    return histogram


def gr_and_iq(positions, box: float):
    dr = RDF_CUTOFF_ANGSTROM / RDF_BINS
    radius = (np.arange(RDF_BINS) + 0.5) * dr
    density = len(positions) / box**3
    histogram = pair_histogram(positions, box, RDF_CUTOFF_ANGSTROM, RDF_BINS)
    normalisation = len(positions) * density * 4.0 * math.pi * radius**2 * dr
    gr = 2.0 * histogram / normalisation
    q = np.arange(600) * 0.05
    iq = np.zeros_like(q)
    integrand = radius * (gr - 1.0)
    for index in range(1, len(q)):
        iq[index] = (
            4.0
            * math.pi
            * density
            * dr
            * np.sum(integrand * np.sin(q[index] * radius))
            / q[index]
        )
    return list(zip(radius.tolist(), gr.tolist())), list(zip(q.tolist(), iq.tolist()))


def quadratic_peak(curve, lower: float, upper: float):
    points = np.asarray([(row[0], row[1]) for row in curve if lower <= row[0] <= upper])
    peak = int(np.argmax(points[:, 1]))
    if peak == 0 or peak == len(points) - 1:
        return {
            "position_angstrom": float(points[peak, 0]),
            "height": float(points[peak, 1]),
        }
    coefficients = np.polyfit(
        points[peak - 1 : peak + 2, 0], points[peak - 1 : peak + 2, 1], 2
    )
    position = -coefficients[1] / (2.0 * coefficients[0])
    height = np.polyval(coefficients, position)
    return {"position_angstrom": float(position), "height": float(height)}


def rdf_metrics(curve, density: float):
    selected = [
        (row[0], row[1]) for row in curve if row[0] <= FIRST_SHELL_CUTOFF_ANGSTROM
    ]
    radius = np.asarray([row[0] for row in selected])
    gr = np.asarray([row[1] for row in selected])
    coordination = (
        4.0 * math.pi * density * float(np.trapezoid(radius * radius * gr, radius))
    )
    return {
        "first_peak": quadratic_peak(curve, 2.2, 3.3),
        "first_shell_cutoff_angstrom": FIRST_SHELL_CUTOFF_ANGSTROM,
        "coordination_integral": coordination,
    }


def configuration_metrics(positions, box: float):
    nearest = np.empty(len(positions), dtype=np.float64)
    coordination = np.empty(len(positions), dtype=np.int64)
    block = 192
    for start in range(0, len(positions), block):
        stop = min(start + block, len(positions))
        displacement = positions[start:stop, None, :] - positions[None, :, :]
        displacement -= box * np.rint(displacement / box)
        distances_squared = np.einsum("ijk,ijk->ij", displacement, displacement)
        distances_squared[np.arange(stop - start), np.arange(start, stop)] = np.inf
        nearest[start:stop] = np.sqrt(np.min(distances_squared, axis=1))
        coordination[start:stop] = np.sum(
            distances_squared < FIRST_SHELL_CUTOFF_ANGSTROM**2, axis=1
        )
    return {
        "nearest_neighbor_angstrom": {
            "mean": float(np.mean(nearest)),
            "sample_stddev": float(np.std(nearest, ddof=1)),
            "q05": float(np.quantile(nearest, 0.05)),
            "median": float(np.median(nearest)),
            "q95": float(np.quantile(nearest, 0.95)),
            "min": float(np.min(nearest)),
        },
        "coordination_within_3_93_angstrom": {
            "mean": float(np.mean(coordination)),
            "sample_stddev": float(np.std(coordination, ddof=1)),
            "min": int(np.min(coordination)),
            "max": int(np.max(coordination)),
        },
        "fraction_nearest_neighbor_below_2_3_angstrom": float(np.mean(nearest < 2.3)),
    }


def native_energy_per_atom(run_dir: Path):
    rows = [
        line.split()
        for line in (run_dir / "LiquidGa50C.EPSR.erg").read_text().splitlines()
    ]
    rows = [row for row in rows if row]
    return float(rows[-1][0])


def rsmith_energy_per_atom(run_dir: Path):
    matches = re.findall(
        r"Energy MC complete:.*?E = ([+-]?[0-9.eE]+) eV",
        (run_dir / "rsmith.log").read_text(),
    )
    if not matches:
        raise ValueError(f"could not find final energy in {run_dir / 'rsmith.log'}")
    return float(matches[-1]) * KJ_PER_MOL_PER_EV / ATOM_COUNT


def score_checkpoint(case_root: Path, seed: int, epoch: int):
    source = case_root / "reference/local/upstream"
    root = case_root / "results/independent-starts"
    start_ato = root / f"starts/seed-{seed}.ato"
    target_total = read_curve(source / "LiquidGa50C.EPSR.t01")
    target_gr = block_average(read_curve(source / "LiquidGa50C.EPSR.g01"), 0.12)
    _, box = read_ato_positions(start_ato)

    if epoch == 0:
        positions, _ = read_ato_positions(start_ato)
        gr_raw, iq = gr_and_iq(positions, box)
        gr = block_average(gr_raw, 0.12)
        total_fit = curve_metrics(
            target_total, iq, 1.3, 29.95, scale=NEUTRON_TOTAL_WEIGHT
        )
        structural = rdf_metrics(gr_raw, len(positions) / box**3)
        configuration = configuration_metrics(positions, box)
        shared = {
            "wall_seconds": 0.0,
            "total_fit": total_fit,
            "gr_vs_worked_ensemble": curve_metrics(target_gr, gr, 1.5, 10.0),
            "rdf": structural,
            "configuration": configuration,
        }
        return {
            "seed": seed,
            "epoch": 0,
            "moves": 0,
            "native": shared,
            "rsmith": shared,
        }

    native = root / f"checkpoints/native/seed-{seed}/epoch-{epoch}"
    rsmith = root / f"checkpoints/rsmith/seed-{seed}/epoch-{epoch}"
    native_gr_raw = read_curve(native / "LiquidGa50C.EPSR.g01")
    native_gr = block_average(native_gr_raw, 0.12)
    native_total = read_curve(native / "LiquidGa50C.EPSR.u01")
    rsmith_gr_raw, rsmith_iq = rsmith_epoch_ensemble(rsmith, box, epoch)
    rsmith_gr = block_average(rsmith_gr_raw, 0.12)

    result = {"seed": seed, "epoch": epoch, "moves": epoch * MOVES_PER_ITERATION}
    for method, total, gr_raw, gr, scale, run_dir in (
        ("native", native_total, native_gr_raw, native_gr, 1.0, native),
        ("rsmith", rsmith_iq, rsmith_gr_raw, rsmith_gr, NEUTRON_TOTAL_WEIGHT, rsmith),
    ):
        result[method] = {
            "wall_seconds": float((run_dir / "wall_seconds.txt").read_text()),
            "total_fit": curve_metrics(target_total, total, 1.3, 29.95, scale=scale),
            "gr_vs_worked_ensemble": curve_metrics(target_gr, gr, 1.5, 10.0),
            "rdf": rdf_metrics(gr_raw, ATOM_COUNT / box**3),
        }
    if epoch == max(DEFAULT_CHECKPOINTS):
        native_positions, _ = read_ato_positions(native / "LiquidGa50C.ato")
        rsmith_positions = read_xyz_positions(rsmith / "refined.xyz")
        result["native"]["configuration"] = configuration_metrics(native_positions, box)
        result["rsmith"]["configuration"] = configuration_metrics(rsmith_positions, box)
        result["native"]["combined_energy_kj_mol_per_atom"] = native_energy_per_atom(
            native
        )
        result["rsmith"]["combined_energy_kj_mol_per_atom"] = rsmith_energy_per_atom(
            rsmith
        )
    return result


def sustained_time(trace, method: str, target: float):
    for index, checkpoint in enumerate(trace):
        remaining = trace[index:]
        if all(
            row[method]["total_fit"]["rms_over_dynamic_range"] <= target
            for row in remaining
        ):
            return (
                checkpoint[method]["wall_seconds"],
                checkpoint["moves"],
                checkpoint["epoch"],
            )
    return None


def summarise(case_root: Path, seeds, checkpoints):
    root = case_root / "results/independent-starts"
    traces = []
    for seed in seeds:
        trace = [score_checkpoint(case_root, seed, epoch) for epoch in checkpoints]
        traces.append({"seed": seed, "checkpoints": trace})

    summary = {
        "protocol": {
            "seeds": list(seeds),
            "replicates": len(seeds),
            "checkpoints": list(checkpoints),
            "moves_per_iteration": MOVES_PER_ITERATION,
            "start_generator": "periodic random-sequential hard spheres",
            "hard_sphere_minimum_angstrom": HARD_SPHERE_MIN_ANGSTROM,
            "first_shell_cutoff_angstrom": FIRST_SHELL_CUTOFF_ANGSTROM,
            "time_to_target_rule": "first observed checkpoint at or below target with every later checkpoint also at or below target",
            "target_residuals": list(TARGET_RESIDUALS),
        },
        "start_sha256": {
            str(seed): sha256(root / f"starts/seed-{seed}.ato") for seed in seeds
        },
        "traces": traces,
        "checkpoint_distributions": {},
        "time_to_target": {},
        "final_structural_distributions": {},
    }

    for epoch in checkpoints:
        rows = [
            trace["checkpoints"][list(checkpoints).index(epoch)] for trace in traces
        ]
        summary["checkpoint_distributions"][str(epoch)] = {}
        for method in ("native", "rsmith"):
            summary["checkpoint_distributions"][str(epoch)][method] = {
                "total_fit_rms_over_dynamic_range": distribution(
                    [row[method]["total_fit"]["rms_over_dynamic_range"] for row in rows]
                ),
                "wall_seconds": distribution(
                    [row[method]["wall_seconds"] for row in rows]
                ),
            }

    for target in TARGET_RESIDUALS:
        target_result = {}
        for method in ("native", "rsmith"):
            reached = [
                sustained_time(trace["checkpoints"], method, target) for trace in traces
            ]
            reached = [value for value in reached if value is not None]
            target_result[method] = {
                "successes": len(reached),
                "replicates": len(traces),
                "success_fraction": len(reached) / len(traces),
            }
            if reached:
                target_result[method]["wall_seconds"] = distribution(
                    [value[0] for value in reached]
                )
                target_result[method]["moves"] = distribution(
                    [value[1] for value in reached]
                )
                target_result[method]["epochs"] = distribution(
                    [value[2] for value in reached]
                )
        if target_result["native"]["successes"] == len(traces) and target_result[
            "rsmith"
        ]["successes"] == len(traces):
            target_result["rsmith_speedup_over_native_median"] = (
                target_result["native"]["wall_seconds"]["median"]
                / target_result["rsmith"]["wall_seconds"]["median"]
            )
        summary["time_to_target"][f"{target:.6f}"] = target_result

    final_rows = [trace["checkpoints"][-1] for trace in traces]
    for method in ("native", "rsmith"):
        summary["final_structural_distributions"][method] = {
            "rdf_peak_position_angstrom": distribution(
                [
                    row[method]["rdf"]["first_peak"]["position_angstrom"]
                    for row in final_rows
                ]
            ),
            "rdf_peak_height": distribution(
                [row[method]["rdf"]["first_peak"]["height"] for row in final_rows]
            ),
            "rdf_coordination_integral": distribution(
                [row[method]["rdf"]["coordination_integral"] for row in final_rows]
            ),
            "configuration_coordination_mean": distribution(
                [
                    row[method]["configuration"]["coordination_within_3_93_angstrom"][
                        "mean"
                    ]
                    for row in final_rows
                ]
            ),
            "nearest_neighbor_mean_angstrom": distribution(
                [
                    row[method]["configuration"]["nearest_neighbor_angstrom"]["mean"]
                    for row in final_rows
                ]
            ),
            "nearest_neighbor_min_angstrom": distribution(
                [
                    row[method]["configuration"]["nearest_neighbor_angstrom"]["min"]
                    for row in final_rows
                ]
            ),
            "close_contact_fraction": distribution(
                [
                    row[method]["configuration"][
                        "fraction_nearest_neighbor_below_2_3_angstrom"
                    ]
                    for row in final_rows
                ]
            ),
            "combined_energy_kj_mol_per_atom": distribution(
                [row[method]["combined_energy_kj_mol_per_atom"] for row in final_rows]
            ),
        }

    output = root / "summary.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "time_to_target": summary["time_to_target"],
                "final_structural_distributions": summary[
                    "final_structural_distributions"
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )
    print(f"Full independent-start summary: {output}")
    return summary


def prepare_starts(case_root: Path, seeds, force: bool):
    template = case_root / "reference/local/upstream/LiquidGa50C.ato"
    starts = case_root / "results/independent-starts/starts"
    starts.mkdir(parents=True, exist_ok=True)
    for seed in seeds:
        ato = starts / f"seed-{seed}.ato"
        xyz = starts / f"seed-{seed}.xyz"
        if ato.exists() and xyz.exists() and not force:
            continue
        positions, attempts = generate_hard_sphere_positions(
            seed, ATOM_COUNT, BOX_ANGSTROM, HARD_SPHERE_MIN_ANGSTROM
        )
        write_ato_with_positions(template, ato, positions, seed)
        write_xyz(xyz, positions, BOX_ANGSTROM)
        observed, _ = read_ato_positions(ato)
        minimum = configuration_metrics(observed, BOX_ANGSTROM)[
            "nearest_neighbor_angstrom"
        ]["min"]
        if minimum + 1.0e-8 < HARD_SPHERE_MIN_ANGSTROM:
            raise RuntimeError(
                f"generated seed {seed} violates the hard-sphere minimum"
            )
        (starts / f"seed-{seed}.json").write_text(
            json.dumps(
                {
                    "seed": seed,
                    "attempts": attempts,
                    "atom_count": ATOM_COUNT,
                    "box_angstrom": BOX_ANGSTROM,
                    "minimum_distance_angstrom": minimum,
                    "ato_sha256": sha256(ato),
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        print(f"seed {seed}: generated start in {attempts} attempts", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--seeds", default=",".join(str(seed) for seed in DEFAULT_SEEDS)
    )
    parser.add_argument(
        "--checkpoints", default=",".join(str(value) for value in DEFAULT_CHECKPOINTS)
    )
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--summarize-only", action="store_true")
    parser.add_argument("--binary", type=Path)
    args = parser.parse_args()
    seeds = tuple(int(value) for value in args.seeds.split(",") if value)
    checkpoints = tuple(
        sorted(set(int(value) for value in args.checkpoints.split(",")))
    )
    if not seeds or not checkpoints or checkpoints[0] != 0 or checkpoints[-1] != 40:
        raise SystemExit(
            "seeds must be nonempty and checkpoints must span 0 through 40"
        )
    if args.jobs <= 0:
        raise SystemExit("--jobs must be positive")
    if args.jobs != 1:
        print(
            "Warning: concurrent checkpoint timings are throughput-contaminated and must not be frozen",
            flush=True,
        )

    case_root = Path(__file__).resolve().parents[1]
    repo_root = case_root.parents[1]
    local = case_root / "reference/local"
    source = local / "upstream"
    record_file = local / "IMPORT.txt"
    if not source.is_dir() or not record_file.is_file():
        raise SystemExit(
            "run import_local_reference.sh --accept-local-testing-terms first"
        )
    record = dict(
        line.split("=", 1)
        for line in record_file.read_text().splitlines()
        if "=" in line
    )
    epsr_binary = Path(record["binary"])
    rsmith_binary = args.binary or repo_root / "target/release/rsmith"
    if not epsr_binary.is_file() or not rsmith_binary.is_file():
        raise SystemExit("native EPSR or release rsmith executable is missing")

    prepare_starts(case_root, seeds, args.force)

    def run_checkpoint(seed: int, epoch: int):
        if epoch == 0:
            print(f"seed {seed}: shared start scored", flush=True)
            return
        root = case_root / "results/independent-starts"
        start_ato = root / f"starts/seed-{seed}.ato"
        print(f"seed {seed}: checkpoint {epoch}, native EPSR", flush=True)
        run_native(
            source,
            epsr_binary,
            root / f"checkpoints/native/seed-{seed}/epoch-{epoch}",
            seed,
            epoch,
            args.force,
            start_ato=start_ato,
        )
        print(f"seed {seed}: checkpoint {epoch}, rsmith", flush=True)
        run_rsmith(
            start_ato,
            source / "LiquidGa50C.EPSR.q01",
            rsmith_binary,
            root / f"checkpoints/rsmith/seed-{seed}/epoch-{epoch}",
            seed,
            epoch,
            MOVES_PER_ITERATION,
            args.force,
        )

    if not args.summarize_only:
        tasks = [(seed, epoch) for seed in seeds for epoch in checkpoints]
        if args.jobs == 1:
            for seed, epoch in tasks:
                run_checkpoint(seed, epoch)
        else:
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=args.jobs
            ) as executor:
                futures = [
                    executor.submit(run_checkpoint, seed, epoch)
                    for seed, epoch in tasks
                ]
                for future in concurrent.futures.as_completed(futures):
                    future.result()

    summarise(case_root, seeds, checkpoints)


if __name__ == "__main__":
    main()
