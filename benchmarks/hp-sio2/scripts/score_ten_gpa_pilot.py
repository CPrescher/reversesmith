#!/usr/bin/env python3
"""Score the frozen 0-to-10 GPa pilot against held-out target structure."""

from __future__ import annotations

import argparse
import json
import math
import re
import tomllib
from collections import Counter
from pathlib import Path

import numpy as np


PAIR_TYPES = ((1, 1), (1, 2), (2, 2))
PAIR_LABELS = {(1, 1): "Si-Si", (1, 2): "Si-O", (2, 2): "O-O"}
SYMBOL_TYPE = {"Si": 1, "SI": 1, "si": 1, "O": 2, "o": 2}
NEUTRON_B = {1: 4.1491, 2: 5.803}
XRAY = {
    1: (
        [6.2915, 3.0353, 1.9891, 1.541],
        [2.4386, 32.3337, 0.6785, 81.6937],
        1.1407,
    ),
    2: (
        [3.0485, 2.2868, 1.5463, 0.867],
        [13.2771, 5.7011, 0.3239, 32.9089],
        0.2508,
    ),
}

parser = argparse.ArgumentParser()
parser.add_argument("--quiet", action="store_true")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
root = case_root / "results/ten-gpa-pilot"
protocol = tomllib.loads((case_root / "expected/ten-gpa-pilot.toml").read_text())
fixture = protocol["fixture"]


def read_lammps(path: Path):
    bounds, records = {}, []
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
        records.append(
            (
                int(fields[0]),
                int(fields[1]),
                [float(value) for value in fields[offset : offset + 3]],
            )
        )
    records.sort()
    return (
        np.asarray([record[2] for record in records]),
        np.asarray([record[1] for record in records], dtype=np.int8),
        np.asarray([bounds[axis][1] - bounds[axis][0] for axis in "xyz"]),
    )


def read_xyz(path: Path):
    lines = path.read_text().splitlines()
    count = int(lines[0])
    match = re.search(r'Lattice="([^"]+)"', lines[1])
    if match is None:
        raise ValueError(f"missing extended-XYZ lattice in {path}")
    lattice = [float(value) for value in match.group(1).split()]
    rows = [line.split() for line in lines[2 : 2 + count]]
    return (
        np.asarray([[float(value) for value in row[1:4]] for row in rows]),
        np.asarray([SYMBOL_TYPE[row[0]] for row in rows], dtype=np.int8),
        np.asarray([lattice[0], lattice[4], lattice[8]]),
    )


def read_epsr_ato(path: Path):
    lines = path.read_text().splitlines()
    count, length = int(lines[0].split()[0]), float(lines[0].split()[1])
    positions, types = [], []
    for index in range(count):
        atom = lines[2 + 5 * index].split()
        species = lines[3 + 5 * index].split()[0]
        positions.append(
            [(float(atom[axis]) + 0.5 * length) % length for axis in (1, 2, 3)]
        )
        types.append(SYMBOL_TYPE[species])
    return (
        np.asarray(positions),
        np.asarray(types, dtype=np.int8),
        np.asarray([length] * 3),
    )


def structure_metrics(path: Path):
    if path.suffix == ".xyz":
        positions, types, box = read_xyz(path)
    elif path.suffix == ".ato":
        positions, types, box = read_epsr_ato(path)
    else:
        positions, types, box = read_lammps(path)
    if len(positions) != protocol["source"]["atoms"]:
        raise ValueError(f"unexpected atom count in {path}: {len(positions)}")
    dr = fixture["rdf_bin_width_a"]
    rmax = fixture["rdf_cutoff_a"]
    cutoff = fixture["si_o_coordination_cutoff_a"]
    nbins = int(round(rmax / dr))
    histograms = {pair: np.zeros(nbins, dtype=np.int64) for pair in PAIR_TYPES}
    minimum = {pair: math.inf for pair in PAIR_TYPES}
    si_coord = np.zeros(int(np.sum(types == 1)), dtype=int)
    o_coord = np.zeros(int(np.sum(types == 2)), dtype=int)
    si_index = np.cumsum(types == 1) - 1
    o_index = np.cumsum(types == 2) - 1
    for first in range(len(positions) - 1):
        delta = positions[first + 1 :] - positions[first]
        delta -= box * np.rint(delta / box)
        distances = np.sqrt(np.einsum("ij,ij->i", delta, delta))
        second_types = types[first + 1 :]
        for other_type in (1, 2):
            pair = tuple(sorted((int(types[first]), other_type)))
            selected = distances[second_types == other_type]
            if len(selected):
                minimum[pair] = min(minimum[pair], float(np.min(selected)))
                bins = np.floor(selected[selected < rmax] / dr).astype(int)
                histograms[pair] += np.bincount(bins, minlength=nbins)
        unlike = np.flatnonzero(
            (second_types != types[first]) & (distances < cutoff)
        ) + first + 1
        if types[first] == 1:
            si_coord[si_index[first]] += len(unlike)
            for oxygen in unlike:
                o_coord[o_index[oxygen]] += 1
        else:
            o_coord[o_index[first]] += len(unlike)
            for silicon in unlike:
                si_coord[si_index[silicon]] += 1
    counts = Counter(int(value) for value in types)
    volume = float(np.prod(box))
    rdf = {}
    for pair in PAIR_TYPES:
        values = []
        for index, count in enumerate(histograms[pair]):
            inner, outer = index * dr, (index + 1) * dr
            shell = 4.0 * math.pi * (outer**3 - inner**3) / 3.0
            denominator = counts[pair[0]] * counts[pair[1]] / volume * shell
            factor = 2.0 if pair[0] == pair[1] else 1.0
            values.append(factor * count / denominator)
        rdf[PAIR_LABELS[pair]] = np.asarray(values)
    return {
        "box_a": box.tolist(),
        "minimum_distance_a": {
            PAIR_LABELS[pair]: minimum[pair] for pair in PAIR_TYPES
        },
        "si_coordination": dict(sorted(Counter(map(int, si_coord)).items())),
        "o_coordination": dict(sorted(Counter(map(int, o_coord)).items())),
        "mean_si_coordination": float(np.mean(si_coord)),
        "rdf": rdf,
    }


def coordination_tv(target, model):
    keys = set(target) | set(model)
    target_total = sum(target.values())
    model_total = sum(model.values())
    return 0.5 * sum(
        abs(target.get(key, 0) / target_total - model.get(key, 0) / model_total)
        for key in keys
    )


def xray_form_factor(atom_type, q):
    amplitudes, decays, constant = XRAY[atom_type]
    s2 = (q / (4.0 * math.pi)) ** 2
    return constant + sum(
        amplitude * math.exp(-decay * s2)
        for amplitude, decay in zip(amplitudes, decays)
    )


def reciprocal_curves(model, q_values):
    dr = fixture["rdf_bin_width_a"]
    radius = (np.arange(len(next(iter(model["rdf"].values())))) + 0.5) * dr
    density = protocol["source"]["atoms"] / np.prod(model["box_a"])
    partials = {}
    for pair, label in PAIR_LABELS.items():
        integrand = radius * (model["rdf"][label] - 1.0)
        partials[pair] = np.asarray(
            [
                1.0
                + 4.0
                * math.pi
                * density
                * dr
                * np.sum(integrand * np.sin(q * radius))
                / q
                for q in q_values
            ]
        )
    concentrations = {1: 1.0 / 3.0, 2: 2.0 / 3.0}
    mean_b = sum(concentrations[k] * NEUTRON_B[k] for k in concentrations)
    neutron = np.asarray(
        [
            sum(
                (1.0 if pair[0] == pair[1] else 2.0)
                * concentrations[pair[0]]
                * concentrations[pair[1]]
                * NEUTRON_B[pair[0]]
                * NEUTRON_B[pair[1]]
                * partials[pair][index]
                / mean_b**2
                for pair in PAIR_TYPES
            )
            - 1.0
            for index in range(len(q_values))
        ]
    )
    xray = []
    for index, q in enumerate(q_values):
        form = {kind: xray_form_factor(kind, q) for kind in concentrations}
        mean_f = sum(concentrations[kind] * form[kind] for kind in concentrations)
        xray.append(
            sum(
                (1.0 if pair[0] == pair[1] else 2.0)
                * concentrations[pair[0]]
                * concentrations[pair[1]]
                * form[pair[0]]
                * form[pair[1]]
                * partials[pair][index]
                / mean_f**2
                for pair in PAIR_TYPES
            )
            - 1.0
        )
    return neutron, np.asarray(xray)


def read_curve(path: Path):
    return np.asarray(
        [
            [float(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
    )


target_curves = {
    "neutron": read_curve(root / "inputs/target-neutron-iq.dat"),
    "xray": read_curve(root / "inputs/target-xray-iq.dat"),
}


def score_model(model, target):
    rdf_rms = {
        pair: float(np.sqrt(np.mean((model["rdf"][pair] - target["rdf"][pair]) ** 2)))
        for pair in target["rdf"]
    }
    q_values = target_curves["neutron"][:, 0]
    neutron, xray = reciprocal_curves(model, q_values)
    neutron_rms = float(
        np.sqrt(np.mean((neutron - target_curves["neutron"][:, 1]) ** 2))
    )
    xray_rms = float(np.sqrt(np.mean((xray - target_curves["xray"][:, 1]) ** 2)))
    return {
        "joint_iq_rms": math.sqrt((neutron_rms**2 + xray_rms**2) / 2.0),
        "mean_partial_rdf_rms": sum(rdf_rms.values()) / len(rdf_rms),
        "minimum_distance_a": model["minimum_distance_a"],
        "neutron_iq_rms": neutron_rms,
        "partial_rdf_rms": rdf_rms,
        "si_coordination_total_variation": coordination_tv(
            target["si_coordination"], model["si_coordination"]
        ),
        "mean_si_coordination": model["mean_si_coordination"],
        "mean_si_coordination_error": abs(
            model["mean_si_coordination"] - target["mean_si_coordination"]
        ),
        "xray_iq_rms": xray_rms,
    }


def run_metadata(run: Path):
    result = {
        "wall_seconds": float((run / "wall-seconds.txt").read_text())
        if (run / "wall-seconds.txt").is_file()
        else None
    }
    log_path = run / "rsmith.log"
    log = log_path.read_text() if log_path.is_file() else ""
    accepted = re.search(r"accepted (\d+)/(\d+)", log)
    if accepted:
        result["accepted_moves"] = int(accepted.group(1))
        result["attempted_moves"] = int(accepted.group(2))
    initial = re.search(
        r"Initial (?:pair potentials|GAP/QUIP|PACE/native) energy =\s*([+\-0-9.eE]+)",
        log,
    )
    final = re.findall(r"\[E:\s*([+\-0-9.eE]+)\]", log)
    if initial:
        result["initial_energy_ev"] = float(initial.group(1))
    if final:
        result["final_energy_ev"] = float(final[-1])
        if initial:
            result["energy_change_ev_per_atom"] = (
                result["final_energy_ev"] - result["initial_energy_ev"]
            ) / protocol["source"]["atoms"]
    return result


target = structure_metrics(root / "inputs/hidden-target.data")
start = score_model(structure_metrics(root / "inputs/cross-start.data"), target)
scores = {"start": start, "runs": {}, "target": {}}
scores["target"] = score_model(target, target)
for method in protocol["pilot"]["methods"]:
    scores["runs"][method] = {}
    for moves in protocol["pilot"]["checkpoints_moves"]:
        run = root / "runs" / method / f"moves-{int(moves):06d}"
        refined = run / "refined.xyz"
        if not refined.is_file():
            continue
        score = score_model(structure_metrics(refined), target)
        score.update(run_metadata(run))
        scores["runs"][method][str(moves)] = score

informative = (
    start["mean_partial_rdf_rms"] >= 0.02
    or start["si_coordination_total_variation"] >= 0.02
)
recoverable = any(
    score["joint_iq_rms"] < start["joint_iq_rms"]
    and score["mean_partial_rdf_rms"] <= 0.9 * start["mean_partial_rdf_rms"]
    for method in scores["runs"].values()
    for score in method.values()
)
pace_usable = False
for moves in map(str, protocol["pilot"]["checkpoints_moves"]):
    pure = scores["runs"].get("rsmith-rmc", {}).get(moves)
    pace = scores["runs"].get("rsmith-pace-w3", {}).get(moves)
    if pure and pace:
        pure_progress = start["joint_iq_rms"] - pure["joint_iq_rms"]
        pace_progress = start["joint_iq_rms"] - pace["joint_iq_rms"]
        if (
            pure_progress > 0.0
            and pace_progress >= 0.5 * pure_progress
            and pace["mean_partial_rdf_rms"] <= start["mean_partial_rdf_rms"]
        ):
            pace_usable = True

scores["decision"] = {
    "add_native_epsr_and_five_seed_replication": informative and recoverable,
    "informative_gap": informative,
    "pace_usable": pace_usable,
    "recoverability": recoverable,
}

comparison_protocol_path = case_root / "expected/ten-gpa-comparison.toml"
comparison_root = case_root / "results/ten-gpa-comparison"
if comparison_protocol_path.is_file() and comparison_root.is_dir():
    comparison_protocol = tomllib.loads(comparison_protocol_path.read_text())
    primary_seed = int(comparison_protocol["design"]["primary_seed"])
    endpoint = int(comparison_protocol["design"]["endpoint_moves"])
    seeds = [int(seed) for seed in comparison_protocol["design"]["seeds"]]
    comparison = {"primary_seed": {}, "endpoints": {}, "aggregate": {}}
    for method in comparison_protocol["design"]["methods"]:
        comparison["primary_seed"][method] = {}
        comparison["endpoints"][method] = {}
        for seed in seeds:
            checkpoints = (
                comparison_protocol["design"]["primary_seed_checkpoints_moves"]
                if seed == primary_seed
                else [endpoint]
            )
            for moves in checkpoints:
                if method == "rsmith-rmc" and seed == primary_seed:
                    run = root / "runs" / method / f"moves-{int(moves):06d}"
                    structure = run / "refined.xyz"
                else:
                    run = (
                        comparison_root
                        / "runs"
                        / method
                        / f"seed-{seed}"
                        / f"moves-{int(moves):06d}"
                    )
                    structure = (
                        run / "Cross.ato"
                        if method == "native-epsr26"
                        else run / "refined.xyz"
                    )
                if not structure.is_file():
                    continue
                score = score_model(structure_metrics(structure), target)
                score.update(run_metadata(run))
                if seed == primary_seed:
                    comparison["primary_seed"][method][str(moves)] = score
                if int(moves) == endpoint:
                    comparison["endpoints"][method][str(seed)] = score

    def distribution(values):
        array = np.asarray(values, dtype=float)
        return {
            "maximum": float(np.max(array)),
            "median": float(np.median(array)),
            "minimum": float(np.min(array)),
            "q1": float(np.quantile(array, 0.25)),
            "q3": float(np.quantile(array, 0.75)),
        }

    floors = comparison_protocol["metrics"]["safety_floor_a"]
    for method, endpoints in comparison["endpoints"].items():
        if len(endpoints) != len(seeds):
            continue
        comparison["aggregate"][method] = {
            metric: distribution([score[metric] for score in endpoints.values()])
            for metric in (
                "joint_iq_rms",
                "mean_partial_rdf_rms",
                "si_coordination_total_variation",
                "wall_seconds",
            )
        }
        comparison["aggregate"][method]["safety_passes"] = sum(
            all(
                score["minimum_distance_a"][pair] >= float(minimum)
                for pair, minimum in floors.items()
            )
            for score in endpoints.values()
        )

    pace = comparison["endpoints"].get("rsmith-pace-w30", {})
    epsr = comparison["endpoints"].get("native-epsr26", {})
    if len(pace) == len(seeds) and len(epsr) == len(seeds):
        wins = sum(
            pace[str(seed)]["mean_partial_rdf_rms"]
            < epsr[str(seed)]["mean_partial_rdf_rms"]
            for seed in seeds
        )
        pace_aggregate = comparison["aggregate"]["rsmith-pace-w30"]
        epsr_aggregate = comparison["aggregate"]["native-epsr26"]
        superior = (
            wins >= 4
            and pace_aggregate["mean_partial_rdf_rms"]["median"]
            < epsr_aggregate["mean_partial_rdf_rms"]["median"]
            and pace_aggregate["joint_iq_rms"]["median"]
            <= epsr_aggregate["joint_iq_rms"]["median"]
            and pace_aggregate["safety_passes"] == len(seeds)
        )
        comparison["decision"] = {
            "paired_hidden_rdf_wins": wins,
            "rsmith_pace_superior_to_epsr_for_frozen_task": superior,
            "safety_passes_required": len(seeds),
        }
    forward = comparison_root / "native-forward-summary.json"
    if forward.is_file():
        comparison["native_forward"] = json.loads(forward.read_text())
    scores["comparison"] = comparison

(root / "score-summary.json").write_text(
    json.dumps(scores, indent=2, sort_keys=True) + "\n"
)
if not args.quiet:
    print(json.dumps(scores, indent=2, sort_keys=True))
