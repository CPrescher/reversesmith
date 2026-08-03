#!/usr/bin/env python3
"""Score all available cross-recovery outputs with one independent analyzer."""

from __future__ import annotations

import csv
import json
import math
import re
from collections import Counter, deque
from itertools import combinations
from pathlib import Path

import numpy as np


PAIR_TYPES = ((1, 1), (1, 2), (2, 2))
PAIR_LABELS = {(1, 1): "Si-Si", (1, 2): "Si-O", (2, 2): "O-O"}
SYMBOL_TYPE = {"Si": 1, "SI": 1, "si": 1, "O": 2, "o": 2}
NEUTRON_B = {1: 4.1491, 2: 5.803}
XRAY = {
    1: ([6.2915, 3.0353, 1.9891, 1.541], [2.4386, 32.3337, 0.6785, 81.6937], 1.1407),
    2: ([3.0485, 2.2868, 1.5463, 0.867], [13.2771, 5.7011, 0.3239, 32.9089], 0.2508),
}


def read_lammps(path: Path):
    lines = path.read_text().splitlines()
    bounds, records = {}, []
    atom_style, in_atoms = None, False
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
        records.append(
            (
                int(fields[0]),
                int(fields[1]),
                [float(value) for value in fields[offset : offset + 3]],
            )
        )
    records.sort()
    box = np.asarray([bounds[axis][1] - bounds[axis][0] for axis in "xyz"])
    positions = np.asarray([record[2] for record in records])
    types = np.asarray([record[1] for record in records], dtype=np.int8)
    return positions, types, box


def read_xyz(path: Path):
    lines = path.read_text().splitlines()
    count = int(lines[0])
    match = re.search(r'Lattice="([^"]+)"', lines[1])
    if match is None:
        raise ValueError(f"missing extended-XYZ lattice in {path}")
    lattice = [float(value) for value in match.group(1).split()]
    box = np.asarray([lattice[0], lattice[4], lattice[8]])
    rows = [line.split() for line in lines[2 : 2 + count]]
    return (
        np.asarray([[float(value) for value in row[1:4]] for row in rows]),
        np.asarray([SYMBOL_TYPE[row[0]] for row in rows], dtype=np.int8),
        box,
    )


def read_rmc6f(path: Path):
    lines = path.read_text().splitlines()
    count = int(
        next(line for line in lines if line.startswith("Number of atoms:")).split()[-1]
    )
    cell = next(line for line in lines if line.startswith("Cell (Ang/deg):"))
    box = np.asarray([float(value) for value in cell.split(":", 1)[1].split()[:3]])
    start = lines.index("Atoms:") + 1
    rows = [line.split() for line in lines[start : start + count]]
    types = np.asarray([SYMBOL_TYPE[row[1]] for row in rows], dtype=np.int8)
    fractional = np.asarray([[float(value) for value in row[3:6]] for row in rows])
    return fractional * box, types, box


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


def read_structure(path: Path):
    if path.suffix == ".xyz":
        return read_xyz(path)
    if path.suffix == ".rmc6f":
        return read_rmc6f(path)
    if path.suffix == ".ato":
        return read_epsr_ato(path)
    return read_lammps(path)


def distribution(values):
    values = np.asarray(values, dtype=float)
    return {
        "count": int(len(values)),
        "mean": float(np.mean(values)),
        "standard_deviation": float(np.std(values)),
    }


def shortest_alternate_path(adjacency, start, target, forbidden):
    queue, visited = deque([(start, 0)]), {start}
    while queue:
        node, depth = queue.popleft()
        for neighbor in adjacency[node]:
            if {node, neighbor} == {forbidden[0], forbidden[1]}:
                continue
            if neighbor == target:
                return depth + 1
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, depth + 1))
    return None


def structure_metrics(path: Path, dr=0.02, rmax=16.0, bond_cutoff=2.2):
    positions, types, box = read_structure(path)
    if len(positions) != 3000:
        raise ValueError(f"unexpected atom count in {path}: {len(positions)}")
    nbins = int(round(rmax / dr))
    histograms = {pair: np.zeros(nbins, dtype=np.int64) for pair in PAIR_TYPES}
    minimum = {pair: math.inf for pair in PAIR_TYPES}
    si_to_o = {index: [] for index in np.flatnonzero(types == 1)}
    o_to_si = {index: [] for index in np.flatnonzero(types == 2)}
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
        if types[first] == 1:
            neighbors = (
                np.flatnonzero((second_types == 2) & (distances < bond_cutoff))
                + first
                + 1
            )
            for oxygen in neighbors:
                si_to_o[first].append(int(oxygen))
                o_to_si[int(oxygen)].append(first)
        else:
            neighbors = (
                np.flatnonzero((second_types == 1) & (distances < bond_cutoff))
                + first
                + 1
            )
            for silicon in neighbors:
                o_to_si[first].append(int(silicon))
                si_to_o[int(silicon)].append(first)
    counts = Counter(int(value) for value in types)
    volume = float(np.prod(box))
    rdf = {}
    for pair in PAIR_TYPES:
        values = []
        for index, count in enumerate(histograms[pair]):
            inner, outer = index * dr, (index + 1) * dr
            shell = 4.0 * math.pi * (outer**3 - inner**3) / 3.0
            denominator = counts[pair[0]] * counts[pair[1]] / volume * shell
            values.append((2.0 if pair[0] == pair[1] else 1.0) * count / denominator)
        rdf[PAIR_LABELS[pair]] = values

    def vector(center, neighbor):
        value = positions[neighbor] - positions[center]
        return value - box * np.rint(value / box)

    def angle(first, second):
        cosine = float(
            np.dot(first, second)
            / math.sqrt(np.dot(first, first) * np.dot(second, second))
        )
        return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))

    osi = [
        angle(vector(si, a), vector(si, b))
        for si, neighbors in si_to_o.items()
        for a, b in combinations(neighbors, 2)
    ]
    sio = [
        angle(vector(o, a), vector(o, b))
        for o, neighbors in o_to_si.items()
        for a, b in combinations(neighbors, 2)
    ]
    graph = {index: set() for index in si_to_o}
    for neighbors in o_to_si.values():
        for first, second in combinations(neighbors, 2):
            graph[first].add(second)
            graph[second].add(first)
    edges = [
        (first, second) for first in graph for second in graph[first] if first < second
    ]
    rings = Counter()
    for edge in edges:
        alternate = shortest_alternate_path(graph, edge[0], edge[1], edge)
        if alternate is not None:
            rings[alternate + 1] += 1
    si_coord, o_coord = (
        Counter(map(len, si_to_o.values())),
        Counter(map(len, o_to_si.values())),
    )
    return {
        "path": str(path),
        "box_a": box.tolist(),
        "number_density_atoms_a3": len(positions) / volume,
        "minimum_distance_a": {PAIR_LABELS[pair]: minimum[pair] for pair in PAIR_TYPES},
        "si_coordination": dict(sorted(si_coord.items())),
        "o_coordination": dict(sorted(o_coord.items())),
        "si_defect_fraction": 1.0 - si_coord[4] / counts[1],
        "o_defect_fraction": 1.0 - o_coord[2] / counts[2],
        "o_si_o_degrees": distribution(osi),
        "si_o_si_degrees": distribution(sio),
        "rings": dict(sorted(rings.items())),
        "rdf": rdf,
    }


def structural_distance(target, model):
    rdf_rms = {}
    for pair in target["rdf"]:
        residual = np.asarray(model["rdf"][pair]) - np.asarray(target["rdf"][pair])
        rdf_rms[pair] = float(np.sqrt(np.mean(residual * residual)))
    keys = set(target["rings"]) | set(model["rings"])
    target_total, model_total = (
        sum(target["rings"].values()),
        sum(model["rings"].values()),
    )
    ring_tv = 0.5 * sum(
        abs(
            target["rings"].get(key, 0) / target_total
            - model["rings"].get(key, 0) / model_total
        )
        for key in keys
    )
    return {
        "mean_partial_rdf_rms": sum(rdf_rms.values()) / len(rdf_rms),
        "partial_rdf_rms": rdf_rms,
        "si_defect_fraction_error": abs(
            model["si_defect_fraction"] - target["si_defect_fraction"]
        ),
        "o_defect_fraction_error": abs(
            model["o_defect_fraction"] - target["o_defect_fraction"]
        ),
        "o_si_o_mean_error_degrees": abs(
            model["o_si_o_degrees"]["mean"] - target["o_si_o_degrees"]["mean"]
        ),
        "si_o_si_mean_error_degrees": abs(
            model["si_o_si_degrees"]["mean"] - target["si_o_si_degrees"]["mean"]
        ),
        "ring_total_variation": ring_tv,
    }


def xray_form_factor(atom_type, q):
    amplitudes, decays, constant = XRAY[atom_type]
    s2 = (q / (4.0 * math.pi)) ** 2
    return constant + sum(
        amplitude * math.exp(-decay * s2)
        for amplitude, decay in zip(amplitudes, decays)
    )


def reciprocal_curves(model, q_values, dr=0.02):
    radius = (np.arange(len(next(iter(model["rdf"].values())))) + 0.5) * dr
    density = 3000.0 / np.prod(model["box_a"])
    partials = {}
    for pair, label in PAIR_LABELS.items():
        integrand = radius * (np.asarray(model["rdf"][label]) - 1.0)
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
    neutron_weights = {
        pair: (1.0 if pair[0] == pair[1] else 2.0)
        * concentrations[pair[0]]
        * concentrations[pair[1]]
        * NEUTRON_B[pair[0]]
        * NEUTRON_B[pair[1]]
        / mean_b**2
        for pair in PAIR_TYPES
    }
    neutron = np.asarray(
        [
            sum(neutron_weights[pair] * partials[pair][i] for pair in PAIR_TYPES) - 1.0
            for i in range(len(q_values))
        ]
    )
    xray = []
    for i, q in enumerate(q_values):
        form = {
            atom_type: xray_form_factor(atom_type, q) for atom_type in concentrations
        }
        mean_f = sum(concentrations[k] * form[k] for k in concentrations)
        total = sum(
            (1.0 if pair[0] == pair[1] else 2.0)
            * concentrations[pair[0]]
            * concentrations[pair[1]]
            * form[pair[0]]
            * form[pair[1]]
            * partials[pair][i]
            / mean_f**2
            for pair in PAIR_TYPES
        )
        xray.append(total - 1.0)
    return partials, neutron, np.asarray(xray)


def common_reciprocal_distance(case: Path, model):
    neutron_rows = read_space_curve(case / "target-neutron-iq.dat")
    xray_rows = read_space_curve(case / "target-xray-iq.dat")
    partial_rows = read_space_curve(case / "hidden-partial-sq.dat")
    q_values = np.asarray([row[0] for row in neutron_rows])
    partials, neutron, xray = reciprocal_curves(model, q_values)
    neutron_target = np.asarray([row[1] for row in neutron_rows])
    xray_target = np.asarray([row[1] for row in xray_rows])
    result = {
        "neutron_iq_rms": float(np.sqrt(np.mean((neutron - neutron_target) ** 2))),
        "xray_iq_rms": float(np.sqrt(np.mean((xray - xray_target) ** 2))),
        "partial_sq_rms": {},
    }
    for index, pair in enumerate(PAIR_TYPES):
        expected = np.asarray([row[index + 1] for row in partial_rows])
        result["partial_sq_rms"][PAIR_LABELS[pair]] = float(
            np.sqrt(np.mean((partials[pair] - expected) ** 2))
        )
    result["mean_partial_sq_rms"] = sum(result["partial_sq_rms"].values()) / len(
        PAIR_TYPES
    )
    return result


def read_space_curve(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            try:
                rows.append(tuple(float(value) for value in line.split()))
            except ValueError:
                pass
    return rows


def rms_column(rows, column):
    values = [row[column] for row in rows if 0.5 <= row[0] < 25.0]
    return math.sqrt(sum(value * value for value in values) / len(values))


def csv_fit_rms(path: Path):
    residual = []
    with path.open(newline="") as stream:
        for row in csv.reader(stream):
            try:
                residual.append(float(row[1]) - float(row[2]))
            except (ValueError, IndexError):
                pass
    return math.sqrt(sum(value * value for value in residual) / len(residual))


def rsmith_fit(run_dir: Path, target_root: Path | None = None):
    result = {
        "wall_seconds": float((run_dir / "wall-seconds.txt").read_text())
        if (run_dir / "wall-seconds.txt").is_file()
        else None
    }
    for stage in ("start", "refined"):
        for kind in ("neutron", "xray"):
            path = run_dir / f"{stage}_{kind}_sq.dat"
            target = (target_root or run_dir.parent) / f"target-{kind}-iq.dat"
            if path.is_file() and target.is_file():
                actual = {round(row[0], 8): row[1] for row in read_space_curve(path)}
                expected = read_space_curve(target)
                residual = [actual[round(row[0], 8)] - row[1] for row in expected]
                result[f"{stage}_{kind}_rms"] = math.sqrt(
                    sum(value * value for value in residual) / len(residual)
                )
    log = (
        (run_dir / "rsmith.log").read_text()
        if (run_dir / "rsmith.log").is_file()
        else ""
    )
    accepted = re.search(r"accepted (\d+)/(\d+)", log)
    if accepted:
        result["accepted_moves"] = int(accepted.group(1))
        result["attempted_moves"] = int(accepted.group(2))
    initial_energy = re.search(
        r"Initial (?:pair potentials|GAP/QUIP|PACE/native) energy =\s*([+\-0-9.eE]+) eV",
        log,
    )
    final_energies = re.findall(r"\[E:\s*([+\-0-9.eE]+)\]", log)
    suggested = re.search(
        r"Suggested weight \(chi2 ≈ energy influence\) =\s*([+\-0-9.eE]+)", log
    )
    if initial_energy:
        result["initial_energy_ev"] = float(initial_energy.group(1))
    if final_energies:
        result["final_energy_ev"] = float(final_energies[-1])
        if initial_energy:
            result["energy_change_ev_per_atom"] = (
                result["final_energy_ev"] - result["initial_energy_ev"]
            ) / 3000.0
    if suggested:
        result["suggested_balance_weight"] = float(suggested.group(1))
    return result


def rmcprofile_fit(run_dir: Path):
    result = {
        "wall_seconds": float((run_dir / "wall-seconds.txt").read_text())
        if (run_dir / "wall-seconds.txt").is_file()
        else None
    }
    for kind, filename in (("neutron", "cross_SQ1.csv"), ("xray", "cross_FQ1.csv")):
        path = run_dir / filename
        if path.is_file():
            result[f"{kind}_rms"] = csv_fit_rms(path)
    config = run_dir / "cross.rmc6f"
    if config.is_file():
        text = config.read_text()
        for label, key in (
            ("generated", "snapshot_generated_moves"),
            ("accepted", "snapshot_accepted_moves"),
        ):
            match = re.search(rf"Number of moves {label}:\s+(\d+)", text)
            if match:
                result[key] = int(match.group(1))
    log = (
        (run_dir / "driver.log").read_text()
        if (run_dir / "driver.log").is_file()
        else ""
    )
    generated = re.findall(r"Moves generated\s+::\s+(\d+)", log)
    accepted = re.findall(r"Moves accepted\s+::\s+(\d+)", log)
    if generated:
        result["completed_generated_moves"] = int(generated[-1])
    if accepted:
        result["completed_accepted_moves"] = int(accepted[-1])
    return result


def epsr_fit(run_dir: Path):
    result = {
        "wall_seconds": float((run_dir / "wall-seconds.txt").read_text())
        if (run_dir / "wall-seconds.txt").is_file()
        else None
    }
    residual = run_dir / "DTBsilica.EPSR.v01"
    if residual.is_file():
        rows = read_space_curve(residual)
        result["neutron_rms"] = rms_column(rows, 1)
        result["xray_rms"] = rms_column(rows, 3)
    return result


case_root = Path(__file__).resolve().parents[1]
fixture_root = case_root / "results/cross-recovery"
summary = {
    "status": "available_cross_recovery_outputs_scored",
    "scope": "smoke adapter validation, not production comparison",
    "cases": {},
}
for case in sorted(fixture_root.glob("target-*_*")):
    target = structure_metrics(case / "hidden-target.data")
    structure_paths = {"cross_start": case / "cross-start.data"}
    candidates = {
        "rsmith_joint": case / "rsmith-joint/refined.xyz",
        "rsmith_pedone_joint": case / "rsmith-pedone-joint/refined.xyz",
        "rsmith_gap_joint": case / "rsmith-gap-joint/refined.xyz",
        "rsmith_epsr_joint": case / "rsmith-epsr-joint/refined.xyz",
        "rmcprofile_joint": case / "rmcprofile-joint/cross.rmc6f",
        "native_epsr_joint": case / "native-epsr-joint/Cross.ato",
    }
    structure_paths.update(
        {key: path for key, path in candidates.items() if path.is_file()}
    )
    structures = {key: structure_metrics(path) for key, path in structure_paths.items()}
    fits = {}
    for name in (
        "rsmith-cross-zero-move",
        "rsmith-neutron_only",
        "rsmith-xray_only",
        "rsmith-joint",
        "rsmith-pedone-joint",
        "rsmith-gap-joint",
        "rsmith-epsr-joint",
    ):
        if (case / name).is_dir() and (case / name / "rsmith.log").is_file():
            fits[name] = rsmith_fit(case / name)
    for name in (
        "rmcprofile-cross-zero-move",
        "rmcprofile-neutron_only",
        "rmcprofile-xray_only",
        "rmcprofile-joint",
    ):
        if (case / name).is_dir() and (case / name / "driver.log").is_file():
            fits[name] = rmcprofile_fit(case / name)
    for name in ("native-epsr-cross-zero-move", "native-epsr-joint"):
        if (case / name).is_dir() and (case / name / "DTBsilica.EPSR.v01").is_file():
            fits[name] = epsr_fit(case / name)
    summary["cases"][case.name] = {
        "target_metrics": {key: value for key, value in target.items() if key != "rdf"},
        "fits": fits,
        "structural_distance_from_hidden_target": {
            key: structural_distance(target, value) for key, value in structures.items()
        },
        "common_reciprocal_distance_from_hidden_target": {
            "hidden_target_self_check": common_reciprocal_distance(case, target),
            **{
                key: common_reciprocal_distance(case, value)
                for key, value in structures.items()
            },
        },
        "model_metrics": {
            key: {metric: value for metric, value in model.items() if metric != "rdf"}
            for key, model in structures.items()
        },
    }

(fixture_root / "score-summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(summary, indent=2, sort_keys=True))


def score_hrmc_root(root: Path, status: str, scope: str):
    if not root.is_dir():
        return
    hrmc_summary = {"status": status, "scope": scope, "cases": {}}
    for hrmc_case in sorted(root.glob("target-*_*")):
        cross_case = fixture_root / hrmc_case.name
        target = structure_metrics(cross_case / "hidden-target.data")
        runs = {}
        for run in sorted(hrmc_case.iterdir()):
            refined = run / "refined.xyz"
            if run.is_dir() and refined.is_file():
                model = structure_metrics(refined)
                runs[run.name] = {
                    "fit": rsmith_fit(run, cross_case),
                    "common_reciprocal_distance_from_hidden_target": common_reciprocal_distance(
                        cross_case, model
                    ),
                    "structural_distance_from_hidden_target": structural_distance(
                        target, model
                    ),
                    "model_metrics": {
                        key: value for key, value in model.items() if key != "rdf"
                    },
                }
        hrmc_summary["cases"][hrmc_case.name] = runs
    (root / "score-summary.json").write_text(
        json.dumps(hrmc_summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(hrmc_summary, indent=2, sort_keys=True))


score_hrmc_root(
    case_root / "results/hrmc-weight-sweep",
    "available_hrmc_weight_pilot_outputs_scored",
    "one-third move per atom pilot, not production comparison",
)
score_hrmc_root(
    case_root / "results/hrmc-pace2024-weight-pilot",
    "available_pace2024_weight_pilot_outputs_scored",
    "one-third move per atom PACE pilot, not production comparison",
)
score_hrmc_root(
    case_root / "results/hrmc-production-bracket",
    "available_hrmc_production_bracket_outputs_scored",
    "two moves per atom, single-seed bracket; not the final multi-seed comparison",
)
score_hrmc_root(
    case_root / "results/hrmc-production-unscreened",
    "available_joint_acceptance_hrmc_production_outputs_scored",
    "two moves per atom, single-seed Pedone/GAP/PACE bracket; not the final multi-seed comparison",
)


def score_multiseed_root(root: Path):
    if not root.is_dir():
        return
    result = {
        "status": "available_multiseed_cross_program_outputs_scored",
        "scope": "ten-seed, 6000-move fixed-budget comparison with independent common scoring",
        "cases": {},
    }
    for ensemble_case in sorted(root.glob("target-*_*")):
        cross_case = fixture_root / ensemble_case.name
        target = structure_metrics(cross_case / "hidden-target.data")
        methods = {}
        for method in sorted(path for path in ensemble_case.iterdir() if path.is_dir()):
            seeds = {}
            for run in sorted(method.glob("seed-*")):
                if method.name == "native-epsr26":
                    structure = run / "Cross.ato"
                    fit = epsr_fit(run)
                elif method.name == "rmcprofile":
                    structure = run / "cross.rmc6f"
                    fit = rmcprofile_fit(run)
                    audit_path = run / "move-audit.json"
                    if audit_path.is_file():
                        audit = json.loads(audit_path.read_text())
                        fit["attempted_moves"] = audit["total_moves"]
                        fit["accepted_moves"] = audit["total_accepted"]
                        fit["sampling_wall_seconds"] = audit[
                            "sampling_wall_seconds"
                        ]
                        fit["exact_tail_wall_seconds"] = audit[
                            "exact_tail_wall_seconds"
                        ]
                else:
                    structure = run / "refined.xyz"
                    fit = rsmith_fit(run, cross_case)
                if not structure.is_file():
                    continue
                model = structure_metrics(structure)
                seeds[run.name] = {
                    "fit": fit,
                    "common_reciprocal_distance_from_hidden_target": common_reciprocal_distance(
                        cross_case, model
                    ),
                    "structural_distance_from_hidden_target": structural_distance(
                        target, model
                    ),
                    "model_metrics": {
                        key: value for key, value in model.items() if key != "rdf"
                    },
                }
            methods[method.name] = seeds
        result["cases"][ensemble_case.name] = methods
    (root / "raw-score-summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(result, indent=2, sort_keys=True))


score_multiseed_root(case_root / "results/multiseed-comparison")


def score_epsr_convergence_root(root: Path):
    if not root.is_dir():
        return
    result = {
        "status": "available_epsr_convergence_pilot_outputs_scored",
        "scope": "single-seed independent deterministic prefixes through 100 empirical-potential refinements",
        "cases": {},
    }
    for pilot_case in sorted(root.glob("target-*_*")):
        cross_case = fixture_root / pilot_case.name
        target = structure_metrics(cross_case / "hidden-target.data")
        methods = {}
        for method_path in sorted(path for path in pilot_case.iterdir() if path.is_dir()):
            method = method_path.name
            checkpoints = {}
            for run in sorted(method_path.glob("seed-*/iter-*")):
                if method == "native-epsr26":
                    structure = run / "Cross.ato"
                    fit = epsr_fit(run)
                else:
                    structure = run / "refined.xyz"
                    fit = rsmith_fit(run, cross_case)
                if not structure.is_file():
                    continue
                model = structure_metrics(structure)
                checkpoints[run.name] = {
                    "seed": run.parent.name,
                    "fit": fit,
                    "common_reciprocal_distance_from_hidden_target": common_reciprocal_distance(
                        cross_case, model
                    ),
                    "structural_distance_from_hidden_target": structural_distance(
                        target, model
                    ),
                    "model_metrics": {
                        key: value for key, value in model.items() if key != "rdf"
                    },
                }
            methods[method] = checkpoints
        result["cases"][pilot_case.name] = methods
    (root / "raw-score-summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(result, indent=2, sort_keys=True))


score_epsr_convergence_root(case_root / "results/epsr-convergence-pilot")
