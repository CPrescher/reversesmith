#!/usr/bin/env python3
"""Independently characterize the local GAP and Pedone SiO2 endpoints."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, deque
from itertools import combinations
from pathlib import Path


TYPE_NAMES = {1: "Si", 2: "O"}
PAIR_TYPES = ((1, 1), (1, 2), (2, 2))
PAIR_LABELS = {(1, 1): "Si-Si", (1, 2): "Si-O", (2, 2): "O-O"}
PEAK_WINDOWS = {(1, 1): (2.5, 3.8), (1, 2): (1.3, 2.2), (2, 2): (2.0, 3.2)}
AVOGADRO = 6.02214076e23
MASS = {1: 28.0855, 2: 15.9994}


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_lammps_data(path: Path):
    lines = path.read_text(encoding="ascii", errors="replace").splitlines()
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
        atom_id, atom_type = int(fields[0]), int(fields[1])
        offset = 3 if atom_style == "charge" else 2
        position = tuple(float(value) for value in fields[offset : offset + 3])
        atoms.append((atom_id, atom_type, position))
    atoms.sort()
    if len(atoms) != 3000 or set(bounds) != {"x", "y", "z"}:
        raise ValueError(f"unexpected structure content in {path}")
    box = tuple(bounds[axis][1] - bounds[axis][0] for axis in "xyz")
    return atoms, box


def displacement(first, second, box):
    return tuple(
        first[axis] - second[axis] - box[axis] * round((first[axis] - second[axis]) / box[axis])
        for axis in range(3)
    )


def norm(vector):
    return math.sqrt(sum(value * value for value in vector))


def distribution(values):
    ordered = sorted(values)
    if not ordered:
        return {"count": 0}

    def quantile(fraction):
        coordinate = fraction * (len(ordered) - 1)
        low = int(math.floor(coordinate))
        high = int(math.ceil(coordinate))
        weight = coordinate - low
        return ordered[low] * (1.0 - weight) + ordered[high] * weight

    mean = sum(ordered) / len(ordered)
    variance = sum((value - mean) ** 2 for value in ordered) / len(ordered)
    return {
        "count": len(ordered),
        "mean": mean,
        "standard_deviation": math.sqrt(variance),
        "q05": quantile(0.05),
        "median": quantile(0.5),
        "q95": quantile(0.95),
    }


def angle_degrees(first, second):
    denominator = norm(first) * norm(second)
    cosine = sum(a * b for a, b in zip(first, second)) / denominator
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def shortest_alternate_path(adjacency, start, target, forbidden):
    queue = deque([(start, 0)])
    visited = {start}
    while queue:
        node, depth = queue.popleft()
        for neighbor in adjacency[node]:
            if (node == forbidden[0] and neighbor == forbidden[1]) or (
                node == forbidden[1] and neighbor == forbidden[0]
            ):
                continue
            if neighbor == target:
                return depth + 1
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, depth + 1))
    return None


def analyze(path: Path, dr: float, rmax: float, bond_cutoff: float):
    atoms, box = read_lammps_data(path)
    volume = math.prod(box)
    counts = Counter(atom_type for _, atom_type, _ in atoms)
    positions = [record[2] for record in atoms]
    types = [record[1] for record in atoms]
    nbins = int(round(rmax / dr))
    histograms = {pair: [0] * nbins for pair in PAIR_TYPES}
    minimum_distance = {pair: math.inf for pair in PAIR_TYPES}
    si_to_o = {index: [] for index, atom_type in enumerate(types) if atom_type == 1}
    o_to_si = {index: [] for index, atom_type in enumerate(types) if atom_type == 2}

    for first in range(len(atoms) - 1):
        for second in range(first + 1, len(atoms)):
            vector = displacement(positions[second], positions[first], box)
            distance = norm(vector)
            pair = tuple(sorted((types[first], types[second])))
            minimum_distance[pair] = min(minimum_distance[pair], distance)
            if distance < rmax:
                histograms[pair][int(distance / dr)] += 1
            if pair == (1, 2) and distance < bond_cutoff:
                si, oxygen = (first, second) if types[first] == 1 else (second, first)
                si_to_o[si].append(oxygen)
                o_to_si[oxygen].append(si)

    rdf = {}
    for pair in PAIR_TYPES:
        first_type, second_type = pair
        values = []
        for index, count in enumerate(histograms[pair]):
            inner, outer = index * dr, (index + 1) * dr
            radius = 0.5 * (inner + outer)
            shell = 4.0 * math.pi * (outer**3 - inner**3) / 3.0
            denominator = counts[first_type] * (counts[second_type] / volume) * shell
            factor = 2.0 if first_type == second_type else 1.0
            values.append((radius, factor * count / denominator if denominator else 0.0))
        rdf[pair] = values

    peak = {}
    for pair in PAIR_TYPES:
        lower, upper = PEAK_WINDOWS[pair]
        selected = [point for point in rdf[pair] if lower <= point[0] <= upper]
        radius, value = max(selected, key=lambda point: point[1])
        peak[PAIR_LABELS[pair]] = {"r_a": radius, "g": value}

    osi_angles = []
    for silicon, neighbors in si_to_o.items():
        vectors = [displacement(positions[oxygen], positions[silicon], box) for oxygen in neighbors]
        osi_angles.extend(angle_degrees(a, b) for a, b in combinations(vectors, 2))
    sio_angles = []
    for oxygen, neighbors in o_to_si.items():
        vectors = [displacement(positions[silicon], positions[oxygen], box) for silicon in neighbors]
        sio_angles.extend(angle_degrees(a, b) for a, b in combinations(vectors, 2))

    si_graph = {index: set() for index in si_to_o}
    for neighbors in o_to_si.values():
        for first, second in combinations(neighbors, 2):
            si_graph[first].add(second)
            si_graph[second].add(first)
    edges = [(first, second) for first in si_graph for second in si_graph[first] if first < second]
    edge_cycle_sizes = []
    acyclic_edges = 0
    for edge in edges:
        alternate = shortest_alternate_path(si_graph, edge[0], edge[1], edge)
        if alternate is None:
            acyclic_edges += 1
        else:
            edge_cycle_sizes.append(alternate + 1)

    total_mass_g = sum(counts[t] * MASS[t] for t in counts) / AVOGADRO
    density = total_mass_g / (volume * 1.0e-24)
    si_coordination = Counter(len(neighbors) for neighbors in si_to_o.values())
    o_coordination = Counter(len(neighbors) for neighbors in o_to_si.values())
    return {
        "path": str(path),
        "sha256": sha256(path),
        "atom_count": len(atoms),
        "counts": {TYPE_NAMES[key]: value for key, value in sorted(counts.items())},
        "box_a": box,
        "volume_a3": volume,
        "number_density_atoms_a3": len(atoms) / volume,
        "mass_density_g_cm3": density,
        "minimum_distance_a": {PAIR_LABELS[key]: value for key, value in minimum_distance.items()},
        "rdf_first_peak": peak,
        "bond_cutoff_a": bond_cutoff,
        "coordination": {
            "Si_by_O": {str(key): value for key, value in sorted(si_coordination.items())},
            "O_by_Si": {str(key): value for key, value in sorted(o_coordination.items())},
            "Si_defect_fraction": 1.0 - si_coordination[4] / counts[1],
            "O_defect_fraction": 1.0 - o_coordination[2] / counts[2],
        },
        "angles_degrees": {
            "O-Si-O": distribution(osi_angles),
            "Si-O-Si": distribution(sio_angles),
        },
        "shortest_ring_diagnostic": {
            "definition": "edge-weighted smallest alternate-path cycle in the Si network",
            "network_edges": len(edges),
            "acyclic_edges": acyclic_edges,
            "cycle_size_by_edge": {
                str(key): value for key, value in sorted(Counter(edge_cycle_sizes).items())
            },
        },
        "rdf": rdf,
    }


parser = argparse.ArgumentParser()
parser.add_argument("--force", action="store_true")
parser.add_argument("--dr", type=float, default=0.02)
parser.add_argument("--rmax", type=float, default=8.0)
parser.add_argument("--bond-cutoff", type=float, default=2.2)
args = parser.parse_args()
if args.dr <= 0 or args.rmax <= args.dr or args.bond_cutoff <= 0:
    raise SystemExit("invalid analysis grid or bond cutoff")

case_root = Path(__file__).resolve().parents[1]
source = case_root / "reference/local/ambient-models"
paths = {"gap": source / "gap-300K.data", "pedone": source / "pedone-300K.data"}
for path in paths.values():
    if not path.is_file():
        raise SystemExit("run import_ambient_models.sh with a local SiO2_glass checkout first")
output = case_root / "results/ambient-models"
if output.exists() and not args.force:
    raise SystemExit(f"output exists: {output} (pass --force to replace files)")
output.mkdir(parents=True, exist_ok=True)

analyses = {
    label: analyze(path, args.dr, args.rmax, args.bond_cutoff)
    for label, path in paths.items()
}
with (output / "partial-rdf.csv").open("w", newline="") as stream:
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(
        [
            "r_A",
            "gap_Si-Si",
            "pedone_Si-Si",
            "gap_Si-O",
            "pedone_Si-O",
            "gap_O-O",
            "pedone_O-O",
        ]
    )
    for index in range(int(round(args.rmax / args.dr))):
        row = [analyses["gap"]["rdf"][(1, 1)][index][0]]
        for pair in PAIR_TYPES:
            row.extend(
                (
                    analyses["gap"]["rdf"][pair][index][1],
                    analyses["pedone"]["rdf"][pair][index][1],
                )
            )
        writer.writerow(row)

upstream_rdf = source / "upstream-partial-rdf.csv"
with (output / "partial-rdf.csv").open(newline="") as generated_stream, upstream_rdf.open(
    newline=""
) as upstream_stream:
    generated_rows = list(csv.reader(generated_stream))[1:]
    upstream_rows = list(csv.reader(upstream_stream))[1:]
if len(generated_rows) != len(upstream_rows):
    raise SystemExit("generated RDF grid does not match the pinned upstream analysis")
rdf_max_absolute_difference = max(
    abs(float(generated) - float(upstream))
    for generated_row, upstream_row in zip(generated_rows, upstream_rows)
    for generated, upstream in zip(generated_row, upstream_row)
)

summary = {
    "status": "ambient_model_endpoints_characterized",
    "source_repository": "CPrescher/SiO2_glass",
    "source_commit": (source / "IMPORT.txt")
    .read_text()
    .split("commit=", 1)[1]
    .splitlines()[0],
    "scope": "single 300 K endpoint per model; not an equilibrium ensemble or ground truth",
    "analysis": {
        "rdf_dr_a": args.dr,
        "rdf_rmax_a": args.rmax,
        "si_o_bond_cutoff_a": args.bond_cutoff,
        "upstream_rdf_sha256": sha256(upstream_rdf),
        "rdf_max_absolute_difference_vs_upstream": rdf_max_absolute_difference,
    },
    "models": {},
}
for label, analysis in analyses.items():
    summary["models"][label] = {key: value for key, value in analysis.items() if key != "rdf"}
(output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(json.dumps(summary, indent=2, sort_keys=True))
print(f"Ambient-model analysis output: {output}")
