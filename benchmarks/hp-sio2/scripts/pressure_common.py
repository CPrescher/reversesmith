"""Common structure and scattering routines for incremental pressure steps."""

from __future__ import annotations

import math
import re
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
        np.asarray([bounds[axis][0] for axis in "xyz"]),
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
        np.zeros(3),
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
        np.zeros(3),
    )


def read_structure(path: Path):
    if path.suffix == ".xyz":
        return read_xyz(path)
    if path.suffix == ".ato":
        return read_epsr_ato(path)
    return read_lammps(path)


def map_to_box(positions, old_box, old_origin, new_box):
    return ((positions - old_origin) / old_box % 1.0) * new_box


def write_lammps(path: Path, title: str, positions, types, box):
    path.write_text(
        f"{title}\n\n{len(positions)} atoms\n2 atom types\n\n"
        + "".join(
            f"0.0 {length:.15g} {axis}lo {axis}hi\n"
            for axis, length in zip("xyz", box)
        )
        + "\nMasses\n\n1 28.0855\n2 15.9994\n\nAtoms # charge\n\n"
        + "".join(
            f"{index + 1} {int(atom_type)} 0.0 "
            + " ".join(f"{value:.15g}" for value in position)
            + "\n"
            for index, (atom_type, position) in enumerate(zip(types, positions))
        )
    )


def structure_metrics(path: Path, fixture):
    positions, types, box, _ = read_structure(path)
    dr, rmax = fixture["rdf_bin_width_a"], fixture["rdf_cutoff_a"]
    cutoff = fixture["si_o_coordination_cutoff_a"]
    caps = fixture["short_range_caps_a"]
    quantiles = fixture["lower_tail_quantiles"]
    nbins = int(round(rmax / dr))
    histograms = {pair: np.zeros(nbins, dtype=np.int64) for pair in PAIR_TYPES}
    distances_by_pair = {pair: [] for pair in PAIR_TYPES}
    si_coord = np.zeros(int(np.sum(types == 1)), dtype=int)
    o_coord = np.zeros(int(np.sum(types == 2)), dtype=int)
    si_index, o_index = np.cumsum(types == 1) - 1, np.cumsum(types == 2) - 1
    for first in range(len(positions) - 1):
        delta = positions[first + 1 :] - positions[first]
        delta -= box * np.rint(delta / box)
        distances = np.sqrt(np.einsum("ij,ij->i", delta, delta))
        second_types = types[first + 1 :]
        for other_type in (1, 2):
            pair = tuple(sorted((int(types[first]), other_type)))
            selected = distances[second_types == other_type]
            bins = np.floor(selected[selected < rmax] / dr).astype(int)
            histograms[pair] += np.bincount(bins, minlength=nbins)
            short = selected[selected < float(caps[PAIR_LABELS[pair]])]
            if len(short):
                distances_by_pair[pair].append(short)
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
    counts, volume = Counter(map(int, types)), float(np.prod(box))
    rdf = {}
    for pair in PAIR_TYPES:
        values = []
        for index, count in enumerate(histograms[pair]):
            inner, outer = index * dr, (index + 1) * dr
            shell = 4.0 * math.pi * (outer**3 - inner**3) / 3.0
            denominator = counts[pair[0]] * counts[pair[1]] / volume * shell
            values.append((2.0 if pair[0] == pair[1] else 1.0) * count / denominator)
        rdf[PAIR_LABELS[pair]] = np.asarray(values)
    lower_tail = {}
    for pair in PAIR_TYPES:
        values = np.concatenate(distances_by_pair[pair])
        lower_tail[PAIR_LABELS[pair]] = {
            "count_within_cap": int(len(values)),
            "minimum_a": float(np.min(values)),
            "quantiles_a": {
                f"{float(level):.3f}": float(np.quantile(values, float(level)))
                for level in quantiles
            },
        }
    return {
        "_lower_tail_values": {
            PAIR_LABELS[pair]: np.concatenate(distances_by_pair[pair])
            for pair in PAIR_TYPES
        },
        "box_a": box.tolist(),
        "lower_tail": lower_tail,
        "mean_si_coordination": float(np.mean(si_coord)),
        "o_coordination": dict(sorted(Counter(map(int, o_coord)).items())),
        "rdf": rdf,
        "si_coordination": dict(sorted(Counter(map(int, si_coord)).items())),
    }


def xray_form_factor(atom_type, q):
    amplitudes, decays, constant = XRAY[atom_type]
    s2 = (q / (4.0 * math.pi)) ** 2
    return constant + sum(
        amplitude * math.exp(-decay * s2)
        for amplitude, decay in zip(amplitudes, decays)
    )


def reciprocal_curves(model, q_values, dr):
    radius = (np.arange(len(next(iter(model["rdf"].values())))) + 0.5) * dr
    density = 3000.0 / np.prod(model["box_a"])
    partials = {}
    for pair, label in PAIR_LABELS.items():
        integrand = radius * (model["rdf"][label] - 1.0)
        partials[pair] = np.asarray(
            [
                1.0
                + 4.0 * math.pi * density * dr * np.sum(integrand * np.sin(q * radius)) / q
                for q in q_values
            ]
        )
    concentrations = {1: 1.0 / 3.0, 2: 2.0 / 3.0}
    mean_b = sum(concentrations[k] * NEUTRON_B[k] for k in concentrations)
    neutron, xray = [], []
    for index, q in enumerate(q_values):
        neutron.append(
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
        )
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
    return partials, np.asarray(neutron), np.asarray(xray)


def write_curve(path: Path, q_values, values, sigma):
    path.write_text(
        "# Q i(Q)=S(Q)-1 sigma\n"
        + "".join(
            f"{q:.12g} {value:.12g} {sigma:.12g}\n"
            for q, value in zip(q_values, values)
        )
    )
