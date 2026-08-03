#!/usr/bin/env python3
"""Prepare symmetric GAP<->Pedone synthetic recovery fixtures.

Raw structures remain local and ignored.  Targets are calculated independently
from minimum-image pair histograms, then converted to rsmith input conventions.
The opposite model is isotropically rescaled to the hidden target box without
relaxation or random displacement.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from collections import Counter
from pathlib import Path


TYPE_NAMES = {1: "Si", 2: "O"}
PAIR_TYPES = ((1, 1), (1, 2), (2, 2))
PAIR_LABELS = {(1, 1): "Si-Si", (1, 2): "Si-O", (2, 2): "O-O"}
NEUTRON_B = {1: 4.1491, 2: 5.803}
XRAY = {
    1: ([6.2915, 3.0353, 1.9891, 1.541], [2.4386, 32.3337, 0.6785, 81.6937], 1.1407),
    2: ([3.0485, 2.2868, 1.5463, 0.867], [13.2771, 5.7011, 0.3239, 32.9089], 0.2508),
}


def sha256(path: Path) -> str:
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
    origin = tuple(bounds[axis][0] for axis in "xyz")
    if len(atoms) != 3000 or set(bounds) != {"x", "y", "z"}:
        raise ValueError(f"unexpected structure content in {path}")
    return atoms, box, origin


def rescale(atoms, old_box, old_origin, new_box):
    scaled = []
    for atom_id, atom_type, position in atoms:
        fractional = tuple(
            (position[axis] - old_origin[axis]) / old_box[axis] for axis in range(3)
        )
        scaled.append(
            (
                atom_id,
                atom_type,
                tuple((fractional[axis] % 1.0) * new_box[axis] for axis in range(3)),
            )
        )
    return scaled


def write_lammps(path: Path, title: str, atoms, box):
    counts = Counter(atom_type for _, atom_type, _ in atoms)
    with path.open("w") as stream:
        stream.write(f"{title}\n\n{len(atoms)} atoms\n{len(counts)} atom types\n\n")
        for axis, length in zip("xyz", box):
            stream.write(f"0.0 {length:.15g} {axis}lo {axis}hi\n")
        stream.write("\nAtoms # charge\n\n")
        for atom_id, atom_type, position in atoms:
            charge = 2.4 if atom_type == 1 else -1.2
            stream.write(
                f"{atom_id} {atom_type} {charge:.1f} "
                + " ".join(f"{value:.15g}" for value in position)
                + "\n"
            )


def write_rmc6f(path: Path, title: str, atoms, box):
    counts = Counter(atom_type for _, atom_type, _ in atoms)
    volume = math.prod(box)
    with path.open("w") as stream:
        stream.write("(Version 6f format configuration file)\n")
        stream.write(f"Metadata title:     {title}\nMetadata material:  SiO2\n")
        stream.write(
            "Number of types of atoms:           2\nAtom types present:                 Si O\n"
        )
        stream.write(f"Number of each atom type:           {counts[1]} {counts[2]}\n")
        stream.write(
            "Number of moves generated:          0\nNumber of moves tried:              0\n"
        )
        stream.write(
            "Number of moves accepted:           0\nNumber of prior configuration saves: 0\n"
        )
        stream.write(f"Number of atoms:                     {len(atoms)}\n")
        stream.write(
            f"Number density (Ang^-3):             {len(atoms) / volume:.15g}\n"
        )
        stream.write("Supercell dimensions:                1 1 1\n")
        stream.write(
            "Cell (Ang/deg):    "
            + " ".join(f"{value:.15g}" for value in box)
            + " 90.0 90.0 90.0\n"
        )
        stream.write("Lattice vectors (Ang):\n")
        stream.write(
            f" {box[0]:.15g} 0.0 0.0\n 0.0 {box[1]:.15g} 0.0\n 0.0 0.0 {box[2]:.15g}\nAtoms:\n"
        )
        for atom_id, atom_type, position in atoms:
            fractional = [position[axis] / box[axis] for axis in range(3)]
            stream.write(
                f"{atom_id:6d} {TYPE_NAMES[atom_type]:>2s} [{atom_type}] "
                + " ".join(f"{value:.15g}" for value in fractional)
                + "\n"
            )


def displacement(first, second, box):
    return tuple(
        first[axis]
        - second[axis]
        - box[axis] * round((first[axis] - second[axis]) / box[axis])
        for axis in range(3)
    )


def partial_rdf(atoms, box, dr, rmax):
    counts = Counter(atom_type for _, atom_type, _ in atoms)
    positions = [position for _, _, position in atoms]
    types = [atom_type for _, atom_type, _ in atoms]
    nbins = int(round(rmax / dr))
    histograms = {pair: [0] * nbins for pair in PAIR_TYPES}
    for first in range(len(atoms) - 1):
        for second in range(first + 1, len(atoms)):
            vector = displacement(positions[second], positions[first], box)
            distance = math.sqrt(sum(value * value for value in vector))
            if distance < rmax:
                pair = tuple(sorted((types[first], types[second])))
                histograms[pair][int(distance / dr)] += 1
    volume = math.prod(box)
    radius = [(index + 0.5) * dr for index in range(nbins)]
    result = {}
    for pair in PAIR_TYPES:
        first_type, second_type = pair
        values = []
        for index, count in enumerate(histograms[pair]):
            inner, outer = index * dr, (index + 1) * dr
            shell = 4.0 * math.pi * (outer**3 - inner**3) / 3.0
            denominator = counts[first_type] * (counts[second_type] / volume) * shell
            factor = 2.0 if first_type == second_type else 1.0
            values.append(factor * count / denominator)
        result[pair] = values
    return radius, result


def partial_sq(radius, rdf, number_density, q_values, dr):
    result = {}
    prefactor = 4.0 * math.pi * number_density * dr
    for pair, values in rdf.items():
        integrand = [r * (g - 1.0) for r, g in zip(radius, values)]
        result[pair] = [
            1.0
            + prefactor
            * sum(value * math.sin(q * r) for value, r in zip(integrand, radius))
            / q
            for q in q_values
        ]
    return result


def xray_form_factor(atom_type, q):
    amplitudes, decays, constant = XRAY[atom_type]
    s2 = (q / (4.0 * math.pi)) ** 2
    return constant + sum(
        amplitude * math.exp(-decay * s2)
        for amplitude, decay in zip(amplitudes, decays)
    )


def totals(partials, q_values):
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
    neutron = [
        sum(neutron_weights[pair] * partials[pair][index] for pair in PAIR_TYPES)
        for index in range(len(q_values))
    ]
    xray = []
    for index, q in enumerate(q_values):
        form = {
            atom_type: xray_form_factor(atom_type, q) for atom_type in concentrations
        }
        mean_f = sum(concentrations[k] * form[k] for k in concentrations)
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
        )
    return neutron, xray


def write_curve(path, q_values, values, sigma):
    path.write_text(
        "# Q i(Q)=S(Q)-1 sigma\n"
        + "".join(
            f"{q:.12g} {value - 1.0:.12g} {sigma:.12g}\n"
            for q, value in zip(q_values, values)
        )
    )


def write_rsmith_config(
    path: Path,
    structure: Path,
    neutron: Path,
    xray: Path,
    output: Path,
    fit_mode: str,
    moves: int,
    seed: int,
    rmax: float,
    dr: float,
    qmin: float,
    qmax: float,
    qstep: float,
):
    datasets = []
    if fit_mode in {"neutron_only", "joint"}:
        datasets.append(
            f'''[data.neutron_sq]\nfile = "{neutron}"\nsigma_column = 3\nconvention = "iq"\nfit_min = {qmin}\nfit_max = {qmax}\nscattering_lengths = {{ Si = 4.1491, O = 5.803 }}\n'''
        )
    if fit_mode in {"xray_only", "joint"}:
        datasets.append(
            f'''[data.xray_sq]\nfile = "{xray}"\nsigma_column = 3\nconvention = "iq"\nfit_min = {qmin}\nfit_max = {qmax}\n'''
        )
    nq = int(round((qmax - qmin) / qstep))
    nbins = int(round(rmax / dr))
    path.write_text(
        f'''[system]\nstructure = "{structure}"\nformat = "lammps"\ntypes = {{ "1" = "Si", "2" = "O" }}\n\n[data]\n{"".join(datasets)}\n[rmc]\nmax_moves = {moves}\nmax_step = 0.1\ncheckpoint_every = 1000000000\nseed = {seed}\nprint_every = {max(1, moves)}\ntarget_acceptance = 0.25\nadjust_step_every = 1000\n\n[sq]\nqmin = {qmin}\nqmax = {qmax}\nnq = {nq}\nlorch = false\nrdf_cutoff = {rmax}\nrdf_nbins = {nbins}\n\n[constraints]\nmin_distance = {{ "Si-Si" = 2.0, "Si-O" = 1.35, "O-O" = 2.0 }}\n'''
    )


def append_rsmith_epsr(path: Path, tables: dict[str, Path], moves: int):
    with path.open("a") as stream:
        stream.write(
            f'''\n[potential]\ncutoff = 12.0\nweight = 1.0\n\n[[potential.tabulated]]\npair = "Si-Si"\nfile = "{tables["Si-Si"]}"\n\n[[potential.tabulated]]\npair = "Si-O"\nfile = "{tables["Si-O"]}"\n\n[[potential.tabulated]]\npair = "O-O"\nfile = "{tables["O-O"]}"\n\n[epsr]\nmode = "pure"\niterations = 1\nfeedback = 0.9\nsmooth_sigma = 0.03\nmoves_per_iteration = {moves}\ntemperature = 0.025852\nmin_r = 1.0\n'''
        )


def append_rsmith_pedone(path: Path):
    with path.open("a") as stream:
        stream.write(
            """\n[potential]\ncutoff = 15.0\nweight = 0.001\n\n[[potential.pedone]]\npair = "Si-Si"\nD0 = 0.0\nalpha = 1.0\nr0 = 1.0\nC0 = 0.0\n\n[[potential.pedone]]\npair = "Si-O"\nD0 = 0.340554\nalpha = 2.006700\nr0 = 2.100000\nC0 = 1.0\n\n[[potential.pedone]]\npair = "O-O"\nD0 = 0.042395\nalpha = 1.379316\nr0 = 3.618701\nC0 = 22.0\n\n[potential.coulomb]\nalpha = 0.25\ncharges = { Si = 2.4, O = -1.2 }\n"""
        )


def append_rsmith_gap(path: Path, model: Path):
    with path.open("a") as stream:
        stream.write(
            f'''\n[ml_potential]\nbackend = "gap_quip"\nmodel = "{model.resolve()}"\ninit_args = "Potential xml_label=GAP_2021_4_19_120_7_32_55_336"\nweight = 0.001\ncutoff = 5.0\ndelta = "local"\n'''
        )


parser = argparse.ArgumentParser()
parser.add_argument("--force", action="store_true")
parser.add_argument(
    "--moves", type=int, default=6000, help="smoke-stage attempted moves"
)
parser.add_argument("--seed", type=int, default=20260802)
parser.add_argument(
    "--gap-model",
    type=Path,
    default=None,
    help="Erhard silica GAP XML; when supplied, prepare the GAP-HRMC arm",
)
args = parser.parse_args()
if args.moves < 0:
    raise SystemExit("--moves must be non-negative")
if args.gap_model is not None and not args.gap_model.is_file():
    raise SystemExit(f"GAP model not found: {args.gap_model}")

case_root = Path(__file__).resolve().parents[1]
source = case_root / "reference/local/ambient-models"
sources = {"gap": source / "gap-300K.data", "pedone": source / "pedone-300K.data"}
for required in sources.values():
    if not required.is_file():
        raise SystemExit("run import_ambient_models.sh before preparing cross recovery")

output = case_root / "results/cross-recovery"
if output.exists():
    if not args.force:
        raise SystemExit(f"output exists: {output} (pass --force to replace it)")
else:
    output.mkdir(parents=True)

settings = {
    "qmin": 0.5,
    "qmax": 25.0,
    "qstep": 0.05,
    "rmax": 16.0,
    "dr": 0.02,
    "sigma": 0.02,
}
q_values = [
    settings["qmin"] + index * settings["qstep"]
    for index in range(
        int(round((settings["qmax"] - settings["qmin"]) / settings["qstep"]))
    )
]
loaded = {label: read_lammps_data(path) for label, path in sources.items()}
potential_root = case_root / "results/reference-potential-smoke"
potential_tables = {
    pair: potential_root / f"epsr-reference-{pair}.dat"
    for pair in ("Si-Si", "Si-O", "O-O")
}
have_reference_potential = all(path.is_file() for path in potential_tables.values())
summary = {
    "status": "cross_recovery_fixtures_prepared",
    "settings": settings,
    "seed": args.seed,
    "moves": args.moves,
    "rsmith_epsr_reference_potential_available": have_reference_potential,
    "rsmith_gap_model_available": args.gap_model is not None,
    "rsmith_gap_model_sha256": sha256(args.gap_model) if args.gap_model else None,
    "cases": {},
}

for target_label, start_label in (("gap", "pedone"), ("pedone", "gap")):
    case = output / f"target-{target_label}_from-{start_label}"
    case.mkdir(exist_ok=True)
    for old in case.glob("rsmith-*"):
        if old.is_dir():
            shutil.rmtree(old)
    target_atoms, target_box, target_origin = loaded[target_label]
    start_atoms, start_box, start_origin = loaded[start_label]
    target_atoms = rescale(target_atoms, target_box, target_origin, target_box)
    cross_atoms = rescale(start_atoms, start_box, start_origin, target_box)
    target_data = case / "hidden-target.data"
    cross_data = case / "cross-start.data"
    write_lammps(
        target_data, f"hidden {target_label} endpoint", target_atoms, target_box
    )
    write_lammps(
        cross_data,
        f"{start_label} endpoint rescaled to {target_label} density",
        cross_atoms,
        target_box,
    )
    write_rmc6f(
        case / "hidden-target.rmc6f",
        f"hidden {target_label} endpoint",
        target_atoms,
        target_box,
    )
    write_rmc6f(
        case / "cross-start.rmc6f",
        f"{start_label} endpoint rescaled to {target_label} density",
        cross_atoms,
        target_box,
    )

    radius, rdf = partial_rdf(
        target_atoms, target_box, settings["dr"], settings["rmax"]
    )
    partials = partial_sq(
        radius, rdf, len(target_atoms) / math.prod(target_box), q_values, settings["dr"]
    )
    neutron_total, xray_total = totals(partials, q_values)
    neutron_path, xray_path = (
        case / "target-neutron-iq.dat",
        case / "target-xray-iq.dat",
    )
    write_curve(neutron_path, q_values, neutron_total, settings["sigma"])
    write_curve(xray_path, q_values, xray_total, settings["sigma"])
    partial_path = case / "hidden-partial-sq.dat"
    partial_path.write_text(
        "# Q S_SiSi S_SiO S_OO\n"
        + "".join(
            f"{q:.12g} "
            + " ".join(f"{partials[pair][index]:.12g}" for pair in PAIR_TYPES)
            + "\n"
            for index, q in enumerate(q_values)
        )
    )

    for fit_mode in ("neutron_only", "xray_only", "joint"):
        run_dir = case / f"rsmith-{fit_mode}"
        run_dir.mkdir()
        write_rsmith_config(
            run_dir / "run.toml",
            cross_data,
            neutron_path,
            xray_path,
            run_dir,
            fit_mode,
            args.moves,
            args.seed,
            settings["rmax"],
            settings["dr"],
            settings["qmin"],
            settings["qmax"],
            settings["qstep"],
        )
    if have_reference_potential:
        epsr_dir = case / "rsmith-epsr-joint"
        epsr_dir.mkdir()
        write_rsmith_config(
            epsr_dir / "run.toml",
            cross_data,
            neutron_path,
            xray_path,
            epsr_dir,
            "joint",
            args.moves,
            args.seed,
            settings["rmax"],
            settings["dr"],
            settings["qmin"],
            settings["qmax"],
            settings["qstep"],
        )
        append_rsmith_epsr(epsr_dir / "run.toml", potential_tables, args.moves)
    pedone_dir = case / "rsmith-pedone-joint"
    pedone_dir.mkdir()
    write_rsmith_config(
        pedone_dir / "run.toml",
        cross_data,
        neutron_path,
        xray_path,
        pedone_dir,
        "joint",
        args.moves,
        args.seed,
        settings["rmax"],
        settings["dr"],
        settings["qmin"],
        settings["qmax"],
        settings["qstep"],
    )
    append_rsmith_pedone(pedone_dir / "run.toml")
    if args.gap_model is not None:
        gap_dir = case / "rsmith-gap-joint"
        gap_dir.mkdir()
        write_rsmith_config(
            gap_dir / "run.toml",
            cross_data,
            neutron_path,
            xray_path,
            gap_dir,
            "joint",
            args.moves,
            args.seed,
            settings["rmax"],
            settings["dr"],
            settings["qmin"],
            settings["qmax"],
            settings["qstep"],
        )
        append_rsmith_gap(gap_dir / "run.toml", args.gap_model)
    zero_dir = case / "rsmith-hidden-zero-move"
    zero_dir.mkdir()
    write_rsmith_config(
        zero_dir / "run.toml",
        target_data,
        neutron_path,
        xray_path,
        zero_dir,
        "joint",
        0,
        args.seed,
        settings["rmax"],
        settings["dr"],
        settings["qmin"],
        settings["qmax"],
        settings["qstep"],
    )
    cross_zero_dir = case / "rsmith-cross-zero-move"
    cross_zero_dir.mkdir()
    write_rsmith_config(
        cross_zero_dir / "run.toml",
        cross_data,
        neutron_path,
        xray_path,
        cross_zero_dir,
        "joint",
        0,
        args.seed,
        settings["rmax"],
        settings["dr"],
        settings["qmin"],
        settings["qmax"],
        settings["qstep"],
    )

    scale = tuple(target_box[index] / start_box[index] for index in range(3))
    summary["cases"][target_label] = {
        "target_source": str(sources[target_label]),
        "target_source_sha256": sha256(sources[target_label]),
        "start_source": str(sources[start_label]),
        "start_source_sha256": sha256(sources[start_label]),
        "target_box_a": target_box,
        "start_native_box_a": start_box,
        "cross_start_scale": scale,
        "target_number_density_atoms_a3": len(target_atoms) / math.prod(target_box),
        "cross_start_number_density_atoms_a3": len(cross_atoms) / math.prod(target_box),
        "files": {
            "hidden_target_sha256": sha256(target_data),
            "cross_start_sha256": sha256(cross_data),
            "neutron_target_sha256": sha256(neutron_path),
            "xray_target_sha256": sha256(xray_path),
        },
    }

(output / "fixture-summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(summary, indent=2, sort_keys=True))
print(f"Cross-recovery fixtures: {output}")
