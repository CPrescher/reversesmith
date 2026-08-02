#!/usr/bin/env python3
"""Generate a thermally sampled fluorite CaF2 supercell in LAMMPS format."""

from __future__ import annotations

import argparse
import math
import random
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--cells", type=int, default=15)
parser.add_argument("--lattice", type=float, default=5.462)
parser.add_argument("--uiso-ca", type=float, default=0.005)
parser.add_argument("--uiso-f", type=float, default=0.007)
parser.add_argument("--seed", type=int, default=20260802)
parser.add_argument("--output", type=Path, default=Path("rsmith/caf2_15x15x15.data"))
args = parser.parse_args()
if args.cells <= 0 or args.lattice <= 0.0 or args.uiso_ca < 0.0 or args.uiso_f < 0.0:
    raise SystemExit("cells and lattice must be positive; Uiso values must be non-negative")

ca_basis = (
    (0.0, 0.0, 0.0),
    (0.0, 0.5, 0.5),
    (0.5, 0.0, 0.5),
    (0.5, 0.5, 0.0),
)
f_basis = tuple(
    (x, y, z)
    for x in (0.25, 0.75)
    for y in (0.25, 0.75)
    for z in (0.25, 0.75)
)
box = args.cells * args.lattice
rng = random.Random(args.seed)
atoms = []
for i in range(args.cells):
    for j in range(args.cells):
        for k in range(args.cells):
            for type_id, basis, uiso in (
                (1, ca_basis, args.uiso_ca),
                (2, f_basis, args.uiso_f),
            ):
                sigma = math.sqrt(uiso)
                for bx, by, bz in basis:
                    ideal = (args.lattice * (i + bx), args.lattice * (j + by), args.lattice * (k + bz))
                    position = tuple((value + rng.gauss(0.0, sigma)) % box for value in ideal)
                    atoms.append((type_id, position))

args.output.parent.mkdir(parents=True, exist_ok=True)
with args.output.open("w") as stream:
    stream.write("fluorite CaF2 generated for the RMCProfile Qmax reproduction\n\n")
    stream.write(f"{len(atoms)} atoms\n2 atom types\n\n")
    stream.write(f"0.0 {box:.12f} xlo xhi\n0.0 {box:.12f} ylo yhi\n0.0 {box:.12f} zlo zhi\n\n")
    stream.write("Atoms # charge\n\n")
    for atom_id, (type_id, (x, y, z)) in enumerate(atoms, 1):
        stream.write(f"{atom_id} {type_id} 0.0 {x:.12f} {y:.12f} {z:.12f}\n")
print(f"wrote {len(atoms)} atoms in a {box:.3f} A box to {args.output}")
