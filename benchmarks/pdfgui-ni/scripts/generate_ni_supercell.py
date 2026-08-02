#!/usr/bin/env python3
"""Generate a thermally sampled fcc Ni supercell in LAMMPS data format."""

import argparse
import math
import random
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--cells", type=int, default=12)
parser.add_argument("--lattice", type=float, default=3.52)
parser.add_argument("--uiso", type=float, default=0.00126651)
parser.add_argument("--seed", type=int, default=20260802)
parser.add_argument("--output", type=Path, default=Path("rsmith/ni_12x12x12.data"))
args = parser.parse_args()
if args.cells <= 0 or args.lattice <= 0.0 or args.uiso < 0.0:
    raise SystemExit("cells and lattice must be positive; uiso must be non-negative")

basis = ((0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0))
box = args.cells * args.lattice
rng = random.Random(args.seed)
sigma = math.sqrt(args.uiso)
positions = []
for i in range(args.cells):
    for j in range(args.cells):
        for k in range(args.cells):
            for bx, by, bz in basis:
                fractional = (i + bx, j + by, k + bz)
                position = [
                    (args.lattice * fractional[axis] + rng.gauss(0.0, sigma)) % box
                    for axis in range(3)
                ]
                positions.append(position)

args.output.parent.mkdir(parents=True, exist_ok=True)
with args.output.open("w") as stream:
    stream.write("fcc Ni generated for the PDFgui benchmark\n\n")
    stream.write(f"{len(positions)} atoms\n1 atom types\n\n")
    stream.write(f"0.0 {box:.12f} xlo xhi\n0.0 {box:.12f} ylo yhi\n0.0 {box:.12f} zlo zhi\n\n")
    stream.write("Atoms # charge\n\n")
    for atom_id, (x, y, z) in enumerate(positions, 1):
        stream.write(f"{atom_id} 1 0.0 {x:.12f} {y:.12f} {z:.12f}\n")
print(f"wrote {len(positions)} atoms to {args.output}")
