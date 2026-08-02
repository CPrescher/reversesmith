#!/usr/bin/env python3
"""Run deterministic thermal supercell replicas and average their PDFs."""

import argparse
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--replicas", type=int, default=10)
parser.add_argument("--first-seed", type=int, default=20260802)
parser.add_argument("--cells", type=int, default=12)
parser.add_argument("--binary", type=Path, default=Path("../../target/release/rsmith"))
args = parser.parse_args()
if args.replicas < 1:
    raise SystemExit("replicas must be positive")

root = Path.cwd()
generator = root / "scripts/generate_ni_supercell.py"
config = root / "rsmith/forward.toml"
curves = []
for replica in range(args.replicas):
    seed = args.first_seed + replica
    output = root / f"results/seed-{seed}"
    subprocess.run(
        [
            "python3",
            str(generator),
            "--cells",
            str(args.cells),
            "--lattice",
            "3.52",
            "--uiso",
            "0.00126651",
            "--seed",
            str(seed),
        ],
        check=True,
    )
    subprocess.run(
        [str(args.binary), str(config), "--output-dir", str(output), "--quiet"],
        check=True,
    )
    rows = []
    for raw in (output / "start_total_fr.dat").read_text().splitlines():
        if raw.startswith("#") or not raw.strip():
            continue
        r, value = map(float, raw.split()[:2])
        rows.append((r, value))
    curves.append(rows)

if any(len(curve) != len(curves[0]) for curve in curves):
    raise SystemExit("replica curves have different lengths")
averaged = []
for points in zip(*curves):
    r = points[0][0]
    if any(abs(point[0] - r) > 1e-12 for point in points[1:]):
        raise SystemExit("replica curves use different r grids")
    averaged.append((r, sum(point[1] for point in points) / len(points)))
target = root / "results/ensemble_mean_fr.dat"
target.parent.mkdir(parents=True, exist_ok=True)
target.write_text(
    f"# mean of {args.replicas} deterministic thermal-supercell replicas\n"
    + "".join(f"{r:.8f} {value:.12g}\n" for r, value in averaged)
)
print(f"wrote {len(averaged)} averaged points to {target}")
