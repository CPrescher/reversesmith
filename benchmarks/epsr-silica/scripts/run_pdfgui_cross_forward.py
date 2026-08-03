#!/usr/bin/env python3
"""Use PDFgui's PDFfit2 engine as a forward-only cross-start PDF control."""

from __future__ import annotations

import json
import math
from pathlib import Path

try:
    from diffpy.pdffit2 import PdfFit
except ImportError as error:
    raise SystemExit(
        "run this script with the Python environment used by PDFgui (diffpy.pdffit2 required)"
    ) from error


def read_lammps(path: Path):
    lines = path.read_text().splitlines()
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
                int(fields[1]),
                tuple(float(value) for value in fields[offset : offset + 3]),
            )
        )
    box = tuple(bounds[axis][1] - bounds[axis][0] for axis in "xyz")
    return atoms, box


def write_stru(path: Path, title: str, atoms, box):
    with path.open("w") as stream:
        stream.write(f"title {title}\nformat pdffit\nscale 1.0\n")
        stream.write("sharp 0.0, 0.0, 1.0, 0.0\nspcgr P1\n")
        stream.write(
            "cell " + ", ".join(str(value) for value in (*box, 90.0, 90.0, 90.0)) + "\n"
        )
        stream.write("dcell 0.0, 0.0, 0.0, 0.0, 0.0, 0.0\n")
        stream.write(f"ncell 1, 1, 1, {len(atoms)}\natoms\n")
        for atom_type, position in atoms:
            symbol = "SI" if atom_type == 1 else "O"
            fractional = [position[axis] / box[axis] for axis in range(3)]
            stream.write(
                f"{symbol:<4s} "
                + " ".join(f"{value:.15g}" for value in fractional)
                + " 1.0\n"
            )
            # PDFfit2 needs a finite displacement width to represent discrete
            # atom-pair peaks on its real-space grid.  Apply the same small,
            # fixed Uiso to target and cross-start structures.
            stream.write(
                "0.0 0.0 0.0 0.0\n0.001 0.001 0.001\n0.0 0.0 0.0\n0.0 0.0 0.0\n0.0 0.0 0.0\n"
            )


def calculate(structure: Path, radiation: str):
    calculator = PdfFit(create_intro=False)
    calculator.read_struct(str(structure))
    calculator.alloc(radiation, 25.0, 0.0, 0.02, 16.0, 800)
    calculator.calc()
    return list(zip(calculator.getR(), calculator.getpdf_fit()))


def write_curve(path: Path, rows, radiation: str):
    path.write_text(
        f"# PDFfit2/PDFgui {radiation} G(r), Qmax=25 1/A\n"
        + "".join(f"{r:.8f} {value:.12g}\n" for r, value in rows)
    )


def metrics(reference, candidate):
    residual = [second[1] - first[1] for first, second in zip(reference, candidate)]
    dynamic_range = max(value for _, value in reference) - min(
        value for _, value in reference
    )
    rms = math.sqrt(sum(value * value for value in residual) / len(residual))
    return {
        "points": len(residual),
        "rms": rms,
        "rms_over_target_range": rms / dynamic_range,
        "max_absolute": max(abs(value) for value in residual),
    }


case_root = Path(__file__).resolve().parents[1]
fixture_root = case_root / "results/cross-recovery"
summary = {
    "program": "PDFgui/PDFfit2",
    "role": "forward-only; no atomistic RMC recovery",
    "uiso_a2": 0.001,
    "cases": {},
}
for case in sorted(fixture_root.glob("target-*_*")):
    output = case / "pdfgui-forward"
    output.mkdir(exist_ok=True)
    target_atoms, target_box = read_lammps(case / "hidden-target.data")
    start_atoms, start_box = read_lammps(case / "cross-start.data")
    write_stru(
        output / "hidden-target.stru",
        f"{case.name} hidden target",
        target_atoms,
        target_box,
    )
    write_stru(
        output / "cross-start.stru", f"{case.name} cross start", start_atoms, start_box
    )
    case_summary = {}
    for radiation in ("N", "X"):
        target = calculate(output / "hidden-target.stru", radiation)
        start = calculate(output / "cross-start.stru", radiation)
        write_curve(output / f"hidden-target-{radiation.lower()}.gr", target, radiation)
        write_curve(output / f"cross-start-{radiation.lower()}.gr", start, radiation)
        case_summary[radiation] = metrics(target, start)
    summary["cases"][case.name] = case_summary

(fixture_root / "pdfgui-forward-summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(summary, indent=2, sort_keys=True))
