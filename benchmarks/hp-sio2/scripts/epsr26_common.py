"""Local EPSR26 adapter shared by incremental high-pressure fixtures."""

from __future__ import annotations

import math
import re
import shutil
import subprocess
import time
from pathlib import Path

from pressure_common import read_lammps


def numeric_rows(path: Path, minimum_columns: int = 2):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        try:
            values = tuple(float(value) for value in line.split())
        except ValueError:
            continue
        if len(values) >= minimum_columns:
            rows.append(values)
    return rows


def replace_setting(text: str, name: str, value: str):
    pattern = re.compile(rf"^(\s*{re.escape(name)}\s+)([^\r\n]*)(\r?$)", re.MULTILINE)
    match = pattern.search(text)
    if match is None:
        raise ValueError(f"setting {name!r} not found in EPSR input")
    old = match.group(2)
    separator = "               "
    comment = separator + old.split(separator, 1)[1] if separator in old else ""
    replacement = match.group(1) + value + comment + match.group(3)
    return text[: match.start()] + replacement + text[match.end() :]


def write_ato(path: Path, structure: Path, template: Path, seed: int):
    positions, types, box, origin = read_lammps(structure)
    if len(positions) != 3000 or max(box) - min(box) > 1.0e-9:
        raise ValueError("EPSR comparison requires 3000 atoms in a cubic box")
    template_lines = template.read_text().splitlines()
    template_count = int(template_lines[0].split()[0])
    tail = template_lines[2 + 5 * template_count :]
    for index, line in enumerate(tail):
        fields = line.split()
        if len(fields) != 35:
            continue
        try:
            [int(value) for value in fields]
        except ValueError:
            continue
        tail[index] = " " + " ".join(
            str(value) for value in [0, -seed, 0, *([0] * 32)]
        )
        break
    with path.open("w") as stream:
        stream.write(f"{len(positions):8d} {box[0]:15.8E} {300.0:15.8E}\n")
        stream.write(template_lines[1] + "\n")
        for atom_id, (atom_type, position) in enumerate(zip(types, positions), 1):
            centered = position - origin - 0.5 * box
            stream.write(
                " 1 "
                + " ".join(f"{value:15.8E}" for value in centered)
                + "  0.00000000E+00  0.00000000E+00  0.00000000E+00"
                + f"  F      0    1.00000 {atom_id:5d}\n"
            )
            stream.write(f" {('Si' if atom_type == 1 else 'O'):<8s} 1      0\n")
            stream.write(
                "     0.00000000     0.00000000     0.00000000\n    0\n    0\n"
            )
        stream.write("\n".join(tail) + "\n")


def read_iq(path: Path):
    return [(row[0], row[1], row[2]) for row in numeric_rows(path, 3)]


def xray_form_factor(values, q):
    s2 = (q / (4.0 * math.pi)) ** 2
    return values[0] + sum(
        amplitude * math.exp(-decay * s2)
        for amplitude, decay in zip(values[1::2], values[2::2])
    )


def provisional_targets(input_root: Path, source: Path):
    neutron = read_iq(input_root / "target-neutron-iq.dat")
    xray = read_iq(input_root / "target-xray-iq.dat")
    neutron_factor = sum(
        row[2]
        for row in numeric_rows(source / "DTBsilica.NWTStot.wts", 3)
        if len(row) == 3
    )
    xparams = [
        row
        for row in numeric_rows(source / "DTBsilica.XWTS.wts", 11)
        if len(row) == 11
    ]
    neutron_native = [
        (q, value * neutron_factor, sigma * neutron_factor)
        for q, value, sigma in neutron
    ]
    xray_native = []
    for q, value, sigma in xray:
        f_si, f_o = xray_form_factor(xparams[0], q), xray_form_factor(xparams[1], q)
        mean_f = (f_si + 2.0 * f_o) / 3.0
        mean_f2 = (f_si * f_si + 2.0 * f_o * f_o) / 3.0
        factor = mean_f2 / (mean_f * mean_f)
        xray_native.append((q, value / factor, sigma / factor))
    return neutron_native, xray_native


def write_epsr_data(path: Path, title: str, rows):
    step = rows[-1][0] - rows[-2][0]
    extended = [*rows, (rows[-1][0] + step, rows[-1][1], rows[-1][2])]
    path.write_text(
        f"# {title}\n"
        + "".join(
            f"{q:.12g} {value:.12g} {sigma:.12g}\n"
            for q, value, sigma in extended
        )
    )


def configure_input(path: Path, moves: int, refinements: int, feedback: float):
    text = path.read_text()
    settings = (
        ("feedback", f"{feedback:.8f}"),
        ("potfac", "0.0 0.0" if moves == 0 else "1.0 1.0"),
        ("num_threds", "1"),
        ("nq", "500"),
        ("qstep", "0.05"),
        ("ireset", "1" if moves == 0 else "2"),
        ("iinit", "1"),
        ("ntimes", "0" if moves == 0 else str(max(1, moves // 3000))),
        ("niter", str(refinements)),
        ("nsumt", "0"),
        ("rho", "0.0"),
        ("cellst", "0.02"),
        ("rmaxgr", "16.0"),
        ("ngrsamples", "3000"),
        ("fwhm", "0.0"),
        ("fwhmq", "0.0"),
        ("nsmoop", "0"),
        ("fnameato", "Cross.ato"),
        ("qmin", "0.5 0.0"),
    )
    for name, value in settings:
        text = replace_setting(text, name, value)
    text = text.replace("SilicaGlassRT.mint01", "target-neutron.dat")
    text = text.replace("SiO2XRD.int01", "target-xray.dat")
    text = re.sub(r"(\bnrtype\s+)5", r"\g<1>6", text)
    path.write_text(text)


def prepare_run(
    source: Path,
    run: Path,
    structure: Path,
    moves: int,
    refinements: int,
    seed: int,
    targets,
    feedback: float,
):
    if run.exists():
        shutil.rmtree(run)
    shutil.copytree(source, run)
    for output in run.glob("DTBsilica.EPSR.*"):
        if output.suffix not in {".inp", ".inpa"}:
            output.unlink()
    write_ato(run / "Cross.ato", structure, source / "DTBsilica.ato", seed)
    write_epsr_data(run / "target-neutron.dat", "EPSR-native neutron i(Q)", targets[0])
    write_epsr_data(run / "target-xray.dat", "EPSR-native X-ray i(Q)", targets[1])
    configure_input(run / "DTBsilica.EPSR.inp", moves, refinements, feedback)


def run_epsr(binary: Path, run: Path):
    started = time.perf_counter()
    completed = subprocess.run(
        [str(binary), str(run) + "/", "epsr", "DTBsilica"],
        cwd=run,
        text=True,
        input="",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    wall = time.perf_counter() - started
    (run / "driver.log").write_text(completed.stdout)
    (run / "wall-seconds.txt").write_text(f"{wall:.9f}\n")
    if completed.returncode != 0 or not (run / "DTBsilica.EPSR.v01").is_file():
        raise RuntimeError(f"native EPSR failed in {run}; see driver.log")
    return wall


def calculated_targets(v01: Path, input_targets, sigma: float = 0.02):
    rows = [row for row in numeric_rows(v01, 5) if 0.5 <= row[0] < 25.0]
    if len(rows) != len(input_targets[0]) or len(rows) != len(input_targets[1]):
        raise ValueError("EPSR residual and target grids differ")
    return (
        [
            (row[0], input_targets[0][index][1] - row[1], sigma)
            for index, row in enumerate(rows)
        ],
        [
            (row[0], input_targets[1][index][1] - row[3], sigma)
            for index, row in enumerate(rows)
        ],
    )


def rms(first, second):
    residual = [a[1] - b[1] for a, b in zip(first, second)]
    return math.sqrt(sum(value * value for value in residual) / len(residual))
