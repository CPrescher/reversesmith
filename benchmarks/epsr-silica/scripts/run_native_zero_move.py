#!/usr/bin/env python3
"""Create a one-configuration, zero-move EPSR26 DTBsilicaNX calculation."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path


def replace_setting(text: str, name: str, value: str) -> str:
    pattern = re.compile(rf"^(\s*{re.escape(name)}\s+)([^\r\n]*)(\r?$)", re.MULTILINE)
    match = pattern.search(text)
    if match is None:
        raise ValueError(f"setting {name!r} not found in EPSR input")
    old_value_and_comment = match.group(2)
    separator = "               "
    comment = ""
    if separator in old_value_and_comment:
        comment = separator + old_value_and_comment.split(separator, 1)[1]
    replacement = match.group(1) + value + comment + match.group(3)
    return text[: match.start()] + replacement + text[match.end() :]


parser = argparse.ArgumentParser()
parser.add_argument("--epsr-root", type=Path)
parser.add_argument("--force", action="store_true")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
local = case_root / "reference/local"
source = local / "upstream"
record_path = local / "IMPORT.txt"
if not source.is_dir() or not record_path.is_file():
    raise SystemExit("run import_local_reference.sh --accept-local-testing-terms first")

record = dict(
    line.split("=", 1) for line in record_path.read_text().splitlines() if "=" in line
)
epsr_root = args.epsr_root or Path(record["epsr_root"])
binary = epsr_root / "bin/epsr"
if not binary.is_file():
    raise SystemExit(f"EPSR executable not found: {binary}")

run_dir = case_root / "results/native-zero-move"
if run_dir.exists():
    if not args.force:
        raise SystemExit(
            f"native run already exists: {run_dir} (pass --force to replace it)"
        )
    shutil.rmtree(run_dir)
shutil.copytree(source, run_dir)

input_path = run_dir / "DTBsilica.EPSR.inp"
text = input_path.read_text()
for name, value in (
    ("potfac", "0.0 0.0"),
    ("ireset", "1"),
    ("iinit", "1"),
    ("ntimes", "0"),
    ("niter", "1"),
    ("nsumt", "0"),
    ("ngrsamples", "6000"),
    ("fwhm", "0.0"),
    ("fwhmq", "0.0"),
):
    text = replace_setting(text, name, value)
input_path.write_text(text)

for path in run_dir.glob("DTBsilica.EPSR.*"):
    if path.suffix not in {".inp", ".inpa"}:
        path.unlink()

completed = subprocess.run(
    [str(binary), str(run_dir) + "/", "epsr", "DTBsilica"],
    cwd=run_dir,
    text=True,
    input="",
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    check=False,
)
(run_dir / "native-zero-move.log").write_text(completed.stdout)
if completed.returncode != 0:
    raise SystemExit(
        f"EPSR failed with exit code {completed.returncode}; see native-zero-move.log"
    )

required = [
    run_dir / f"DTBsilica.EPSR.{extension}" for extension in ("f01", "g01", "u01")
]
missing = [str(path) for path in required if not path.is_file()]
if missing:
    raise SystemExit("EPSR did not create required outputs: " + ", ".join(missing))
report = (run_dir / "DTBsilica.EPSR.out").read_text()
if (
    "Trial no.     0" not in report
    or "No. of configurations in sum     1" not in report
):
    raise SystemExit(
        "native output is not the requested one-configuration, zero-move calculation"
    )
print(f"Native zero-move EPSR output: {run_dir}")
