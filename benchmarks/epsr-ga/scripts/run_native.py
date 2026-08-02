#!/usr/bin/env python3
"""Create a zero-move native EPSR26 forward calculation in an ignored sandbox."""

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
    comment_separator = "               "
    comment = ""
    if comment_separator in old_value_and_comment:
        comment = comment_separator + old_value_and_comment.split(comment_separator, 1)[1]
    replacement = match.group(1) + value + comment + match.group(3)
    return text[: match.start()] + replacement + text[match.end() :]


parser = argparse.ArgumentParser()
parser.add_argument("--epsr-root", type=Path)
parser.add_argument("--force", action="store_true", help="replace the previous ignored native run")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
local = case_root / "reference/local"
source = local / "upstream"
import_record = local / "IMPORT.txt"
if not source.is_dir() or not import_record.is_file():
    raise SystemExit("run import_local_reference.sh --accept-local-testing-terms first")

record = dict(
    line.split("=", 1) for line in import_record.read_text().splitlines() if "=" in line
)
epsr_root = args.epsr_root or Path(record["epsr_root"])
binary = epsr_root / "bin/epsr"
if not binary.is_file():
    raise SystemExit(f"EPSR executable not found: {binary}")

run_dir = case_root / "results/native-zero-move"
if run_dir.exists():
    if not args.force:
        raise SystemExit(f"native run already exists: {run_dir} (pass --force to replace it)")
    shutil.rmtree(run_dir)
shutil.copytree(source, run_dir)

input_path = run_dir / "LiquidGa50C.EPSR.inp"
text = input_path.read_text()
for name, value in (
    ("potfac", "0.0 0.0"),
    ("ireset", "1"),
    ("iinit", "1"),
    ("ntimes", "0"),
    ("niter", "1"),
    ("nsumt", "0"),
    ("ngrsamples", "5000"),
    ("fwhm", "0.0"),
    ("fwhmq", "0.0"),
):
    text = replace_setting(text, name, value)
input_path.write_text(text)

# Remove copied worked-example outputs so stale files cannot make a failed run
# look successful. Keep only EPSR inputs and the configuration/potential files.
for path in run_dir.glob("LiquidGa50C.EPSR.*"):
    if path.suffix not in {".inp", ".inpa"}:
        path.unlink()

completed = subprocess.run(
    [str(binary), str(run_dir) + "/", "epsr", "LiquidGa50C"],
    cwd=run_dir,
    text=True,
    input="",
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    check=False,
)
(run_dir / "native-zero-move.log").write_text(completed.stdout)
if completed.returncode != 0:
    raise SystemExit(f"EPSR failed with exit code {completed.returncode}; see native-zero-move.log")

required = [run_dir / f"LiquidGa50C.EPSR.{extension}" for extension in ("f01", "g01", "u01")]
missing = [str(path) for path in required if not path.is_file()]
if missing:
    raise SystemExit("EPSR did not create required outputs: " + ", ".join(missing))
native_report = (run_dir / "LiquidGa50C.EPSR.out").read_text()
if "Trial no.     0" not in native_report or "No. of configurations in sum     1" not in native_report:
    raise SystemExit("native output is not the requested one-configuration, zero-move calculation")
print(f"Native zero-move EPSR output: {run_dir}")
