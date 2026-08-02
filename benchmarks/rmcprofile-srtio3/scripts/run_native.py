#!/usr/bin/env python3
"""Generate native RMCProfile7 outputs from the official SrTiO3 tutorial."""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import tarfile
from pathlib import Path


EXPECTED_SHA256 = "828b9292167a6691e02d1a9073006a42be7300ce770309702a827e5886b105c3"
IMAGE = "arm64v8/ubuntu:24.04"


def docker(package: Path, run_dir: Path, executable: str, *arguments: str) -> None:
    subprocess.run(
        [
            "docker",
            "run",
            "--rm",
            "--platform",
            "linux/arm64",
            "-v",
            f"{package}:/opt/rmc",
            "-v",
            f"{run_dir}:/work",
            "-w",
            "/work",
            "-e",
            "LD_LIBRARY_PATH=/opt/rmc/exe/libs:/usr/lib/aarch64-linux-gnu",
            IMAGE,
            f"/opt/rmc/exe/{executable}",
            *arguments,
        ],
        check=True,
    )


parser = argparse.ArgumentParser()
parser.add_argument("--archive", type=Path, required=True)
args = parser.parse_args()
case_root = Path(__file__).resolve().parents[1]
native_root = case_root / "results/native"
package_root = native_root / "RMCProfile_package"
run_dir = native_root / "run"

digest = hashlib.sha256(args.archive.read_bytes()).hexdigest()
if digest != EXPECTED_SHA256:
    raise SystemExit(f"archive checksum mismatch: expected {EXPECTED_SHA256}, got {digest}")
if native_root.exists():
    shutil.rmtree(native_root)
native_root.mkdir(parents=True)
with tarfile.open(args.archive, "r:gz") as archive:
    archive.extractall(native_root, filter="data")
shutil.copytree(package_root / "tutorial/ex_4/01small", run_dir)

docker(package_root, run_dir, "RMCCreate", "-rmc7", "-supercell", "[1 1 1]", "SRTIO3_5K.rmc6f")

control = (run_dir / "SRTIO3_5K_RMC7.dat").read_text()
control = control.replace("PRINT_PERIOD :: 10000", "PRINT_PERIOD :: 1")
control = control.replace("SAVE_PERIOD  :: 10000", "SAVE_PERIOD  :: 0")
control = control.replace("TIME_LIMIT   :: 10 MINUTES", "TIME_LIMIT   :: 0 MINUTES")
control = control.replace("%%ITERATION_LIMIT :: 0", "ITERATION_LIMIT :: 0")
control = control.replace("SRTIO3_5K.rmc6f", "SRTIO3_5K.rmc7")
control = control.replace("> ATOMS :: Sr Ti O", "> ATOMS :: O Ti Sr")
control = control.replace("> END_POINT :: 2500", "> END_POINT :: 780")
minimum_start = control.index("MINIMUM_DISTANCES ::")
flags_start = control.index("FLAGS ::")
control = control[:minimum_start] + control[flags_start:]
control = control.replace("  > CSSR\n", "").replace("  > COORDINATION_NUMBER\n", "")
bragg_start = control.index("\nBRAGG ::")
control = control[:bragg_start].rstrip() + "\n"
(run_dir / "SRTIO3_5K_RMC7.dat").write_text(control)

docker(package_root, run_dir, "rmcprofile", "SRTIO3_5K_RMC7")
required = (
    "SRTIO3_5K.rmc7",
    "SRTIO3_5K.rmc7_grpartials.csv",
    "SRTIO3_5K.rmc7_n_sqpartials.csv",
    "srtio3_wk_5k_rmc.fq_1.csv",
    "srtio3_wk_5k_rmc.gr_1.csv",
)
missing = [name for name in required if not (run_dir / name).exists()]
if missing:
    raise SystemExit(f"RMCProfile did not produce: {', '.join(missing)}")
print(f"native RMCProfile7 outputs are in {run_dir}")
