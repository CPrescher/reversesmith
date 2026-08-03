#!/usr/bin/env python3
"""Fetch or import the public Erhard silica GAP model with checksum gates."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import urllib.request
import zipfile
from pathlib import Path


URL = "https://zenodo.org/records/6353684/files/sio2_potential_data.zip?download=1"
EXPECTED_MD5 = "6a16de69b5e17fd18160d9b55972a3e1"
ARCHIVE_PREFIX = "sio2_potential_data/potential/"


def digest(path: Path, algorithm: str) -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


parser = argparse.ArgumentParser()
parser.add_argument("--archive", type=Path, help="use an existing downloaded ZIP")
parser.add_argument("--force", action="store_true")
args = parser.parse_args()

case_root = Path(__file__).resolve().parents[1]
output = case_root / "reference/local/public-gap"
archive = args.archive or output / "sio2_potential_data.zip"
if output.exists() and not args.force:
    raise SystemExit(f"output exists: {output} (pass --force to replace it)")
if output.exists():
    shutil.rmtree(output)
output.mkdir(parents=True)
if args.archive is None:
    print(f"Downloading {URL}")
    urllib.request.urlretrieve(URL, archive)
if not archive.is_file():
    raise SystemExit(f"archive not found: {archive}")
actual_md5 = digest(archive, "md5")
if actual_md5 != EXPECTED_MD5:
    raise SystemExit(f"archive MD5 mismatch: {actual_md5} != {EXPECTED_MD5}")

with zipfile.ZipFile(archive) as bundle:
    members = [
        name
        for name in bundle.namelist()
        if name.startswith(ARCHIVE_PREFIX) and not name.endswith("/")
    ]
    if not any(name.endswith("silica_gap.xml") for name in members):
        raise SystemExit("silica_gap.xml missing from verified archive")
    for member in members:
        destination = output / Path(member).name
        with bundle.open(member) as source, destination.open("wb") as target:
            shutil.copyfileobj(source, target)

files = {
    path.name: {"bytes": path.stat().st_size, "sha256": digest(path, "sha256")}
    for path in sorted(output.glob("silica_gap.xml*"))
}
record = {
    "source": URL,
    "doi": "10.5281/zenodo.6353684",
    "archive_md5": actual_md5,
    "archive_sha256": digest(archive, "sha256"),
    "files": files,
}
(output / "IMPORT.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
print(json.dumps(record, indent=2, sort_keys=True))
print(f"GAP model: {output / 'silica_gap.xml'}")
