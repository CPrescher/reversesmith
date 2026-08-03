#!/usr/bin/env python3
"""Fetch or import the public Erhard 2024 Si-O ACE potential with checksum gates."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import urllib.request
import zipfile
from pathlib import Path


URL = "https://zenodo.org/api/records/10419194/files/potential.zip/content"
DOI = "10.5281/zenodo.10419194"
EXPECTED_MD5 = "c8eb4cd111af60d10c15bc3f1de9adbb"
EXPECTED_ARCHIVE_SHA256 = (
    "28b29becd2c3185c6a44e872f304af7689b30b22842f21ec91e52f3641dd72cb"
)
MODEL_MEMBER = "potential/SiOx_potential.yace"
EXPECTED_MODEL_SHA256 = (
    "c8f00d8f0cbc131b0298b79260ba8098975624363c9d178223e51e48f025e97a"
)


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
output = case_root / "reference/local/public-ace2024"
archive = args.archive or output / "potential.zip"
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
actual_sha256 = digest(archive, "sha256")
if actual_md5 != EXPECTED_MD5 or actual_sha256 != EXPECTED_ARCHIVE_SHA256:
    raise SystemExit(
        f"archive checksum mismatch: md5={actual_md5}, sha256={actual_sha256}"
    )

model = output / "SiOx_potential.yace"
with zipfile.ZipFile(archive) as bundle, bundle.open(MODEL_MEMBER) as source:
    with model.open("wb") as target:
        shutil.copyfileobj(source, target)
model_sha256 = digest(model, "sha256")
if model_sha256 != EXPECTED_MODEL_SHA256:
    raise SystemExit(f"model SHA-256 mismatch: {model_sha256}")

record = {
    "source": URL,
    "doi": DOI,
    "article_doi": "10.1038/s41467-024-45840-9",
    "archive_md5": actual_md5,
    "archive_sha256": actual_sha256,
    "model": {
        "name": model.name,
        "bytes": model.stat().st_size,
        "sha256": model_sha256,
    },
}
(output / "IMPORT.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
print(json.dumps(record, indent=2, sort_keys=True))
print(f"ACE model: {model}")
