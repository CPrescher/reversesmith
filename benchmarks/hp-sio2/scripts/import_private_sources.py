#!/usr/bin/env python3
"""Import pinned private 0/10 GPa SiO2 sources through authenticated gh."""

from __future__ import annotations

import hashlib
import json
import subprocess
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
protocol_path = case_root / "expected/ten-gpa-pilot.toml"
protocol = tomllib.loads(protocol_path.read_text())
source = protocol["source"]
output = case_root / "reference/local/ten-gpa"
output.mkdir(parents=True, exist_ok=True)

files = {
    "compression-input.in": (
        source["compression_input"], source["compression_input_sha256"]
    ),
    "ambient-0gpa.data": (
        source["ambient_structure"], source["ambient_structure_sha256"]
    ),
    "target-10gpa.data": (
        source["target_structure"], source["target_structure_sha256"]
    ),
}


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


manifest_files = {}
for local_name, (repository_path, expected_hash) in files.items():
    endpoint = (
        f"repos/{source['repository']}/contents/{repository_path}"
        f"?ref={source['commit']}"
    )
    completed = subprocess.run(
        [
            "gh",
            "api",
            "-H",
            "Accept: application/vnd.github.raw+json",
            endpoint,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr.decode(errors="replace"))
    digest = sha256(completed.stdout)
    if digest != expected_hash:
        raise SystemExit(
            f"{repository_path}: downloaded hash {digest} != {expected_hash}"
        )
    target = output / local_name
    target.write_bytes(completed.stdout)
    manifest_files[local_name] = {
        "repository_path": repository_path,
        "sha256": digest,
        "bytes": len(completed.stdout),
    }

manifest = {
    "status": "pinned_private_hp_sio2_sources_imported",
    "repository": source["repository"],
    "commit": source["commit"],
    "protocol_sha256": hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
    "files": manifest_files,
}
(output / "IMPORT.json").write_text(
    json.dumps(manifest, indent=2, sort_keys=True) + "\n"
)
print(output / "IMPORT.json")
