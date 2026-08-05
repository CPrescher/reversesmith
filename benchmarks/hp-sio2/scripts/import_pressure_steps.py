#!/usr/bin/env python3
"""Import pinned 0, 5, and 10 GPa structures through authenticated gh."""

from __future__ import annotations

import hashlib
import json
import subprocess
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
protocol_path = case_root / "expected/pressure-steps.toml"
protocol = tomllib.loads(protocol_path.read_text())
output = case_root / "reference/local/pressure-steps"
output.mkdir(parents=True, exist_ok=True)
records = {}
for structure in protocol["structures"]:
    pressure = int(structure["pressure_gpa"])
    endpoint = (
        f"repos/{protocol['source_repository']}/contents/{structure['path']}"
        f"?ref={protocol['source_commit']}"
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
    digest = hashlib.sha256(completed.stdout).hexdigest()
    if digest != structure["sha256"]:
        raise SystemExit(f"{structure['path']}: {digest} != {structure['sha256']}")
    path = output / f"{pressure:02d}gpa.data"
    path.write_bytes(completed.stdout)
    records[str(pressure)] = {
        "bytes": len(completed.stdout),
        "path": structure["path"],
        "sha256": digest,
    }
manifest = {
    "commit": protocol["source_commit"],
    "files": records,
    "protocol_sha256": hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
    "repository": protocol["source_repository"],
    "status": "hp_sio2_pressure_steps_imported",
}
(output / "IMPORT.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
print(output / "IMPORT.json")
