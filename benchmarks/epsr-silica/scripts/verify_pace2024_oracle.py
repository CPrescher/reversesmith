#!/usr/bin/env python3
"""Verify the pinned Erhard-2024 ACE provenance and stored oracle agreement."""

from __future__ import annotations

import hashlib
import math
import tomllib
from pathlib import Path


case_root = Path(__file__).resolve().parents[1]
expected = tomllib.loads((case_root / "expected/pace2024-oracle.toml").read_text())
model = case_root / "reference/local/public-ace2024/SiOx_potential.yace"
failures: list[str] = []

if model.is_file():
    actual_hash = hashlib.sha256(model.read_bytes()).hexdigest()
    if actual_hash != expected["source"]["model_sha256"]:
        failures.append(f"model SHA-256 mismatch: {actual_hash}")
    if model.stat().st_size != expected["source"]["model_bytes"]:
        failures.append(f"model byte count mismatch: {model.stat().st_size}")
else:
    print(f"note: local model absent; provenance-only checks for {model}")

observed = expected["observed"]
guards = expected["guards"]
for label in ("initial", "moved"):
    error = abs(
        observed[f"oracle_{label}_energy_ev"]
        - observed[f"rsmith_{label}_energy_ev_printed"]
    )
    if not math.isfinite(error) or error > guards["total_energy_abs_error_ev_max"]:
        failures.append(f"{label} total-energy error {error} exceeds guard")
delta_error = abs(
    observed["oracle_move_delta_ev"]
    - observed["rsmith_move_delta_ev_from_printed_totals"]
)
if not math.isfinite(delta_error) or delta_error > guards["move_delta_abs_error_ev_max"]:
    failures.append(f"move-delta error {delta_error} exceeds guard")

if failures:
    raise SystemExit("PACE 2024 oracle verification failed:\n- " + "\n- ".join(failures))
print("passed pinned Erhard-2024 ACE provenance and energy-oracle checks")
