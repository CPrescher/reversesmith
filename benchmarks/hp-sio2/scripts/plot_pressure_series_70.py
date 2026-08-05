#!/usr/bin/env python3
"""Plot the frozen four-program pressure-series metrics."""

from __future__ import annotations

import argparse
import tomllib
from pathlib import Path

import matplotlib.pyplot as plt


plt.rcParams["svg.hashsalt"] = "rsmith-hp-sio2"

parser = argparse.ArgumentParser()
parser.add_argument("--output-dir", type=Path)
args = parser.parse_args()
case_root = Path(__file__).resolve().parents[1]
observed = tomllib.loads(
    (case_root / "expected/pressure-series-70-observed.toml").read_text()
)
output = args.output_dir or case_root / "figures"
output.mkdir(parents=True, exist_ok=True)
pressure = observed["target_pressure_gpa"]
series = [
    ("start", observed["start"], "Mapped start", "#222222", "--", "o"),
    ("rsmith_rmc", observed["methods"]["rsmith_rmc"], "rsmith RMC", "#7a7a7a", "-", "s"),
    ("rsmith_pace_w30", observed["methods"]["rsmith_pace_w30"], "rsmith PACE", "#0072b2", "-", "o"),
    ("native_epsr26", observed["methods"]["native_epsr26"], "EPSR26", "#d55e00", "-", "^"),
    ("native_rmcprofile", observed["methods"]["native_rmcprofile"], "RMCProfile", "#009e73", "-", "D"),
]
metrics = [
    ("mean_partial_rdf_rms", "Held-out partial $g(r)$ RMS", True),
    ("joint_iq_rms", "Common neutron/X-ray $i(Q)$ RMS", True),
    ("si_coordination_tv", "Si-coordination total variation", False),
    ("mean_lower_tail_error_a", "Mean lower-tail error (Å)", True),
]
fig, axes = plt.subplots(2, 2, figsize=(8.0, 6.2), constrained_layout=True)
for panel, (axis, (key, ylabel, log_scale)) in enumerate(zip(axes.flat, metrics)):
    for _, values, label, color, linestyle, marker in series:
        axis.plot(
            pressure,
            values[key],
            label=label,
            color=color,
            linestyle=linestyle,
            marker=marker,
            linewidth=1.6,
            markersize=4.5,
        )
    if log_scale:
        axis.set_yscale("log")
    axis.set_xlabel("Target pressure (GPa)")
    axis.set_ylabel(ylabel)
    axis.grid(alpha=0.25, linewidth=0.6)
    axis.text(0.02, 0.96, chr(ord("a") + panel), transform=axis.transAxes, va="top", fontweight="bold")
axes[0, 0].legend(frameon=False, fontsize=8, ncol=2)
for suffix in ("svg", "png"):
    path = output / f"pressure-series-70.{suffix}"
    fig.savefig(
        path,
        dpi=240,
        bbox_inches="tight",
        metadata={"Date": None},
    )
    if suffix == "svg":
        path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
plt.close(fig)
