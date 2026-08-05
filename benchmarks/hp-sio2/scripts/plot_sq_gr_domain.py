#!/usr/bin/env python3
"""Plot reciprocal, total-real, and held-out partial curves for the domain test."""

from __future__ import annotations

import argparse
import math
import tomllib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pressure_common import reciprocal_curves, structure_metrics, total_gr_from_iq


plt.rcParams["svg.hashsalt"] = "rsmith-hp-sio2"

parser = argparse.ArgumentParser()
parser.add_argument("--output-dir", type=Path)
args = parser.parse_args()
case_root = Path(__file__).resolve().parents[1]
root = case_root / "results/sq-gr-domain"
protocol = tomllib.loads((case_root / "expected/sq-gr-domain.toml").read_text())
fixture = protocol["fixture"]
output = args.output_dir or case_root / "figures"
output.mkdir(parents=True, exist_ok=True)


def read_curve(path):
    return np.asarray(
        [
            [float(value) for value in line.split()[:2]]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
    )


q_target = read_curve(root / "inputs/target-xray-iq.dat")
gr_target = read_curve(root / "inputs/target-xray-gr.dat")
q_values, r_values = q_target[:, 0], gr_target[:, 0]
paths = [
    ("Hidden target", root / "inputs/hidden-target.data", "#111111", "-", 2.0),
    ("Mapped start", root / "inputs/cross-start.data", "#999999", "--", 1.3),
    ("$S(Q)$ only", root / "runs/sq-only/refined.xyz", "#d55e00", "-", 1.3),
    ("local $g(r)$ only", root / "runs/gr-only/refined.xyz", "#009e73", "-", 1.3),
    ("$S(Q)+g(r)$", root / "runs/sq-plus-gr/refined.xyz", "#0072b2", "-", 1.6),
]
curves = {}
for label, path, color, linestyle, linewidth in paths:
    metrics = structure_metrics(path, fixture)
    _, _, iq = reciprocal_curves(metrics, q_values, fixture["rdf_bin_width_a"])
    density = fixture["atoms"] / math.prod(metrics["box_a"])
    total_gr = total_gr_from_iq(
        q_values,
        iq,
        r_values,
        density,
        fixture["gr_qmax_a_inverse"],
        fixture["gr_lorch"],
    )
    curves[label] = (metrics, iq, total_gr, color, linestyle, linewidth)

fig, axes = plt.subplots(2, 3, figsize=(11.0, 6.3), constrained_layout=True)
for label, (_, iq, _, color, linestyle, linewidth) in curves.items():
    axes[0, 0].plot(q_values, iq - q_target[:, 1], color=color, linestyle=linestyle, linewidth=linewidth, label=label)
axes[0, 0].axhline(0.0, color="#222222", linewidth=0.6)
axes[0, 0].set(xlabel=r"$Q$ ($\mathrm{Å}^{-1}$)", ylabel=r"$i_X(Q)-i_{X,target}(Q)$", xlim=(0.5, 17.0))

for label, (_, _, total_gr, color, linestyle, linewidth) in curves.items():
    axes[0, 1].plot(r_values, total_gr, color=color, linestyle=linestyle, linewidth=linewidth, label=label)
axes[0, 1].set(xlabel=r"$r$ (Å)", ylabel=r"total $g_X(r)$", xlim=(1.2, 6.0))

metric_names = ["X-ray $i(Q)$", "total $g_X(r)$", "local partial $g(r)$", "full partial $g(r)$"]
score = tomllib.loads((case_root / "expected/sq-gr-domain-observed.toml").read_text())
method_records = [
    ("$S(Q)$", score["methods"]["sq_only"], "#d55e00"),
    ("$g(r)$", score["methods"]["gr_only"], "#009e73"),
    ("joint", score["methods"]["sq_plus_gr"], "#0072b2"),
]
metric_keys = ["xray_iq_rms", "xray_gr_local_rms", "local_partial_rdf_rms", "full_partial_rdf_rms"]
x = np.arange(len(metric_names))
width = 0.24
for index, (label, record, color) in enumerate(method_records):
    axes[0, 2].bar(x + (index - 1) * width, [record[key] for key in metric_keys], width, label=label, color=color)
axes[0, 2].set_xticks(x, metric_names, rotation=25, ha="right")
axes[0, 2].set_ylabel("RMS")
axes[0, 2].legend(frameon=False, fontsize=8)

rdf_radius = (np.arange(len(next(iter(curves["Hidden target"][0]["rdf"].values())))) + 0.5) * fixture["rdf_bin_width_a"]
for axis, pair in zip(axes[1], ("Si-Si", "Si-O", "O-O")):
    for label, (metrics, _, _, color, linestyle, linewidth) in curves.items():
        axis.plot(rdf_radius, metrics["rdf"][pair], color=color, linestyle=linestyle, linewidth=linewidth, label=label)
    axis.set(xlabel=r"$r$ (Å)", ylabel=fr"$g_{{\mathrm{{{pair.replace('-', '')}}}}}(r)$", xlim=(1.0, 6.0))

for panel, axis in enumerate(axes.flat):
    axis.grid(alpha=0.22, linewidth=0.6)
    axis.text(0.02, 0.96, chr(ord("a") + panel), transform=axis.transAxes, va="top", fontweight="bold")
axes[0, 0].legend(frameon=False, fontsize=7.5, ncol=2)
for suffix in ("svg", "png"):
    path = output / f"sq-gr-domain.{suffix}"
    fig.savefig(
        path,
        dpi=240,
        bbox_inches="tight",
        metadata={"Date": None},
    )
    if suffix == "svg":
        path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
plt.close(fig)
