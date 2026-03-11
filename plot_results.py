#!/usr/bin/env python3
"""Plot rsmith computed S(Q) and g(r) against experimental data."""

import os
import re
import tomllib
import numpy as np
import matplotlib.pyplot as plt

_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")


def read_xyz_header(path):
    """Read atom count and box length from an extended XYZ file."""
    with open(path) as f:
        n_atoms = int(f.readline().strip())
        lattice_line = f.readline()
    m = re.search(r'Lattice="([^"]+)"', lattice_line)
    vals = [float(x) for x in m.group(1).split()]
    box_length = vals[0]  # assume cubic
    return n_atoms, box_length


def compute_xray_gr(Q, sx, rho0, r_out, lorch=True):
    """Fourier transform total X-ray S(Q) -> g(r)."""
    h = (sx - 1.0).copy()
    if lorch:
        Q_max = Q[-1]
        W = np.where(Q > 0, np.sin(np.pi * Q / Q_max) / (np.pi * Q / Q_max), 1.0)
        h *= W
    sinQr = np.sin(np.outer(r_out, Q))
    integrand = Q[np.newaxis, :] * h[np.newaxis, :] * sinQr
    integral = _trapz(integrand, Q, axis=1)
    return np.where(r_out > 0.01, 1.0 + integral / (2.0 * np.pi**2 * rho0 * r_out), 1.0)


# --- Load data ---
with open("config.toml", "rb") as f:
    config = tomllib.load(f)

sq_start = np.loadtxt("start_sq.dat", comments="#")
gr_start = np.loadtxt("start_gr.dat", comments="#")
sq_exp = np.loadtxt(config["data"]["xray_sq"]["file"])
gr_exp = np.loadtxt(config["data"]["xray_gr"]["file"])

has_refined = os.path.exists("refined_sq.dat")
if has_refined:
    sq_refined = np.loadtxt("refined_sq.dat", comments="#")
    gr_refined = np.loadtxt("refined_gr.dat", comments="#")

# Number density from refined (or starting) structure
xyz_path = "refined.xyz" if os.path.exists("refined.xyz") else config["system"]["structure"]
N_TOTAL, box_L = read_xyz_header(xyz_path)
rho0 = N_TOTAL / box_L**3

Q_calc = sq_start[:, 0]

# Parse pair labels from header (e.g. "g_SiO" -> "Si-O")
with open("start_gr.dat") as f:
    header = f.readline().strip().lstrip("# ").split()
pair_labels = [col.replace("g_", "").replace("_", "") for col in header[1:]]
pair_labels = [re.sub(r"([A-Z][a-z]?)([A-Z])", r"\1-\2", lbl) for lbl in pair_labels]

# Fourier transform total S(Q) -> X-ray weighted g(r)
r_gr = np.linspace(0.0, 10.0, 1000)
gx_start = compute_xray_gr(Q_calc, sq_start[:, 1], rho0, r_gr)
if has_refined:
    gx_refined = compute_xray_gr(sq_refined[:, 0], sq_refined[:, 1], rho0, r_gr)

# --- Plot ---
fig, axes = plt.subplots(3, 1, figsize=(8, 12))

# 1. Total X-ray S(Q): starting + refined vs experiment
ax = axes[0]
ax.plot(Q_calc, sq_start[:, 1], "b-", lw=1.5, label="MD starting config")
if has_refined:
    ax.plot(sq_refined[:, 0], sq_refined[:, 1], "g-", lw=1.5, label="RMC refined")
ax.plot(sq_exp[:, 0], sq_exp[:, 1], "r-", lw=1.5, alpha=0.8, label="Experiment")
ax.axhline(1, color="gray", ls="--", lw=0.5)
ax.set(
    xlabel="Q (1/Å)",
    ylabel="S(Q)",
    xlim=(0.3, 18),
    title="X-ray total structure factor — CaSiO₃ glass",
)
ax.legend()

# 2. X-ray weighted g(r): computed vs experiment
ax = axes[1]
ax.plot(r_gr, gx_start, "b-", lw=1.5, label="MD starting config")
if has_refined:
    ax.plot(r_gr, gx_refined, "g-", lw=1.5, label="RMC refined")
ax.plot(gr_exp[:, 0], gr_exp[:, 1], "r-", lw=1.5, alpha=0.8, label="Experiment")
ax.axhline(1, color="gray", ls="--", lw=0.5)
ax.set(
    xlabel="r (Å)",
    ylabel="g(r)",
    xlim=(0, 10),
    title="X-ray weighted pair distribution function — CaSiO₃ glass",
)
ax.legend()

# 3. Partial g(r) (refined if available, else starting)
ax = axes[2]
gr_partials = gr_refined if has_refined else gr_start
r = gr_partials[:, 0]
for i, label in enumerate(pair_labels):
    ax.plot(r, gr_partials[:, 1 + i], lw=1.2, label=label)
ax.axhline(1, color="gray", ls="--", lw=0.5)
ax.set(
    xlabel="r (Å)",
    ylabel="g(r)",
    xlim=(0, 10),
    title="Partial pair distribution functions — rsmith",
)
ax.legend(ncol=2)

plt.tight_layout()
plt.savefig("results_comparison.png", dpi=150, bbox_inches="tight")
print("Saved results_comparison.png")
plt.show()
