#!/usr/bin/env python3
"""Plot reversesmith computed S(Q) and g(r) against experimental data."""

import os
import numpy as np
import matplotlib.pyplot as plt

_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")

N_TOTAL = 1000


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
    return np.where(
        r_out > 0.01, 1.0 + integral / (2.0 * np.pi**2 * rho0 * r_out), 1.0
    )


# --- Load data ---
sq_start = np.loadtxt("computed_sq.dat", comments="#")
gr_calc = np.loadtxt("computed_gr.dat", comments="#")
sq_exp = np.loadtxt("../CaSiO3_ambient_sample.sq")
gr_exp = np.loadtxt("../CaSiO3_ambient_sample.gr")

has_refined = os.path.exists("refined_sq.dat")
if has_refined:
    sq_refined = np.loadtxt("refined_sq.dat", comments="#")

# Number density
rho0 = N_TOTAL / 23.858889664691711**3

Q_calc = sq_start[:, 0]
pair_labels = ["Ca-Ca", "Ca-Si", "Ca-O", "Si-Si", "Si-O", "O-O"]

# Fourier transform total S(Q) -> X-ray weighted g(r)
r_gr = np.linspace(0.0, 10.0, 1000)
gx_start = compute_xray_gr(Q_calc, sq_start[:, 1], rho0, r_gr)
if has_refined:
    gx_refined = compute_xray_gr(sq_refined[:, 0], sq_refined[:, 1], rho0, r_gr)

# --- Plot ---
fig, axes = plt.subplots(3, 1, figsize=(8, 12))

# 1. Partial g(r)
ax = axes[0]
r = gr_calc[:, 0]
for i, label in enumerate(pair_labels):
    ax.plot(r, gr_calc[:, 1 + i], lw=1.2, label=label)
ax.axhline(1, color="gray", ls="--", lw=0.5)
ax.set(
    xlabel="r (Å)", ylabel="g(r)", xlim=(0, 10),
    title="Partial pair distribution functions — reversesmith",
)
ax.legend(ncol=2)

# 2. Total X-ray S(Q): starting + refined vs experiment
ax = axes[1]
ax.plot(Q_calc, sq_start[:, 1], "b-", lw=1.5, label="MD starting config")
if has_refined:
    ax.plot(sq_refined[:, 0], sq_refined[:, 1], "g-", lw=1.5, label="RMC refined")
ax.plot(sq_exp[:, 0], sq_exp[:, 1], "r-", lw=1.5, alpha=0.8, label="Experiment")
ax.axhline(1, color="gray", ls="--", lw=0.5)
ax.set(
    xlabel="Q (1/Å)", ylabel="S(Q)", xlim=(0.3, 18),
    title="X-ray total structure factor — CaSiO₃ glass",
)
ax.legend()

# 3. X-ray weighted g(r): computed vs experiment
ax = axes[2]
ax.plot(r_gr, gx_start, "b-", lw=1.5, label="MD starting config")
if has_refined:
    ax.plot(r_gr, gx_refined, "g-", lw=1.5, label="RMC refined")
ax.plot(gr_exp[:, 0], gr_exp[:, 1], "r-", lw=1.5, alpha=0.8, label="Experiment")
ax.axhline(1, color="gray", ls="--", lw=0.5)
ax.set(
    xlabel="r (Å)", ylabel="g(r)", xlim=(0, 10),
    title="X-ray weighted pair distribution function — CaSiO₃ glass",
)
ax.legend()

plt.tight_layout()
plt.savefig("results_comparison.png", dpi=150, bbox_inches="tight")
print("Saved results_comparison.png")
plt.show()
