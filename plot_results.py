#!/usr/bin/env python3
"""Plot rsmith computed S(Q) and g(r) against experimental data.

Usage: python3 plot_results.py [config.toml] [output_dir]

Defaults to config.toml in the current directory and output in current directory.
Automatically detects which datasets (xray_sq, neutron_sq, xray_gr) are present.
"""

import os
import re
import sys
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
    return n_atoms, vals[0]  # assume cubic


def compute_weighted_gr(Q, sx, rho0, r_out, lorch=True):
    """Fourier transform total S(Q) -> weighted g(r)."""
    h = (sx - 1.0).copy()
    if lorch:
        Q_max = Q[-1]
        W = np.where(Q > 0, np.sin(np.pi * Q / Q_max) / (np.pi * Q / Q_max), 1.0)
        h *= W
    sinQr = np.sin(np.outer(r_out, Q))
    integrand = Q[np.newaxis, :] * h[np.newaxis, :] * sinQr
    integral = _trapz(integrand, Q, axis=1)
    return np.where(r_out > 0.01, 1.0 + integral / (2.0 * np.pi**2 * rho0 * r_out), 1.0)


def find_sq_file(output_dir, kind, data_type):
    """Find S(Q) file, trying new name then old name."""
    for name in [f"{kind}_{data_type}_sq.dat", f"{kind}_sq.dat"]:
        path = os.path.join(output_dir, name)
        if os.path.isfile(path):
            return path
    return None


# --- Parse arguments ---
config_path = sys.argv[1] if len(sys.argv) > 1 else "config.toml"
output_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(config_path) or "."
config_dir = os.path.dirname(config_path) or "."

with open(config_path, "rb") as f:
    config = tomllib.load(f)

# Number density
xyz_refined = os.path.join(output_dir, "refined.xyz")
xyz_path = xyz_refined if os.path.exists(xyz_refined) else os.path.join(
    config_dir, config["system"]["structure"]
)
N_TOTAL, box_L = read_xyz_header(xyz_path)
rho0 = N_TOTAL / box_L**3

# --- Discover datasets ---
datasets = []  # (label, data_type, exp_path, convention)
for data_type, section_key in [("xray", "xray_sq"), ("neutron", "neutron_sq")]:
    section = config.get("data", {}).get(section_key)
    if section is None:
        continue
    exp_path = os.path.join(config_dir, section["file"])
    convention = section.get("convention", "sq")
    datasets.append((section_key, data_type, exp_path, convention))

has_gr_exp = config.get("data", {}).get("xray_gr") is not None

# --- Load partial g(r) ---
gr_start_path = os.path.join(output_dir, "start_gr.dat")
gr_refined_path = os.path.join(output_dir, "refined_gr.dat")
has_gr_start = os.path.isfile(gr_start_path)
has_gr_refined = os.path.isfile(gr_refined_path)

if has_gr_start:
    gr_start = np.loadtxt(gr_start_path, comments="#")
    with open(gr_start_path) as f:
        header = f.readline().strip().lstrip("# ").split()
    pair_labels = [col.replace("g_", "").replace("_", "") for col in header[1:]]
    pair_labels = [re.sub(r"([A-Z][a-z]?)([A-Z])", r"\1-\2", lbl) for lbl in pair_labels]
if has_gr_refined:
    gr_refined = np.loadtxt(gr_refined_path, comments="#")

# --- Determine layout ---
n_sq_panels = len(datasets)
n_gr_panels = 0
if has_gr_exp:
    n_gr_panels += 1  # weighted g(r) vs experiment
if has_gr_start or has_gr_refined:
    n_gr_panels += 1  # partial g(r)
n_panels = n_sq_panels + n_gr_panels

if n_panels == 0:
    print("No data to plot.")
    sys.exit(0)

fig, axes = plt.subplots(n_panels, 1, figsize=(8, 4 * n_panels))
if n_panels == 1:
    axes = [axes]
panel = 0

# --- S(Q) panels ---
for section_key, data_type, exp_path, convention in datasets:
    ax = axes[panel]
    panel += 1

    sq_exp = np.loadtxt(exp_path)
    sq_start_path = find_sq_file(output_dir, "start", data_type)
    sq_refined_path = find_sq_file(output_dir, "refined", data_type)

    q_max_exp = sq_exp[:, 0].max()

    if sq_start_path:
        sq_start = np.loadtxt(sq_start_path, comments="#")
        mask = sq_start[:, 0] <= q_max_exp
        ax.plot(sq_start[mask, 0], sq_start[mask, 1], "b-", lw=1.2, label="Starting config")
    if sq_refined_path:
        sq_ref = np.loadtxt(sq_refined_path, comments="#")
        mask = sq_ref[:, 0] <= q_max_exp
        ax.plot(sq_ref[mask, 0], sq_ref[mask, 1], "g-", lw=1.5, label="RMC refined")
    ax.plot(sq_exp[:, 0], sq_exp[:, 1], "r-", lw=1.5, alpha=0.8, label="Experiment")

    # Convention-dependent baseline
    if convention == "sq":
        ax.axhline(1, color="gray", ls="--", lw=0.5)
        ylabel = "S(Q)"
    elif convention == "iq":
        ax.axhline(0, color="gray", ls="--", lw=0.5)
        ylabel = "i(Q) = S(Q) - 1"
    elif convention == "fq":
        ax.axhline(0, color="gray", ls="--", lw=0.5)
        ylabel = "F(Q) = Q(S(Q) - 1)"
    else:
        ylabel = "S(Q)"

    label_nice = "X-ray" if data_type == "xray" else "Neutron"
    ax.set(xlabel="Q (1/\u00c5)", ylabel=ylabel, xlim=(0.3, q_max_exp))
    ax.set_title(f"{label_nice} structure factor")
    ax.legend()

# --- Weighted g(r) panel (if experimental g(r) available) ---
if has_gr_exp:
    ax = axes[panel]
    panel += 1

    gr_exp_data = np.loadtxt(os.path.join(config_dir, config["data"]["xray_gr"]["file"]))
    r_gr = np.linspace(0.0, 10.0, 1000)

    # Use experimental S(Q) Q range for FT to avoid truncation artifacts
    sq_start_path = find_sq_file(output_dir, "start", "xray")
    sq_refined_path = find_sq_file(output_dir, "refined", "xray")
    # Find experimental S(Q) Q_max to truncate computed S(Q) before FT
    xray_sq_cfg = config.get("data", {}).get("xray_sq")
    if xray_sq_cfg:
        sq_exp_for_qmax = np.loadtxt(os.path.join(config_dir, xray_sq_cfg["file"]))
        q_cut = sq_exp_for_qmax[:, 0].max()
    else:
        q_cut = None
    if sq_start_path:
        sq_s = np.loadtxt(sq_start_path, comments="#")
        if q_cut is not None:
            m = sq_s[:, 0] <= q_cut
            sq_s = sq_s[m]
        gx_start = compute_weighted_gr(sq_s[:, 0], sq_s[:, 1], rho0, r_gr)
        ax.plot(r_gr, gx_start, "b-", lw=1.2, label="Starting config")
    if sq_refined_path:
        sq_r = np.loadtxt(sq_refined_path, comments="#")
        if q_cut is not None:
            m = sq_r[:, 0] <= q_cut
            sq_r = sq_r[m]
        gx_ref = compute_weighted_gr(sq_r[:, 0], sq_r[:, 1], rho0, r_gr)
        ax.plot(r_gr, gx_ref, "g-", lw=1.5, label="RMC refined")
    ax.plot(gr_exp_data[:, 0], gr_exp_data[:, 1], "r-", lw=1.5, alpha=0.8, label="Experiment")

    ax.axhline(1, color="gray", ls="--", lw=0.5)
    ax.set(xlabel="r (\u00c5)", ylabel="g(r)", xlim=(0, 10))
    ax.set_title("X-ray weighted pair distribution function")
    ax.legend()

# --- Partial g(r) panel ---
if has_gr_start or has_gr_refined:
    ax = axes[panel]
    panel += 1

    gr_data = gr_refined if has_gr_refined else gr_start
    r = gr_data[:, 0]
    for i, label in enumerate(pair_labels):
        ax.plot(r, gr_data[:, 1 + i], lw=1.2, label=label)
    ax.axhline(1, color="gray", ls="--", lw=0.5)
    ax.set(xlabel="r (\u00c5)", ylabel="g(r)", xlim=(0, 10))
    ax.set_title("Partial pair distribution functions")
    ax.legend(ncol=2)

plt.tight_layout()
outfile = os.path.join(output_dir, "results_comparison.png")
plt.savefig(outfile, dpi=150, bbox_inches="tight")
print(f"Saved {outfile}")
plt.show()
