#!/usr/bin/env python3
"""Compare results across ensemble runs.

Usage: python3 compare_ensemble.py [config.toml] [run01 run02 ...]

If no run directories are given, auto-discovers run01..run99 in the config
file's directory. Reads each run's rsmith.log and analysis files, then prints
a summary table and optionally overlays S(Q) fits.
"""

import os
import re
import sys
import numpy as np

try:
    import matplotlib.pyplot as plt  # type: ignore[reportPossiblyUnbound]
    HAS_PLT = True
except ImportError:
    plt = None  # type: ignore[assignment]
    HAS_PLT = False


def parse_log(log_path):
    """Extract key metrics from rsmith.log."""
    info = {}
    with open(log_path) as f:
        for line in f:
            # RNG seed
            m = re.search(r"RNG seed \(\w+\): (\d+)", line)
            if m:
                info["seed"] = int(m.group(1))

            # Final chi2 line:
            # "Final chi2 = 1.234567 (best at move 450000), accepted 234567/500000 (46.9%)"
            m = re.search(
                r"Final chi2 = ([\d.]+) \(best at move (\d+)\), accepted (\d+)/(\d+) \(([\d.]+)%\)",
                line,
            )
            if m:
                info["chi2"] = float(m.group(1))
                info["best_move"] = int(m.group(2))
                info["accepted"] = int(m.group(3))
                info["total_moves"] = int(m.group(4))
                info["acceptance"] = float(m.group(5))

            # Density (from rescaling line or initial output)
            m = re.search(r"density = ([\d.]+) g/cm", line)
            if m:
                info["density"] = float(m.group(1))
            # Also pick up from initial volume line if no rescaling
            m = re.search(r"rho0 = ([\d.]+) atoms/A\^3", line)
            if m:
                info["rho0"] = float(m.group(1))

            # Last status line (fallback for chi2 components)
            m = re.match(r"Move \d+/\d+: (.+)", line)
            if m:
                body = m.group(1)
                # Extract sq/gr components if present
                m2 = re.search(r"sq: ([\d.]+)", body)
                if m2:
                    info["sq_chi2"] = float(m2.group(1))
                m2 = re.search(r"gr: ([\d.]+)", body)
                if m2:
                    info["gr_chi2"] = float(m2.group(1))

    return info


def parse_cn_file(cn_path):
    """Extract mean coordination numbers from analysis CN file.

    Reads summary comment lines like:
    # Ca-O: mean = 7.123, std = 1.234, min = 4, max = 11, cutoff = 3.00 A
    """
    cns = {}
    with open(cn_path) as f:
        for line in f:
            m = re.match(r"# (\S+): mean = ([\d.]+)", line)
            if m:
                cns[m.group(1)] = float(m.group(2))
    return cns


# --- Parse arguments ---
config_path = sys.argv[1] if len(sys.argv) > 1 else "config.toml"
config_dir = os.path.dirname(config_path) or "."

if len(sys.argv) > 2:
    run_dirs = sys.argv[2:]
else:
    # Auto-discover run directories: any subdirectory containing rsmith.log
    run_dirs = sorted(
        d
        for d in os.listdir(config_dir)
        if os.path.isdir(os.path.join(config_dir, d))
        and os.path.isfile(os.path.join(config_dir, d, "rsmith.log"))
    )
    run_dirs = [os.path.join(config_dir, d) for d in run_dirs]

if not run_dirs:
    print("No run directories found. Usage: compare_ensemble.py [config.toml] [run01 run02 ...]")
    sys.exit(1)

# --- Collect metrics ---
results = []
for run_dir in run_dirs:
    log_path = os.path.join(run_dir, "rsmith.log")
    if not os.path.isfile(log_path):
        print(f"  Skipping {run_dir}: no rsmith.log")
        continue

    info = parse_log(log_path)
    info["dir"] = os.path.basename(run_dir)

    # Try to read CN data
    for cn_name in ["analysis_refined_cn.dat", "analysis_analysis_cn.dat"]:
        cn_path = os.path.join(run_dir, cn_name)
        if os.path.isfile(cn_path):
            info["cn"] = parse_cn_file(cn_path)
            break

    results.append(info)

if not results:
    print("No valid runs found.")
    sys.exit(1)

# Sort by chi2
results.sort(key=lambda r: r.get("chi2", float("inf")))

# --- Print summary table ---
print(f"\nEnsemble summary ({len(results)} runs):")
print("-" * 80)

# Header
has_components = any("sq_chi2" in r for r in results)
densities = [r.get("density") for r in results]
has_density = any(d is not None for d in densities)
# Show density column if densities vary across runs
show_density = has_density and len(set(d for d in densities if d is not None)) > 1

def fmt_row(r):
    density_str = f" {r.get('density', float('nan')):>8.4f}" if show_density else ""
    if has_components:
        return (
            f"{r['dir']:<10} {r.get('chi2', float('nan')):>10.4f} "
            f"{r.get('sq_chi2', float('nan')):>10.4f} {r.get('gr_chi2', float('nan')):>10.4f}"
            f"{density_str}"
            f" {r.get('acceptance', float('nan')):>7.1f}% {r.get('best_move', 0):>10d} "
            f"{r.get('seed', 0):>12d}"
        )
    else:
        return (
            f"{r['dir']:<10} {r.get('chi2', float('nan')):>10.6f}"
            f"{density_str}"
            f" {r.get('acceptance', float('nan')):>7.1f}% {r.get('best_move', 0):>10d} "
            f"{r.get('seed', 0):>12d}"
        )

density_hdr = f" {'g/cm3':>8}" if show_density else ""
if has_components:
    print(f"{'Run':<10} {'chi2':>10} {'sq':>10} {'gr':>10}{density_hdr} {'accept%':>8} {'best_move':>10} {'seed':>12}")
else:
    print(f"{'Run':<10} {'chi2':>10}{density_hdr} {'accept%':>8} {'best_move':>10} {'seed':>12}")
print("-" * 80)
for r in results:
    print(fmt_row(r))

# Statistics
chi2_vals = [r["chi2"] for r in results if "chi2" in r]
if len(chi2_vals) > 1:
    print("-" * 80)
    print(f"{'chi2:':<10} min = {min(chi2_vals):.6f}, max = {max(chi2_vals):.6f}, "
          f"mean = {np.mean(chi2_vals):.6f}, std = {np.std(chi2_vals):.6f}")

# Best run
best = results[0]
print(f"\nBest run: {best['dir']} (chi2 = {best.get('chi2', float('nan')):.6f})")

# --- CN comparison ---
cn_results = [(r["dir"], r["cn"]) for r in results if "cn" in r]
if cn_results:
    all_pairs = sorted(set(p for _, cn in cn_results for p in cn))
    print(f"\nCoordination numbers:")
    print("-" * (12 + 10 * len(all_pairs)))
    header = f"{'Run':<12}" + "".join(f"{p:>10}" for p in all_pairs)
    print(header)
    print("-" * (12 + 10 * len(all_pairs)))
    for name, cn in cn_results:
        row = f"{name:<12}" + "".join(f"{cn.get(p, float('nan')):>10.3f}" for p in all_pairs)
        print(row)

    # CN statistics
    if len(cn_results) > 1:
        print("-" * (12 + 10 * len(all_pairs)))
        means = {p: np.mean([cn.get(p, np.nan) for _, cn in cn_results]) for p in all_pairs}
        stds = {p: np.std([cn.get(p, np.nan) for _, cn in cn_results]) for p in all_pairs}
        print(f"{'mean':<12}" + "".join(f"{means[p]:>10.3f}" for p in all_pairs))
        print(f"{'std':<12}" + "".join(f"{stds[p]:>10.3f}" for p in all_pairs))

# --- Overlay plots ---
if HAS_PLT:
    import tomllib

    with open(config_path, "rb") as f:
        config = tomllib.load(f)

    sq_types = []
    for data_type, key in [("xray", "xray_sq"), ("neutron", "neutron_sq")]:
        if config.get("data", {}).get(key) is not None:
            sq_types.append((data_type, key))

    # Check which g(r) data is available
    has_partial_gr = any(
        os.path.isfile(os.path.join(d, "refined_gr.dat")) for d in run_dirs
    )
    has_total_gr = any(
        os.path.isfile(os.path.join(d, "refined_total_gr.dat")) for d in run_dirs
    )
    has_gr_exp = config.get("data", {}).get("xray_gr") is not None

    # Count panels: S(Q) datasets + total g(r) + partial g(r)
    n_panels = len(sq_types) + (1 if has_total_gr else 0) + (1 if has_partial_gr else 0)

    if n_panels > 0:
        fig, axes = plt.subplots(n_panels, 1, figsize=(8, 4 * n_panels), squeeze=False)
        panel = 0

        # --- S(Q) panels ---
        for data_type, key in sq_types:
            ax = axes[panel, 0]
            panel += 1

            # Plot experimental
            exp_path = os.path.join(config_dir, config["data"][key]["file"])
            if os.path.isfile(exp_path):
                sq_exp = np.loadtxt(exp_path)
                ax.plot(sq_exp[:, 0], sq_exp[:, 1], "k--", lw=2, alpha=0.7, label="Experiment")
                q_max = sq_exp[:, 0].max()
            else:
                q_max = 25.0

            # Overlay all runs
            for r in results:
                run_dir = next(d for d in run_dirs if os.path.basename(d) == r["dir"])
                for name in [f"refined_{data_type}_sq.dat", "refined_sq.dat"]:
                    sq_path = os.path.join(run_dir, name)
                    if os.path.isfile(sq_path):
                        sq = np.loadtxt(sq_path, comments="#")
                        mask = sq[:, 0] <= q_max
                        chi2_label = f" ({r.get('chi2', 0):.4f})" if "chi2" in r else ""
                        ax.plot(sq[mask, 0], sq[mask, 1], lw=0.8, alpha=0.6,
                                label=f"{r['dir']}{chi2_label}")
                        break

            label_nice = "X-ray" if data_type == "xray" else "Neutron"
            ax.set(xlabel="Q (1/\u00c5)", ylabel="S(Q)", xlim=(0.3, q_max))
            ax.set_title(f"{label_nice} S(Q) \u2014 ensemble comparison")
            ax.legend(fontsize=7, ncol=2)

        # --- Total g(r) panel ---
        if has_total_gr:
            ax = axes[panel, 0]
            panel += 1

            # Plot experimental g(r)
            if has_gr_exp:
                gr_exp_path = os.path.join(config_dir, config["data"]["xray_gr"]["file"])
                if os.path.isfile(gr_exp_path):
                    gr_exp = np.loadtxt(gr_exp_path)
                    ax.plot(gr_exp[:, 0], gr_exp[:, 1], "k--", lw=2, alpha=0.7, label="Experiment")

            # Overlay total g(r) from all runs
            for r in results:
                run_dir = next(d for d in run_dirs if os.path.basename(d) == r["dir"])
                tgr_path = os.path.join(run_dir, "refined_total_gr.dat")
                if not os.path.isfile(tgr_path):
                    continue
                tgr = np.loadtxt(tgr_path, comments="#")
                chi2_label = f" ({r.get('chi2', 0):.4f})" if "chi2" in r else ""
                ax.plot(tgr[:, 0], tgr[:, 1], lw=0.8, alpha=0.6,
                        label=f"{r['dir']}{chi2_label}")

            ax.axhline(1, color="gray", ls="--", lw=0.5)
            ax.set(xlabel="r (\u00c5)", ylabel="g(r)", xlim=(0, 10))
            ax.set_title("Total X-ray g(r) \u2014 ensemble comparison")
            ax.legend(fontsize=7, ncol=2)

        # --- Partial g(r) panel ---
        if has_partial_gr:
            ax = axes[panel, 0]
            panel += 1

            # Get pair labels from first available file
            pair_labels: list[str] = []
            for d in run_dirs:
                gr_path = os.path.join(d, "refined_gr.dat")
                if os.path.isfile(gr_path):
                    with open(gr_path) as fh:
                        header = fh.readline().strip().lstrip("# ").split()
                    pair_labels = [col.replace("g_", "").replace("_", "") for col in header[1:]]
                    pair_labels = [re.sub(r"([A-Z][a-z]?)([A-Z])", r"\1-\2", lbl) for lbl in pair_labels]
                    break

            # Overlay g(r) from all runs
            colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
            n_pairs = len(pair_labels) if pair_labels else 0
            for ri, r in enumerate(results):
                run_dir = next(d for d in run_dirs if os.path.basename(d) == r["dir"])
                gr_path = os.path.join(run_dir, "refined_gr.dat")
                if not os.path.isfile(gr_path):
                    continue
                gr = np.loadtxt(gr_path, comments="#")
                r_vals = gr[:, 0]
                for pi in range(n_pairs):
                    color = colors[pi % len(colors)]
                    label = f"{pair_labels[pi]}" if ri == 0 else None
                    ax.plot(r_vals, gr[:, 1 + pi], lw=0.6, alpha=0.4, color=color, label=label)

            ax.axhline(1, color="gray", ls="--", lw=0.5)
            ax.set(xlabel="r (\u00c5)", ylabel="g(r)", xlim=(0, 10))
            ax.set_title("Partial g(r) \u2014 ensemble comparison")
            if pair_labels:
                ax.legend(fontsize=7, ncol=2)

        plt.tight_layout()
        outfile = os.path.join(config_dir, "ensemble_comparison.png")
        plt.savefig(outfile, dpi=150, bbox_inches="tight")
        print(f"\nSaved ensemble comparison to {outfile}")
        plt.show()
