# Structural Analysis

Post-refinement validation of coordination environments and bond angle distributions. RMC is underdetermined -- it fits the scattering data but can introduce physically unreasonable local structures. The `--analyze` mode checks whether the refined structure preserves known short-range order.

## Usage

```bash
# Compare starting vs. refined structure (automatic)
rsmith config.toml --analyze

# Analyze a specific structure
rsmith config.toml --analyze some_structure.xyz
```

**Without an explicit path**, `--analyze` runs the analysis twice:
1. The starting structure defined in `[system]`
2. `refined.xyz` in the config directory (if it exists)

This gives a side-by-side comparison of how the refinement changed the local structure.

**With an explicit path**, only that structure is analyzed. The format is auto-detected from the file extension (`.xyz` or `.data`/`.lmp`).

The config file is always required -- it defines the pair cutoffs and species.

## Configuration

See [`[analysis]` config reference](config/analysis.md) for all options. If your config already has `[[constraints.coordination]]` blocks, `--analyze` uses their cutoffs automatically.

## What it computes

### Coordination numbers

For each pair A-B with cutoff r_c: count the number of B neighbours within r_c of each A atom.

Reports:
- **Histogram**: count and fraction of A atoms with CN = 0, 1, 2, ...
- **Statistics**: mean, standard deviation, min, max

### Bond angles

1. Build a neighbour list from all defined pairs
2. For each atom M with at least 2 neighbours X and Y: compute the angle X-M-Y using minimum-image PBC
3. Classify by canonical triplet label (sorted end-species, e.g., "O-Si-O")
4. Bin into a histogram, normalize to probability density P(theta)

Reports: number of angles, peak angle, mean angle for each triplet.

## Output files

Without explicit path (dual mode):

| File | Description |
|------|-------------|
| `analysis_starting_cn.dat` | CN histograms for starting structure |
| `analysis_starting_angles.dat` | Angle distributions for starting structure |
| `analysis_refined_cn.dat` | CN histograms for refined structure |
| `analysis_refined_angles.dat` | Angle distributions for refined structure |

With explicit path:

| File | Description |
|------|-------------|
| `analysis_analysis_cn.dat` | CN histograms |
| `analysis_analysis_angles.dat` | Angle distributions |

### CN file format

```
# CN  count_Si-O  frac_Si-O  count_Ca-O  frac_Ca-O
0  0  0.000000  0  0.000000
4  1997  0.998500  44  0.022000
5  3  0.001500  523  0.261500
```

### Angle file format

```
# angle_deg  P_O-Si-O  P_Si-O-Si  P_Ca-O-Si
0.50  0.000000  0.000000  0.000000
109.50  0.025431  0.000000  0.002341
```

Columns are probability densities normalized so the integral over 0--180 degrees equals 1.

## Interpreting results

### Expected values for CaSiO3 glass

| Property | Expected | Problem if... |
|----------|----------|---------------|
| Si-O CN | 4.0 (tetrahedral) | Significant 3- or 5-fold fraction |
| O-Si-O angle | ~109.5 deg | Peak shifts below ~105 deg |
| Si-O-Si angle | ~140--150 deg | Bimodal or very broad |
| Ca-O CN | ~6 (variable) | Unphysically low (<3) |

### Comparing starting vs. refined

- **CN distribution unchanged** -- good, RMC preserved local order
- **CN distribution broadened** -- constraints may be too loose, or potentials too weak
- **Angle distributions sharpened** -- the refinement improved agreement with experimental short-range order
- **Angle distributions broadened significantly** -- the RMC is distorting polyhedra to fit high-Q features; consider tighter constraints or higher potential weight

This comparison is the primary way to validate that the hybrid RMC (potential + constraints) is working correctly.
