# Structural Analysis (`--analyze`)

Post-refinement validation of coordination environments and bond angle distributions. RMC is underdetermined — it fits the scattering data but can introduce physically unreasonable local structures. This mode checks whether the refined structure preserves known short-range order (e.g., tetrahedral SiO4 units).

## Usage

```bash
# Compare starting structure vs. refined structure (both analyzed automatically)
reversesmith config.toml --analyze

# Analyze a specific structure only
reversesmith config.toml --analyze some_structure.xyz
```

**Without an explicit path**, `--analyze` runs the analysis twice:
1. The starting structure defined in `[system]`
2. `refined.xyz` in the config directory (if it exists)

This gives a side-by-side comparison of how the RMC refinement changed the local structure.

**With an explicit path** (`--analyze path.xyz`), only that structure is analyzed. The format is auto-detected from the file extension (`.xyz` or `.data`/`.lmp`).

The config file is always required — it defines which pair cutoffs and species to use.

## Configuration

Add an optional `[analysis]` section to the config TOML:

```toml
[analysis.cutoffs]       # pair cutoffs for CN and angle analysis
"Si-O" = 2.2
"Ca-O" = 3.0
```

If `[analysis.cutoffs]` is not specified, the cutoffs from `[[constraints.coordination]]` are used as fallback. Both sources are merged, with explicit `[analysis.cutoffs]` taking precedence.

### Optional settings

```toml
[analysis]
angle_bins = 180                    # bins for 0-180 degree range (default: 180)
angle_triplets = ["O-Si-O"]        # only compute these triplets (default: all)
```

The `angle_triplets` filter uses canonical labels with sorted end-species: `"O-Si-O"`, not `"Si-O-O"`.

### Minimal example (no `[analysis]` section needed)

If your config already has coordination constraints:

```toml
[[constraints.coordination]]
pair = "Si-O"
min = 3
max = 6
cutoff = 2.2
```

Then `--analyze` will use the Si-O cutoff of 2.2 A automatically. Add `[analysis.cutoffs]` only if you need additional pairs or different cutoffs.

## What it computes

### Coordination numbers

For each pair A-B with cutoff r_c: count the number of B neighbors within r_c of each A atom. Reports:

- **Histogram**: count and fraction of A atoms with CN = 0, 1, 2, ...
- **Statistics**: mean, standard deviation, min, max

### Bond angles

1. Build a neighbor list from all defined pairs (both directions, e.g., Si-O also gives O bonded to Si).
2. For each atom M with at least 2 neighbors X and Y: compute the angle X-M-Y via `acos(dot(dr_MX, dr_MY) / (|dr_MX| * |dr_MY|))` using minimum-image PBC.
3. Classify each angle by a canonical triplet label with sorted end-species (e.g., "O-Si-O", "Ca-O-Si").
4. Bin into a histogram and normalize to probability density P(theta).

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

### File formats

**CN file** (`analysis_*_cn.dat`):

```
# CN  count_Si-O  frac_Si-O  count_Ca-O  frac_Ca-O
0  0  0.000000  0  0.000000
1  0  0.000000  0  0.000000
...
4  1997  0.998500  44  0.022000
5  3  0.001500  523  0.261500
```

**Angle file** (`analysis_*_angles.dat`):

```
# angle_deg  P_O-Si-O  P_Si-O-Si  P_Ca-O-Si  ...
0.50  0.000000  0.000000  0.000000
1.50  0.000000  0.000000  0.000000
...
109.50  0.025431  0.000000  0.002341
```

The angle columns are probability densities normalized so that the integral over 0-180 degrees equals 1.

## Interpreting results

### Expected values for CaSiO3 glass

| Property | Expected | Indicates problem if... |
|----------|----------|------------------------|
| Si-O CN | 4.0 (tetrahedral) | significant 3- or 5-fold fraction |
| O-Si-O angle | ~109.5 deg (tetrahedral) | peak shifts below ~105 deg |
| Si-O-Si angle | ~140-150 deg (bridge angle) | bimodal or very broad |
| Ca-O CN | ~6 (variable) | unphysically low (<3) |

### Comparing starting vs. refined

- If the starting structure (e.g., from MD) shows correct coordination and the refined structure does not, the RMC may be over-fitting noise or the constraints are too loose.
- Large changes in CN distributions suggest the refinement needs tighter coordination constraints.
- Angle distributions broadening significantly may indicate the RMC is distorting polyhedra to fit high-Q features.
