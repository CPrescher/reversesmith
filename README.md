# Reversesmith

A Reverse Monte Carlo (RMC) structure refinement tool written in Rust. Refines atomic structures against experimental X-ray and neutron scattering data while enforcing physical constraints.

## Overview

Reversesmith iteratively displaces atoms in a model structure to minimize the difference between computed and experimental scattering functions. It supports fitting against:

- **S(Q)** — total X-ray or neutron structure factor in reciprocal space
- **g(r)** — total X-ray pair distribution function in real space

Both targets can be fitted simultaneously with independent weights and fit ranges.

## Quick Start

```bash
cargo build --release
./target/release/reversesmith config.toml
```

To compute S(Q) and g(r) from a structure without refinement:

```bash
./target/release/reversesmith config.toml --compute-sq-only
```

## Configuration File

The input is a single TOML file. All file paths are resolved relative to the config file location.

### `[system]` — Structure Input

```toml
[system]
structure = "structure.data"   # Path to atomic structure file
format = "lammps"              # "lammps" or "xyz"

[system.types]                 # Required for LAMMPS format: type ID -> element
1 = "Ca"
2 = "Si"
3 = "O"
```

**LAMMPS format**: Reads the `Atoms # charge` section. Expects columns: `atom_id type charge x y z [ix iy iz]`. Box bounds from `xlo xhi` / `ylo yhi` / `zlo zhi` lines.

**XYZ format**: Standard XYZ with box in the comment line, either as `Lattice="Lx 0 0 0 Ly 0 0 0 Lz"` (extended XYZ) or `Lx Ly Lz`.

### `[data]` — Experimental Data

Each dataset is a two-column whitespace-separated file (lines starting with `#` are skipped).

#### `[data.xray_sq]` — X-ray Structure Factor

```toml
[data.xray_sq]
file = "experimental.sq"    # Two columns: Q (1/A), S(Q)
weight = 1.0                # Relative weight in total chi2 (default: 1.0)
sigma = 0.02                # Experimental uncertainty (default: 0.01)
fit_min = 0.5               # Only fit Q > fit_min (default: 0, fit all)
fit_max = 18.0              # Only fit Q < fit_max (default: inf, fit all)
```

#### `[data.neutron_sq]` — Neutron Structure Factor

Same format as `xray_sq`. Uses coherent neutron scattering lengths instead of X-ray form factors.

#### `[data.xray_gr]` — X-ray Pair Distribution Function

```toml
[data.xray_gr]
file = "experimental.gr"    # Two columns: r (A), g(r)
weight = 0.3                # Relative weight in total chi2 (default: 1.0)
sigma = 0.02                # Experimental uncertainty (default: 0.01)
fit_min = 0.0               # Only fit r > fit_min (default: 0)
fit_max = 7.0               # Only fit r < fit_max (default: inf)
qmax = 17.97                # Q_max for inverse FT (default: experimental S(Q) Qmax)
lorch = true                # Apply Lorch window in Q-space (default: true)
```

**Important**: The model g(r) is computed by inverse Fourier transform of the model's total X-ray S(Q). The `qmax` and `lorch` parameters must match those used to derive the experimental g(r) from experimental S(Q), otherwise peak heights will be inconsistent (see [g(r) Consistency](#gr-consistency) below).

### `[rmc]` — Refinement Parameters

```toml
[rmc]
max_moves = 1_000_000       # Total MC moves to attempt
max_step = 0.1              # Initial maximum atomic displacement (A)
checkpoint_every = 50_000   # Save checkpoint every N moves
seed = 42                   # Random number generator seed
print_every = 1000          # Print status every N moves
target_acceptance = 0.3     # Target acceptance rate for step adaptation
adjust_step_every = 5000    # Adjust step size every N moves

# Simulated annealing (omit for standard RMC at T=1)
anneal_start = 2.0          # Starting temperature
anneal_end = 0.01           # Final temperature
anneal_steps = 200_000      # Moves over which to anneal (rest of run at anneal_end)

# Early stopping
convergence_threshold = 1e-4  # Stop if chi2 improves less than this...
convergence_window = 50_000   # ...over this many moves (0 = disabled)
```

### `[sq]` — S(Q) Computation Parameters

```toml
[sq]
qmin = 0.0                  # Minimum Q for model S(Q) grid (1/A)
qmax = 20.0                 # Maximum Q (1/A)
nq = 500                    # Number of Q grid points
lorch = true                # Apply Lorch window in r-space for S(Q) computation
rdf_cutoff = 11.0           # Maximum r for RDF histogram (A)
rdf_nbins = 550             # Number of RDF histogram bins
```

The model computes S(Q) on a uniform grid from `qmin` to `qmax`. This grid should extend beyond the experimental Q range to avoid edge effects.

### `[constraints]` — Physical Constraints

Moves that violate constraints are rejected before evaluating chi2.

#### Minimum Distances

```toml
[constraints.min_distance]
"Ca-O" = 1.8
"Si-O" = 1.2
"O-O" = 2.0
"Ca-Si" = 2.2
"Ca-Ca" = 2.3
"Si-Si" = 2.3
```

Pair order does not matter (`"Ca-O"` and `"O-Ca"` are equivalent).

#### Coordination Constraints

```toml
[[constraints.coordination]]
pair = "Si-O"
min = 3          # Minimum coordination number
max = 6          # Maximum coordination number
cutoff = 2.2     # Neighbour search radius (A)
```

Multiple `[[constraints.coordination]]` blocks can be defined. After a proposed move, the coordination number of the moved atom (and its neighbours) is checked against the bounds.

## Output Files

| File | Description |
|------|-------------|
| `refined.xyz` | Refined atomic structure (extended XYZ with lattice) |
| `refined_sq.dat` | Final computed total X-ray S(Q) |
| `refined_total_gr.dat` | Final computed total X-ray g(r) via inverse FT |
| `checkpoint.dat` | Checkpoint for resuming (text format) |

With `--compute-sq-only`:

| File | Description |
|------|-------------|
| `computed_sq.dat` | Total X-ray S(Q) of input structure |
| `computed_gr.dat` | All partial g(r) functions |
| `computed_total_gr.dat` | Total X-ray g(r) via inverse FT |

## Algorithms

### RMC Algorithm

The core algorithm follows McGreevy & Pusztai (1988). Each iteration:

1. **Select** a random atom and propose a displacement drawn uniformly from [-max_step, +max_step] in each dimension.
2. **Check constraints**: reject immediately if minimum distance or coordination constraints are violated.
3. **Compute chi2**: incrementally update the structure factor and evaluate the fit to experimental data.
4. **Accept/reject** via the Metropolis criterion:
   - If chi2 decreases: always accept.
   - If chi2 increases by Delta-chi2: accept with probability `P = exp(-Delta-chi2 / (2T))`.

At standard RMC temperature T=1, this samples configurations consistent with the data within the experimental uncertainty. For optimization (driving chi2 as low as possible), use T < 1 via simulated annealing.

### Incremental S(Q) Updates

Recomputing S(Q) from scratch after each atom move would be O(N_atoms^2). Instead, reversesmith uses incremental updates:

1. **Histogram delta**: Compute the RDF histogram contribution of the moved atom at its old and new positions. Only histogram bins that change contribute to the S(Q) update.

2. **Partial S(Q) delta**: For each changed bin, compute the contribution to each partial S_ab(Q) via the precomputed sin(Q*r) lookup table:

```
Delta-S_ab(Q_k) = (4*pi*rho0*dr / Q_k) * sum_i [r_i * W(r_i) * Delta-g_ab(r_i) * sin(Q_k * r_i)]
```

3. **Total S(Q)**: Reweight and sum: `S_X(Q) = sum_{ab} w_ab(Q) * S_ab(Q)`

This reduces the per-move cost from O(N^2 * N_Q) to O(N * Delta-bins * N_Q), where Delta-bins is typically small.

### S(Q) Computation (Forward Transform)

Partial structure factors S_ab(Q) are computed from partial g_ab(r) via:

```
S_ab(Q) = 1 + (4*pi*rho0 / Q) * integral[r * (g_ab(r) - 1) * W(r) * sin(Q*r) dr]
```

where W(r) is the optional Lorch modification function `sin(pi*r/r_max) / (pi*r/r_max)` that suppresses truncation ripples.

For the initial computation and `--compute-sq-only`, a DST-II (Discrete Sine Transform Type II) via FFT is used for efficiency, with 8x zero-padding for accurate interpolation onto the requested Q grid.

During RMC, the same formula is evaluated by direct summation using precomputed sin(Q*r) tables, enabling incremental updates.

### X-ray Weighting (Faber-Ziman Convention)

The total X-ray structure factor combines partial S_ab(Q) with Q-dependent weights:

```
S_X(Q) = sum_{a<=b} (2 - delta_ab) * c_a * c_b * f_a(Q) * f_b(Q) / <f(Q)>^2 * S_ab(Q)
```

where:
- `c_a` = concentration (atomic fraction) of species a
- `f_a(Q)` = X-ray atomic form factor (Cromer-Mann parameterization)
- `<f(Q)>` = concentration-weighted average form factor
- `delta_ab` = 1 if a=b, 0 otherwise

Neutron weighting uses Q-independent coherent scattering lengths `b_a` instead of form factors.

### g(r) from S(Q) (Inverse Transform)

The model g(r) for comparison with experimental data is computed by inverse Fourier transform of the total X-ray S(Q):

```
g_X(r) = 1 + (dq / (2*pi^2 * rho0 * r)) * sum_k Q_k * W(Q_k) * (S_X(Q_k) - 1) * sin(Q_k * r)
```

The sum runs over Q points up to `qmax` (matching the experimental data range). An optional Q-space Lorch window `W(Q) = sin(pi*Q/Q_max) / (pi*Q/Q_max)` damps high-Q contributions to reduce termination ripples.

This is precomputed as a matrix-vector multiplication: `g(r_i) = 1 + sum_k M[i][k] * (S(Q_k) - 1)`, enabling efficient incremental updates when S(Q) changes.

### g(r) Consistency

The model g(r) is derived from the model S(Q) via inverse Fourier transform. The experimental g(r) is typically derived from experimental S(Q) via the same type of transform. For the chi2 comparison to be meaningful, **both transforms must use identical parameters**:

- **Q_max**: The maximum Q used in the FT. Set via the `qmax` option on `[data.xray_gr]`. Defaults to the maximum Q of the experimental S(Q) data.
- **Lorch window**: Whether a Q-space Lorch function is applied. Set via the `lorch` option on `[data.xray_gr]`. Defaults to `true`.
- **Number density rho0**: Computed from the model structure.

If these parameters don't match those used to produce the experimental g(r), peak heights will differ systematically — the model g(r) may have sharper or broader peaks than the experimental g(r), creating a bias that distorts the refinement.

### Simulated Annealing

Standard RMC uses T=1, which samples configurations consistent with the data rather than optimizing. For structure refinement starting from a good initial model (e.g., from MD simulation), simulated annealing is more effective:

**Exponential cooling schedule:**

```
T(n) = T_start * (T_end / T_start) ^ (n / anneal_steps)
```

For n > anneal_steps, T remains at T_end.

| Phase | Temperature | Behavior |
|-------|-------------|----------|
| Early (T > 1) | High | Explores broadly, escapes local minima |
| Middle (T ~ 0.1-1) | Moderate | Selective exploration |
| Late (T < 0.1) | Low | Greedy optimization, only improvements accepted |

**Step size coupling**: During annealing, the maximum displacement scales with temperature (clamped to max_step): `step = max_step * min(T, 1.0)`. After annealing completes, the step adapts based on acceptance rate (targeting ~30%).

**Practical guidelines:**
- `anneal_start = 1.0-2.0`: Higher values explore more aggressively.
- `anneal_end = 0.01-0.1`: Lower values optimize more strictly.
- `anneal_steps`: Controls the fraction of the run spent annealing. Typically 10-30% of `max_moves`. The remaining moves optimize at `anneal_end`.

### Best-Structure Tracking

The RMC tracks the atom positions corresponding to the lowest chi2 encountered during the entire run. At completion, the best structure is restored and saved. This prevents the output from being worse than intermediate states, which can happen when T fluctuates or when the system random-walks at high temperature.

### Adaptive Step Size

The maximum atomic displacement adjusts during the run:

- **During annealing**: Step scales linearly with temperature (capped at `max_step`). Smaller T means smaller, more refined moves.
- **After annealing**: Step adapts based on acceptance rate. If acceptance > target + 5%, step increases by 5%. If acceptance < target - 5%, step decreases by 5%.

### Chi-squared

The total chi2 is the sum of weighted contributions from all datasets:

```
chi2_total = sum_datasets [weight * sum_i ((y_calc(x_i) - y_exp(x_i))^2 / sigma_i^2)]
```

where the sum over i runs only over data points within the specified `fit_min` to `fit_max` range. The `weight` parameter controls the relative influence of each dataset. The `sigma` parameter sets the per-point experimental uncertainty (uniform across points within a dataset).

**Balancing S(Q) and g(r)**: Since g(r) is derived from S(Q), they contain correlated information. The `weight` parameter on g(r) should typically be small (0.01-0.3) to avoid one target overwhelming the other.

### Constraints

Constraints are evaluated before computing chi2, so invalid moves are rejected cheaply:

- **Minimum distances**: No atom pair of a given type can be closer than the specified distance. O(N) per move.
- **Coordination numbers**: The number of type-B neighbors within a cutoff of each type-A atom must be within [min, max]. Only atoms near the moved atom are rechecked.

### Convergence Detection

If `convergence_threshold > 0`, the RMC monitors chi2 improvement over sliding windows of `convergence_window` moves. If the improvement is less than the threshold, the run stops early. During annealing, convergence checking is deferred until T reaches `anneal_end`.

## Supported Elements

X-ray form factors (Cromer-Mann): Ca, Si, O, Na, Al, Mg

Neutron scattering lengths: Ca, Si, O, Na, Al, Mg

Additional elements can be added in `src/xray.rs`.

## Dependencies

- **rand** — Random number generation
- **serde** + **toml** — Configuration file parsing
- **rayon** — Parallel computation (histograms, S(Q))
- **rustfft** — FFT for DST-II based S(Q) computation

## Example Configuration

```toml
[system]
structure = "glass.data"
format = "lammps"

[system.types]
1 = "Ca"
2 = "Si"
3 = "O"

[data.xray_sq]
file = "experimental.sq"
weight = 1.0
sigma = 0.02
fit_min = 0.5

[data.xray_gr]
file = "experimental.gr"
weight = 0.3
sigma = 0.02
fit_max = 7.0

[rmc]
max_moves = 1_000_000
max_step = 0.1
checkpoint_every = 50_000
seed = 42
print_every = 1000
anneal_start = 2.0
anneal_end = 0.01
anneal_steps = 200_000
convergence_threshold = 1e-4
convergence_window = 50_000

[sq]
qmin = 0
qmax = 20.0
nq = 500
lorch = true
rdf_cutoff = 11.0
rdf_nbins = 550

[constraints.min_distance]
"Ca-O" = 1.8
"Si-O" = 1.2
"O-O" = 2.0
"Ca-Si" = 2.2
"Ca-Ca" = 2.3
"Si-Si" = 2.3

[[constraints.coordination]]
pair = "Si-O"
min = 3
max = 6
cutoff = 2.2
```

## References

- McGreevy, R.L. & Pusztai, L. (1988). Reverse Monte Carlo Simulation: A New Technique for the Determination of Disordered Structures. *Mol. Simul.*, 1, 359-367.
- Keen, D.A. (2001). A comparison of various commonly used correlation functions for describing total scattering. *J. Appl. Cryst.*, 34, 172-177.
- Waasmaier, D. & Kirfel, A. (1995). New analytical scattering-factor functions for free atoms and ions. *Acta Cryst.*, A51, 416-431.
