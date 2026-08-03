# EPSR Configuration

For validated native-EPSR comparisons, convergence results, and the limits of
the current scientific claims, see
[Validation and Reference-Code Benchmarks](../validation.md).

The `[epsr]` section enables Empirical Potential Structure Refinement (Soper, 1996). When present, an outer loop wraps the RMC refinement, iteratively refining a perturbation potential so the simulation naturally reproduces the experimental data.

The aim is to sample configurations from a reference-plus-empirical
Hamiltonian while improving agreement with the data, rather than accepting
arbitrary coordinate changes from the data residual alone. Finite sampling,
the chosen reference potential, and the evolving empirical potential still
need independent structural validation; enabling EPSR does not by itself prove
thermodynamic or chemical correctness.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mode` | string | `"hybrid"` | `"hybrid"` combines data and energy acceptance; `"pure"` performs energy-only equilibrium MC followed by an EPSR update |
| `iterations` | integer | 10 | Number of outer EPSR iterations |
| `feedback` | float | 0.2 | Feedback factor for EP update: `EP += feedback * kT * Δg(r)` |
| `smooth_sigma` | float | 0.02 | Gaussian smoothing width (Å) applied to EP |
| `moves_per_iteration` | integer | from `[rmc]` | MC moves per EPSR epoch |
| `temperature` | float | 0.025 | kT in eV for EP update (~300K) |
| `min_r` | float | 1.0 | Zero EP below this distance (Å) |
| `convergence` | float | 0.0 | Stop if relative EP change (max |ΔEP| / max |EP|) falls below this for `convergence_window` consecutive iterations. 0 = run all iterations |
| `convergence_window` | integer | 3 | Number of consecutive iterations that must satisfy the convergence criterion before stopping |
| `ep_restart` | string | — | Directory containing previous `epsr_ep_{pair}.dat` files to seed the EP |

## Example

```toml
[epsr]
mode = "pure"
iterations = 20
feedback = 0.2
smooth_sigma = 0.02
moves_per_iteration = 200_000
temperature = 0.025
min_r = 1.0
convergence = 0.05
convergence_window = 3
```

## Algorithm

Each EPSR outer iteration:

1. Build combined potential: `V_ref(r) + EP(r)`
2. Run the configured hybrid or pure MC epoch for `moves_per_iteration` moves
3. Compute residual: `ΔS(Q) = S_exp(Q) - S_sim(Q)`
4. Decompose to partials via proportional weighting: `ΔS_ab(Q) = w_ab(Q) * ΔS(Q) / Σ w_cd(Q)²`
5. Sine transform each `ΔS_ab(Q)` to `Δg_ab(r)`
6. Update: `EP_ab(r) += feedback * kT * Δg_ab(r)`
7. Gaussian-smooth EP, zero below `min_r`
8. Check convergence and repeat

## Output files

Per iteration:
- `epsr_ep_{pair}.dat` — cumulative empirical potential (r in Å, V in eV)

Log output per iteration:
```
EPSR iter N: chi2 = X, max |ΔEP| = Y eV, acceptance = Z%
```

## Pure and hybrid modes

`mode = "pure"` follows the original EPSR separation: every epoch samples the
combined reference-plus-empirical potential using energy-only Metropolis MC,
then calculates S(Q) from the resulting configuration and updates the empirical
potential. Pure mode uses an exact cutoff-plus-skin Verlet neighbor list for
the analytical potential trials.

`mode = "hybrid"` retains rsmith's data-plus-energy RMC acceptance. Its first
iteration uses the full `[rmc]` settings, including annealing, convergence
checking, and optional best-structure restoration. Subsequent iterations use
equilibrium settings:

- **No best-structure restoration**: The EP update needs the equilibrium S(Q), not a biased low-chi2 snapshot, because the goal is to shift the equilibrium itself.

- **No early convergence stopping**: Each iteration runs the full `moves_per_iteration` moves. The system needs to equilibrate under the combined potential before the residual S(Q) is meaningful.

- **No annealing**: Iterations start from the equilibrated structure at `anneal_end` temperature.

As a result, the chi2 reported for hybrid iterations 2+ is the **equilibrium
chi2** under the current combined potential, not the best fluctuation. In pure
mode, the reported post-epoch scattering residual is diagnostic: it does not
participate directly in accepting the moves within that epoch.

## Notes

- EPSR uses the first X-ray S(Q) dataset when available and otherwise falls
  back to the first neutron S(Q) dataset
- The reference potential from `[potential]` is preserved; EP is added on top
- If no `[potential]` section exists, EP alone drives the simulation
- The `feedback` parameter controls how aggressively the EP adapts; values of 0.1–0.3 are typical
- Use `ep_restart` to continue from a previous EPSR run: the EP files are loaded and further refined

## Restarting from a previous run

To continue refining EP from a previous EPSR run, point `ep_restart` to the output directory:

```toml
[epsr]
iterations = 10
feedback = 0.2
ep_restart = "../previous_run"  # loads epsr_ep_Ca-O.dat, etc.
```

The EP tables are interpolated onto the current grid and accumulated. Missing pair files are skipped (that pair starts from zero).
