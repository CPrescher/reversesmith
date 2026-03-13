# Changelog

All notable changes to rsmith will be documented in this file.

## [1.1.1] - 2026-03-13

### Performance
- **Chirp-Z Transform (CZT)**: replaced O(N_bins × N_Q) sin table lookup with FFT-based
  Bluestein's algorithm. Working set fits in L1 cache (~128 KB) vs sin table spilling to
  L2 (~5 MB). ~25% faster for S(Q)-only, ~16% faster for full configs.
- **Single-pass histogram delta**: combined two cell-list traversals (old/new positions)
  into one pass that computes both distances per neighbor. ~13% faster for S(Q)-only,
  ~7% faster for full configs.
- **Fine constraint cell list**: build a second cell list with cutoff matching the actual
  constraint distances (~4 Å) instead of reusing the RDF cell list (~21 Å). Reduces
  constraint checking from iterating all atoms to only nearby atoms. ~51% faster for
  S(Q)-only, ~35% faster for full configs.

### Added
- New `src/czt.rs` module implementing the Chirp-Z Transform via Bluestein's algorithm
- `rustfft` dependency for FFT computation

## [1.1.0] - 2026-03-12

### Added
- **EPSR (Empirical Potential Structure Refinement)**: outer loop that iteratively refines
  a perturbation potential so the MC simulation naturally reproduces experimental S(Q).
  Produces thermodynamically consistent structures — equilibrium states of a Hamiltonian,
  not just arbitrary arrangements that match the data.
  - New `[epsr]` config section with `iterations`, `feedback`, `smooth_sigma`,
    `moves_per_iteration`, `temperature`, `min_r`, `convergence`, and `ep_restart`
  - Per-iteration output: `epsr_ep_{pair}.dat` (cumulative empirical potential)
  - Per-iteration structure snapshots: `refined_iter_N.xyz`
  - EP restart support: continue refining from a previous EPSR run
- New `src/epsr.rs` module with EP state management, residual decomposition,
  Fourier transform, Gaussian smoothing, and combined potential builder
- `PairPotential::from_vec()` and `PotentialSet::add_potential()` for programmatic
  potential construction
- `Clone` derives on `PairPotential` and `PotentialSet`
- `partial_sq` and `total_sq` fields on `RmcState` for post-run access to structure factors
- `restore_best` flag on `RmcParams` to control best-structure restoration (EPSR needs
  equilibrium S(Q), not biased best-chi2 snapshots)

### Fixed
- Temperature defaulting to 1.0 when annealing was disabled (anneal_start == anneal_end),
  causing ~88% acceptance and structure destruction. Now correctly uses `anneal_end`.
- Premature convergence during annealing: convergence checking is now deferred until after
  the annealing phase completes, with a fresh baseline window.
- Default temperature changed from 1.0 to 0.1 (T=1 almost always increases chi2).

### Documentation
- New EPSR configuration reference (`docs/src/config/epsr.md`)
- Updated RMC docs with temperature guidance, convergence timing, and default changes
- Example `[epsr]` section in `config.toml`

## [1.0.0]

Initial release with RMC refinement, X-ray/neutron S(Q) and g(r) fitting,
Pedone/Buckingham/Coulomb potentials, simulated annealing, coordination
constraints, and structural analysis.
