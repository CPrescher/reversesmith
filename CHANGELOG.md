# Changelog

All notable changes to rsmith will be documented in this file.

## Unreleased

### Added
- Lennard-Jones 12-6 reference pair potentials with explicit pair parameters,
  shifted-cutoff tabulation, validation, and LAMMPS reference tests.
- Pure EPSR refinement mode validated against the EPSR26 LiquidGa50C worked
  example, including deterministic forward parity and a ten-seed stochastic
  comparison.
- Cutoff-plus-skin Verlet neighbor lists for analytical pair-potential Monte
  Carlo, with conservative rebuild tracking and brute-force periodic tests.
- Release-mode GAP/QUIP atom-count benchmark using the same diamond-Si
  structures and rejected single-atom trial path as the SNAP and MACE examples.

### Fixed
- Convert experimental `S(Q)`, `i(Q) = S(Q)-1`, and
  `F(Q) = Q[S(Q)-1]` conventions back to internal `S(Q)` before empirical-
  potential updates.

### Performance
- Reduced the ten-seed LiquidGa50C single-thread median from 44.64 s to
  16.71 s without changing the accepted trajectories or scientific metrics.
  On the same Apple M4 Pro, rsmith is 1.740x faster than native EPSR26 for the
  fixed one-million-move protocol.

### Documentation
- Added the EPSR benchmark ladder, local-only workshop-data import workflow,
  frozen LiquidGa results, and the silica acquisition/publication-permission
  protocol.
- Documented current GAP/QUIP whole-system evaluation scaling and published Si
  GAP benchmark results alongside native SNAP and incremental MACE.

## [1.4.0] - 2026-08-02

### Added
- **Native SNAP backend**: evaluate LAMMPS/FitSNAP `.snapcoeff` and `.snapparam`
  models directly in Rust without a LAMMPS runtime dependency.
- Support for linear, explicit multi-element (`chemflag`), and quadratic SNAP
  energy models, including the standard switching, normalization, self-weight,
  and isolated-atom conventions.
- **MACE/Python backend**: persistent Python worker integration for using MACE
  checkpoints as RMC energy regularizers, with CPU/GPU device selection,
  PyTorch thread control, and optional float32 conversion.
- Three MACE energy-delta strategies: broadly compatible full-system evaluation,
  exact dense-local evaluation, and exact incremental message-passing inference
  for supported short-range `MACE` and `ScaleShiftMACE` checkpoints.
- Release-mode MACE and native-SNAP scaling examples that emit CSV timing data
  for atom-count, thread-count, precision, and delta-mode comparisons.

### Performance
- Native SNAP caches accepted per-atom energies and recomputes only local
  environments affected by a proposed single-atom move.
- Dense-local MACE clusters are built with a dedicated Rust cell list and sent
  to the Python worker as compact atom-index and periodic-image arrays.
- Dense-local MACE caches accepted per-atom energies, requiring only one
  energy-only forward pass over the proposed cluster per trial.
- Incremental MACE caches accepted hidden features and recomputes only the
  causal subgraph changed at each message-passing layer. Its measured cost is
  effectively independent of total atom count once the local environment has
  stabilized.

### Validation
- Added frozen LAMMPS reference fixtures for standard, chemical, and quadratic
  SNAP models, plus rotation-invariance and transactional-cache tests.
- Added real-checkpoint MACE integration coverage comparing full, local, and
  incremental deltas across rejected and accepted moves.
- Hardened SNAP parsing and radial-parameter validation to reject inconsistent
  or non-finite model data instead of silently evaluating it.

### Documentation
- Documented native SNAP setup, supported model conventions, and numerical
  reference testing.
- Added a MACE mode-selection guide covering compatibility, exactness,
  atom-count scaling, precision, compilation, and measured CPU performance.

## [1.3.0] - 2026-07-08

### Added
- **ML-potential-regularized RMC**: new `[ml_potential]` config section for using
  machine-learning interatomic potentials as RMC energy regularizers.
- Optional `gap-quip` feature with a GAP/QUIP backend using local affected-atom
  per-atom energy deltas through a small C ABI shim.
- Generic `EnergyModel` abstraction so pair potentials and ML potentials share the
  RMC energy-regularizer path.
- External QUIP build helper that writes both `env.sh` and pkg-config metadata for
  portable linking against a local QUIP/libAtoms build.
- Feature-gated GAP/QUIP integration test and mock energy-model RMC acceptance tests.

### Changed
- Existing pair potentials remain config-compatible but now implement the generic
  energy regularizer interface.
- Default builds remain QUIP-free; QUIP is only required when building with
  `--features gap-quip`.
- CI now runs the default Rust test suite on both Linux and macOS.

### Documentation
- Added GAP/QUIP build instructions for macOS, Linux, and cluster/module systems.
- Added ML-potential configuration reference and example `[ml_potential]` config.

## [1.2.0] - 2026-03-16

### Added
- **f(r) fitting**: new `[data.xray_fr]` config option for fitting the reduced pair
  distribution function f(r) = (2/π) ∫ Q·[S(Q)-1]·W(Q)·sin(Qr) dQ. Unlike g(r), f(r)
  does not depend on number density ρ₀, making it suitable for comparing ensemble runs
  at different densities (`--density`). Mutually exclusive with `[data.xray_gr]`.
- f(r) panel in `compare_ensemble.py` with experimental f(r) computed from S(Q)
- Documentation for f(r) in config reference, output files, and scattering algorithms

## [1.1.2] - 2026-03-13

### Performance
- **CZT for g(r) inverse FT**: replaced the dense FT matrix multiply (n_r × n_Q FMAs)
  with a CZT sine transform for the Q → r inverse Fourier transform used in g(r) fitting.
  ~38% faster for full configs (S(Q)-only unaffected).

### Fixed
- CZT post-chirp phase formula generalized to handle input grids not starting at dR/2
  (bin-center convention). The previous formula was correct for S(Q) but incorrect for
  the g(r) inverse transform where the Q-grid starts at 0.

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
