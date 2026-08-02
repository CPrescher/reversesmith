# Native SNAP Implementation

The native backend reproduces the energy-evaluation part of LAMMPS SNAP in
Rust. LAMMPS is a numerical reference in the test suite, not a library or
subprocess used by rsmith. The implementation is divided into model-file
loading, descriptor mathematics, coefficient contraction, and a stateful local
RMC runtime.

## Evaluation flow

For initialization, rsmith performs the following operations:

1. Parse the `.snapparam` and `.snapcoeff` files and validate all dimensions,
   finite values, flags, units metadata, and coefficient counts.
2. Map each rsmith atom type to an element block by its chemical symbol.
3. Precompute the independent bispectrum indices and Clebsch--Gordan coupling
   tables for `twojmax`.
4. Derive all species-pair cutoffs and construct a SNAP-specific periodic cell
   list using their maximum.
5. Evaluate and cache every accepted-state atomic energy and their sum.

For every trial move, only local state is evaluated:

1. Find the moved atom and every atom whose environment contains the moved atom
   at either its old or new position.
2. Sum the cached old energies for that affected set.
3. Recompute the same atomic environments with the trial position and subtract
   the old sum.
4. Keep the new per-atom values in a pending transaction.
5. On acceptance, update the position, cell-list bucket, affected energy-cache
   entries, and total energy. On rejection, discard the pending transaction.

The old and new neighbor shells are both required: an atom can leave one local
environment while entering another. Periodic coordinates are wrapped before
cell lookup, and pair distances use the minimum-image convention.

## Descriptors

SNAP maps each three-dimensional neighbor vector onto the three-sphere. The
native implementation builds Cayley--Klein parameters from the neighbor
direction and radial angle, then constructs Wigner `U^j` matrices with a
recurrence. Angular momenta are stored as doubled integers, so half-integer
quantum numbers do not require floating-point comparisons.

The neighbor density for a central atom is accumulated as

```text
u^j = u^j_self + sum_n f_c(r_n) * w_n * U^j(r_n)
```

where `w_n` is the coefficient-file element weight and `f_c` combines the
configured outer and optional inner radial switches. The species-pair cutoff is
`rcutfac * (R_c + R_n)`.

Clebsch--Gordan tables couple pairs of density matrices to produce the
bispectrum components. Their independent `(j1, j2, j)` ordering matches the
LAMMPS pyramidal ordering used by `.snapcoeff` files. `bnormflag` applies the
LAMMPS normalization, and `bzeroflag` removes the isolated-atom contribution.

With `chemflag`, rsmith maintains one neighbor density per model element and
emits the bispectrum inside each ordered `(e1, e2, e3)` channel. This order,
including the `wselfallflag` self-density convention, is the coefficient order
written by FitSNAP and consumed by LAMMPS.

## Atomic and total energy

For `K` descriptors, a linear atomic model is

```text
E_i = beta_0 + sum_k beta_k * B_ik
```

When `quadraticflag` is enabled, the coefficient file stores the upper triangle
of the symmetric quadratic matrix row by row. rsmith evaluates

```text
E_i = beta_0 + beta^T B_i + 1/2 B_i^T A B_i
```

Diagonal stored coefficients retain the factor `1/2`; off-diagonal products
have a factor of one because the symmetric matrix contains them twice. The
total SNAP energy is `E = sum_i E_i`.

## Local RMC cache

SNAP is strictly local, so a single-atom move changes only the atomic energy of
the moved atom and central atoms that include it as a neighbor. If `A` is the
union of the affected old and new cutoff shells, the backend computes

```text
delta_E = sum_(i in A) E_i_trial - sum_(i in A) E_i_accepted
```

All unaffected terms cancel and remain cached. A dedicated cell list supplies a
coarse maximum-cutoff candidate set; exact element-pair cutoffs are applied
before descriptor accumulation. This avoids rebuilding the whole structure or
calling an external calculator for each RMC proposal.

The accept/reject interface is transactional. A second proposal cannot begin
until the previous pending values have either been committed or discarded, and
debug builds verify that the sum of cached atomic energies remains consistent
with the cached total.

## Validation

The implementation is covered at four levels:

- parser and invariant tests for malformed parameters, element mappings,
  coefficient dimensions, and numerical edge cases;
- mathematical tests for bispectrum index counts and order,
  Clebsch--Gordan normalization, Wigner-matrix unitarity, and isolated-atom
  subtraction;
- frozen LAMMPS references for linear Si, explicit chemical InP under both
  self-density conventions, and quadratic Si, each including a displaced
  configuration;
- optional tests against the SNAP potentials distributed with LAMMPS, including
  the linear and quadratic Si Zuo models and the chemical InP model.

The frozen references include descriptors, per-atom energies where needed, and
total energies. Normal tests therefore require no LAMMPS installation while
still guarding compatibility with its numerical convention. The generating
inputs and provenance are kept in `tests/data/snap` so the fixtures can be
reproduced with a LAMMPS build containing ML-SNAP.

An ignored release-mode benchmark compares cached local trials with full model
rebuilds:

```console
cargo test --release --test snap_reference_fixture \
  benchmark_local_trials_against_full_rebuilds -- --ignored --nocapture
```

The source is organized as follows:

- `src/ml_potential/snap_native.rs`: model parsing, validation, descriptor
  assembly, and energy contraction;
- `src/ml_potential/snap_native/math.rs`: angular indices, coupling tables,
  Wigner matrices, and bispectrum contraction;
- `src/ml_potential/snap_native/runtime.rs`: type/cutoff initialization,
  neighbor cells, accepted-state caches, and trial transactions.
