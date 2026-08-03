# Performance Notes

## Benchmark

10,000-atom CaSiO3 system (box 51.5 A), 550 RDF bins, 500 Q points, 6 pair types, fitting S(Q) + g(r) simultaneously:

- ~390 us per move
- ~650k moves in ~250s (release build, Apple Silicon)

## Per-move cost breakdown

| Component | Cost | % |
|-----------|------|---|
| Delta-S(Q) computation | ~275 us | 71% |
| g(r) inverse FT update | ~50 us | 13% |
| RDF histogram (via cell list) | ~13 us | 3% |
| Constraint checking | ~6 us | 2% |
| Chi2, bookkeeping, etc. | ~44 us | 11% |

The Delta-S(Q) computation dominates because nearly all histogram bins change per move (each atom has ~4000 neighbours spread across most bins). The inner loop is a SAXPY (scalar * vector + vector accumulate) of length N_Q, repeated ~N_bins * N_pairs times. LLVM auto-vectorizes this effectively.

## Pair potential overhead

When potentials are active, the energy computation adds O(N_neighbors)
distance calculations and table lookups per move. Pure EPSR mode stores the
candidate indices in a directed cutoff-plus-skin Verlet list and rebuilds the
list before accumulated displacements could invalidate its skin. The potential
evaluation itself is a single linear interpolation for neighbors inside the
physical cutoff.

For the 5,000-atom EPSR26 LiquidGa50C example, a 12 A physical cutoff plus a
3 A skin gives about 738 candidates per atom. The previous 3x3x3 cell stencil
visited essentially all 5,000 atoms per trial. Replacing it with the Verlet
list reduced the ten-seed, one-million-move median from 44.64 s to 16.71 s;
native EPSR26 required 29.09 s under the same one-thread protocol. Acceptance,
S(Q), and g(r) distributions were numerically unchanged.

A separate ten-seed convergence test starts both programs from identical
independently generated hard-sphere liquids. Using a common external residual,
rsmith reaches a sustained 3.0% fit in a median 2.02 s versus 22.40 s for
native EPSR26 (11.11x faster), and 2.5% in 3.72 s versus 22.58 s (6.07x).
All ten rsmith runs reach 2.0% within 7.11 s and 1.5% within 13.88 s; no native
run reaches either target within one million moves. See the
`benchmarks/epsr-ga` record for the frozen checkpoint rule, censoring,
structural diagnostics, and scientific limitations.

## Delayed acceptance

Hybrid RMC can optionally apply the cheap experimental-data Metropolis test
before evaluating the configured energy model. Set
`[rmc] delayed_acceptance = true`; after the initial weight-calibration window,
data-stage rejections skip the potential call entirely. This applies uniformly
to pair potentials, SNAP, GAP, and MACE, although the largest wall-time savings
are expected for GAP and MACE.

If the fraction of proposals passing the data stage is `p`, the approximate
per-attempt cost is `data_cost + p * potential_cost`. For example, a 30% pass
rate would reduce a 21 ms GAP-dominated attempt to roughly 7 ms, or a 56 ms
MACE-dominated attempt to roughly 17 ms, before accounting for changes in the
overall acceptance rate. Delayed acceptance preserves the target distribution
but can mix more slowly when the experimental and energy changes oppose each
other, so benchmark accepted or effectively independent structures per second.

## MACE energy-delta modes

The MACE/Python backend offers three strategies for calculating the energy
change of a single-atom proposal. `full` evaluates the entire system before and
after every move and therefore grows with atom count. `local` evaluates one
dense finite-range cluster and reaches a constant cost once the simulation box
is larger than that cluster. `incremental` caches accepted hidden features and
recomputes only the changed subgraph at each message-passing layer, giving the
lowest constant local cost for supported models.

Median times from the Apple M4 Pro benchmark using the small MACE-MP model,
eight CPU threads, and rejected moves are:

| Atoms | Full float64 | Local float64 | Incremental float64 | Incremental float32 |
|------:|-------------:|--------------:|--------------------:|--------------------:|
| 216 | 123.5 ms | 471.7 ms | 54.8 ms | 36.0 ms |
| 1,000 | 501.5 ms | 577.4 ms | 80.5 ms | 52.8 ms |
| 1,728 | 865.9 ms | 574.8 ms | 80.2 ms | 52.3 ms |
| 4,096 | 2.060 s | 572.8 ms | 81.2 ms | 52.5 ms |
| 8,000 | 4.308 s | 576.9 ms | 81.7 ms | 53.8 ms |

The 216-atom local cluster contains many explicit periodic images, making dense
`local` slower than `full`. This is why the correct mode should be chosen based
on both checkpoint compatibility and the actual system size. See
[ML Potentials](./config/ml-potential.md#choosing-the-mace-delta-mode) for the
compatibility rules and selection workflow.

For incremental float32 at 4,096 atoms, eight threads gave the best measured
time (52.1 ms), compared with 105.3 ms for one thread. Ten and fourteen threads
were slower. Thread scaling is therefore workload- and machine-specific rather
than proportional to core count.

## GAP/QUIP performance

GAP supports `delta = "full"` and `delta = "local"`. Full mode evaluates both
complete periodic systems. Local mode caches accepted per-atom energies,
constructs a bounded explicit-image cluster, and evaluates only affected
central atoms once per trial. The adapter requests local energies without
calculating forces or virials.

The following Apple M4 Pro measurements use the published general-purpose Si
PRX GAP (`GAP_2017_6_17_60_4_3_56_165`): a 5 A SOAP descriptor with 9,000
sparse points plus its repulsive core potential. Full times are medians of
three rejected trials after one warm-up; local times are medians of seven after
two warm-ups. Both use eight OpenMP and OpenBLAS threads. The structures and
displacement sequence match the SNAP and MACE benchmarks.

| Atoms | GAP full | GAP local | Local speedup | SNAP native | MACE incremental float32 |
|------:|---------:|----------:|--------------:|------------:|--------------------------:|
| 216 | 252.6 ms | 21.17 ms | 11.9x | 2.415 ms | 36.0 ms |
| 1,000 | 1.143 s | 21.07 ms | 54.3x | 2.377 ms | 52.8 ms |
| 1,728 | 1.980 s | 21.15 ms | 93.6x | 2.387 ms | 52.3 ms |
| 4,096 | 4.730 s | 21.32 ms | 221.8x | 2.399 ms | 52.5 ms |
| 8,000 | 9.527 s | 21.81 ms | 436.9x | 2.407 ms | 53.8 ms |

Every local case evaluated 29 central atoms in a 413-atom/image cluster. The
small 3% change from 216 to 8,000 atoms shows that steady-state trial cost is
effectively independent of total atom count. Local mode pays a one-time full
evaluation to initialize its accepted-energy cache; this took about 4.7 s for
8,000 atoms in the same process, in addition to model loading.

These potentials are not equivalent in complexity or accuracy. In particular,
the benchmark GAP has 9,000 SOAP sparse points, while the MACE result uses the
small MACE-MP test checkpoint and the SNAP result uses the linear
`Si_Zuo_JPCA2020` model. The table measures current rsmith backend cost on a
common structure, not scientific model quality.

Local CPU scaling at 1,000 atoms was useful but sublinear: 39.38, 27.52, 25.15,
21.25, and 21.00 ms at 1, 2, 4, 8, and 14 threads respectively. Fourteen
threads were 1.88x faster than one, with nearly all practical gain reached by
eight threads. Independent RMC replicas remain the more efficient way to use
additional cores after that point.

## Scaling

- **Atom count**: at fixed density and cutoff, neighbor-based trial cost is
  approximately independent of total atom count; a fixed number of moves per
  atom therefore gives approximately O(N) total wall time.
- **Q points**: O(N_Q) per move. Reducing Q points from 500 to 250 halves the dominant cost.
- **RDF bins**: O(N_bins) per move. Fewer bins (e.g., 275 at 0.04 A resolution) can significantly reduce cost.
- **Potential cutoff**: larger cutoffs increase the Verlet-list population
  approximately with the cutoff-plus-skin volume.

## RDF cutoff and performance

The `rdf_cutoff` has the largest impact on per-move cost after `nq`. Here is an approximate comparison for 10,000 atoms (box ~51 A, 500 Q points):

| `rdf_cutoff` | `rdf_nbins` | Neighbors/atom | Per-move time | Note |
|---|---|---|---|---|
| 11 A | 550 | ~4000 | ~390 us | Fast, FSDP may be inaccurate |
| 15 A | 750 | ~10000 | ~700 us | Good for most applications |
| 20 A | 1000 | ~24000 | ~1.5 ms | Accurate FSDP |
| 25 A | 1250 | ~47000 | ~3 ms | Matches LAMMPS `compute rdf` quality |

The dominant scaling factor is the number of neighbors, which grows as `cutoff^3`. Each neighbor contributes to multiple RDF bins and thus to the incremental S(Q) update.

See [S(Q) Computation](./config/sq.md) for guidance on choosing `rdf_cutoff`.

## Optimization history

This section documents optimization attempts and their outcomes, to avoid re-investigating
dead ends and to record what actually helped.

### Implemented: g(r) fit-range-only FT matrix

**Result: ~10% faster for configs with g(r) data.**

The inverse Fourier transform that converts total S(Q) to model g(r) was computed for all
r-points in the experimental dataset (e.g., 1001 points for a 0–10 Å range). However, chi2
is only evaluated over the fit range (e.g., 0–7 Å = 700 points). By building the FT matrix
and computing the g(r) update only for fit-range points, we eliminate ~30% of the FT work
without any approximation.

Benchmark (1560-atom CaSiO3, 50K moves, S(Q) + g(r) + Pedone potential):
- Before: 48.5s
- After: 43.3s (−10.7%)
- S(Q)-only configs are unaffected

### Attempted and reverted: combined sin table

**Result: 7–8% regression in S(Q) path.**

Folded `rw[i] * prefactor_sq * inv_q[k]` into `sin_table` to eliminate the separate
prefactor pass after the bin loop. Despite doing fewer total operations, the combined table
was consistently slower. The original two-phase approach (accumulate with raw sin values,
then apply prefactors in a separate pass) appears to optimize better under LLVM — likely
because the separate passes are individually simpler and more amenable to auto-vectorization.

### Implemented: Chirp-Z Transform (CZT) via Bluestein's algorithm

**Result: ~25% faster for S(Q)-only, ~16% faster for full configs.**

Replaced the O(nbins × nq) sin table lookup with an FFT-based Chirp-Z Transform. The S(Q)
delta computation is a discrete sine transform with non-conjugate Q and r grids, which maps
naturally to Bluestein's algorithm: pre-chirp the input, zero-pad to length L (next power of
2 ≥ nbins + nq − 1 = 4096), FFT, pointwise multiply with a precomputed kernel, IFFT,
post-chirp and extract imaginary parts.

The key advantage is cache locality: the entire CZT working set (~128 KB) fits in L1 cache,
while the sin table (~5 MB in f32) spills to L2/DRAM. Despite the FFT overhead, the
cache-friendly access pattern wins decisively.

Benchmark (1560-atom CaSiO3, 50K moves):
- S(Q)-only: 25.6s → 19.2s (−25%)
- Full config (S(Q) + g(r) + Pedone potential): 39.4s → 33.1s (−16%)
- Chi2 identical (f64 throughout, no precision loss)

This also eliminates the 5 MB sin table allocation.

### Superseded: f32 sin table (now replaced by CZT)

**Result was: ~9% faster than f64 sin table.**

Stored sin(Q_k * r_i) as f32 instead of f64 to halve memory bandwidth in the S(Q) delta
inner loop. This reduced the sin table from 10 MB to 5 MB and gave a measurable speedup on
the memory-bandwidth-limited SAXPY loop. However, this approach was superseded by the CZT,
which eliminates the sin table entirely.

### Attempted and reverted: Goertzel/Clenshaw recurrence (no sin table)

**Result: 5.5× slower.**

Replaced the 10 MB sin lookup table with on-the-fly computation using the Goertzel
recurrence `sin(Q_k * r) = 2·cos(dQ·r)·sin(Q_{k-1}·r) − sin(Q_{k-2}·r)`. This eliminates
all memory bandwidth for the sin table (which exceeds L2 cache). However, the recurrence has
a sequential data dependency (each value depends on the two previous), which prevents SIMD
vectorization of the inner loop. The table-based approach auto-vectorizes to 4-wide NEON on
Apple Silicon, making it ~5× faster despite the cache misses. Hardware prefetching of the
sequential table access pattern further reduces the memory penalty.

### Investigated but not implemented: affected-pairs-only S(Q) delta

**Estimated savings: <0.5%.**

When moving an atom of type t, only the 3 pair channels involving t (out of 6 for a 3-type
system) have nonzero histogram deltas. Skipping the other 3 channels avoids scanning 2500
zero-valued bins and one prefactor pass of 500 multiplies per skipped channel. However, the
existing `if dh == 0.0 { continue; }` check already short-circuits unaffected bins at
negligible cost. Total savings: ~9000 trivial operations per move = ~0.1s over 50K moves.

### Investigated but not implemented: sparse bin tracking

**Estimated savings: <0.5%.**

Collect indices of nonzero histogram bins during the histogram computation to avoid scanning
all 2500 bins in the S(Q) delta loop. The overhead of collecting and storing sparse indices
offsets the minor savings from skipping zero bins, especially since the `dh == 0.0` branch
is well-predicted by the CPU.

### Investigated but not implemented: CellList precomputed neighbor table

**Measured cost: <1% of total runtime.**

Profiling confirmed that `CellList::neighbor_cells()` accounts for fewer than 0.02% of
samples. Precomputing a neighbor table adds complexity with no measurable benefit.

### Implemented: Fine cell list for constraint checking

**Result: ~51% faster for S(Q)-only, ~35% faster for full configs.**

The constraint checking functions (`check_min_distances_cell`, `check_coordination_cell`) used
the same cell list as the RDF histogram, built with `rdf_cutoff` (e.g., 21 Å). For a ~44 Å
box, this gives `nc=3` per dimension (cell size 14.7 Å), so `neighbor_cells()` returns all
27 cells = the entire box. Every constraint check iterated all 10,000 atoms, even though
min_distance constraints only need ~2 Å and coordination constraints need ~4.4 Å radius.

The fix: build a second, finer cell list with `cutoff = max(max_min_distance, 2 * max_coordination_cutoff)`.
For typical CaSiO3 constraints (max_min_dist=2.0 Å, coord_cutoff=2.2 Å), this gives
cutoff=4.4 Å, nc=10, 1000 cells. Each constraint check now iterates ~27/1000 of the atoms
(~270 atoms in 27 neighbor cells) instead of all 10,000.

Benchmark (10,000-atom CaSiO3, 50K moves):
- S(Q)-only: 16.9s → 8.3s (−51%)
- Full config (S(Q) + g(r) + Pedone potential): 30.3s → 19.6s (−35%)

### Bottleneck analysis

With the CZT, the S(Q) delta computation is now **compute-bound** (FFTs in L1 cache) rather
than memory-bandwidth-limited. The per-move cost scales as O(L log L) where L is the FFT
length (next power of 2 ≥ nbins + nq). For typical parameters (nbins=2500, nq=500), L=4096.
The dominant remaining costs are the histogram computation and g(r) inverse FT update.

## Tips

- Use `cargo build --release` -- debug builds are ~10x slower.
- Set `rdf_cutoff` and potential `cutoff` no larger than needed. For quick convergence tests, start with 11--15 A; for publication-quality fits, use 20+ A.
- Run multiple independent replicas with different seeds (trivially parallel) and pick the best result.
