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

When potentials are active, the energy computation adds O(N_neighbors) distance calculations and table lookups per move. For a 15 A cutoff, this is ~8000 neighbours, adding roughly 20-30 us per move (~5-8% overhead). The potential evaluation itself is a single linear interpolation per neighbour.

## Scaling

- **Atom count**: O(N) per move (via cell list). Wall time scales linearly with system size for fixed number of moves per atom.
- **Q points**: O(N_Q) per move. Reducing Q points from 500 to 250 halves the dominant cost.
- **RDF bins**: O(N_bins) per move. Fewer bins (e.g., 275 at 0.04 A resolution) can significantly reduce cost.
- **Potential cutoff**: Larger cutoff = more neighbours = more overhead for both histograms and potentials.

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

### Bottleneck analysis

With the CZT, the S(Q) delta computation is now **compute-bound** (FFTs in L1 cache) rather
than memory-bandwidth-limited. The per-move cost scales as O(L log L) where L is the FFT
length (next power of 2 ≥ nbins + nq). For typical parameters (nbins=2500, nq=500), L=4096.
The dominant remaining costs are the histogram computation and g(r) inverse FT update.

## Tips

- Use `cargo build --release` -- debug builds are ~10x slower.
- Set `rdf_cutoff` and potential `cutoff` no larger than needed. For quick convergence tests, start with 11--15 A; for publication-quality fits, use 20+ A.
- Use `--compute-sq-only` to check S(Q) quality at different `rdf_cutoff` values before committing to a long refinement.
- Run multiple independent replicas with different seeds (trivially parallel) and pick the best result.
