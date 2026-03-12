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

The S(Q) delta computation (the SAXPY inner loop over Q-points for each nonzero RDF bin) is
**memory-bandwidth-limited**. The sin table (nbins × nq × 8 bytes) typically exceeds L2
cache (e.g., 2500 × 500 × 8 = 10 MB). Each move accesses ~6 MB of sin table data (one 4 KB
row per nonzero bin × ~1500 bins for large-cutoff systems). Hardware prefetching mitigates
this for sequential access, but the fundamental constraint is DRAM bandwidth, not compute.

This means the dominant per-move cost scales with `rdf_nbins × nq × sizeof(f64)` — reducing
either `rdf_nbins` (via smaller `rdf_cutoff`) or `nq` is the most effective way to speed up
the refinement.

## Tips

- Use `cargo build --release` -- debug builds are ~10x slower.
- Set `rdf_cutoff` and potential `cutoff` no larger than needed. For quick convergence tests, start with 11--15 A; for publication-quality fits, use 20+ A.
- Use `--compute-sq-only` to check S(Q) quality at different `rdf_cutoff` values before committing to a long refinement.
- Run multiple independent replicas with different seeds (trivially parallel) and pick the best result.
