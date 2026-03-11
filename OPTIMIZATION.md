# Performance Optimization Notes

## Summary

Optimized the RMC hot loop from ~287s to ~251s (12.5% speedup) for a 10,000-atom
CaSiO3 system with 650k moves. Results are bit-identical to the original code.

## System under test

- 10,000 atoms (Ca, Si, O), box 51.5³ Å
- 550 RDF bins, 500 Q points, 6 atom-pair types
- Fitting against X-ray S(Q) and g(r) simultaneously
- Coordination constraint: Si-O [3,6] within 2.2 Å
- 6 minimum-distance constraints (1.2–2.3 Å)

## Bottleneck analysis

Per-move cost breakdown (~388 μs/move total):

| Component | Est. cost | % of total |
|---|---|---|
| **ΔS(Q) computation** (6 pairs × 550 bins × 500 Q) | ~275 μs | ~71% |
| **g(r) inverse FT** (1001 × 500 matvec) | ~50 μs | ~13% |
| **Histogram** (2 × ~4200 neighbor distances via cell list) | ~13 μs | ~3% |
| **Constraint checking** (~4200 neighbors via cell list) | ~6 μs | ~2% |
| **Chi² evaluation**, bookkeeping, etc. | ~44 μs | ~11% |

The ΔS(Q) computation dominates because nearly all 550 bins change per move
(each atom has ~4000 neighbors spread across most bins). The inner loop is a
SAXPY (scalar × vector + vector accumulate) of length 500, repeated ~550 × 6 =
3300 times. LLVM auto-vectorizes this effectively.

## Optimizations applied

### 1. Precomputed constraint lookup table (high impact)

**Problem:** `check_min_distances()` called `format!("{}-{}", sa, sb)` to build
HashMap keys on every atom pair — ~10,000 string allocations per move × 650k
moves = billions of heap allocations.

**Fix:** `PrecomputedConstraints` stores a flat `min_dist_sq[ti * n_types + tj]`
array built once before the loop. Constraint checks use direct array indexing
instead of string formatting + HashMap lookup. Same for coordination constraints
(type IDs instead of string parsing).

Files: `constraints.rs` (added `PrecomputedConstraints`, `check_min_distances_fast`,
`check_coordination_fast`), `rmc.rs` (use new functions).

### 2. Reused histogram buffers (moderate impact)

**Problem:** `atom_histogram_st()` allocated a new `Vec<f64>` (3300 elements)
on every call — 2 calls per move = 1.3M allocations.

**Fix:** Pre-allocate `old_hist_buf` and `new_hist_buf` once, pass as `&mut [f64]`,
clear with `.fill(0.0)` each move.

### 3. Cell list for spatial lookups (moderate impact)

**Problem:** Constraint checking and histogram computation iterated over ALL
atoms O(N) per move.

**Fix:** Added `cells.rs` implementing a linked-cell list with periodic boundary
conditions. Neighbor lookups only visit 27 adjacent cells. For the current system
(box/cutoff = 51.5/11 = 4.7 → 4³ = 64 cells), this reduces neighbor iteration
from 10,000 to ~4,200 atoms (2.4× reduction).

The cell list is updated on accepted moves via `cell_list.move_atom()`.
`neighbor_cells()` returns a fixed `[usize; 27]` array to avoid per-call
heap allocation.

### 4. Precomputed g(r) delta vector (minor impact)

**Problem:** The g(r) inverse FT inner loop computed `new_total_sq[k] - total_sq[k]`
redundantly for each r point.

**Fix:** Precompute `delta_sq_buf[k]` once per move, enabling better
auto-vectorization of the dot product loop.

### 5. Precomputed annealing log ratio (minor impact)

**Problem:** `powf()` called every move for annealing temperature.

**Fix:** Precompute `ln(T_end/T_start)` once, use `exp(frac * log_ratio)` per move.

## Approaches tried but not adopted

### FFT-based DST for ΔS(Q)

The ΔS(Q) computation IS a Discrete Sine Transform. Replacing the brute-force
O(nbins × nq) sum with FFT-based DST should theoretically help.

**Why it didn't work:** For 550 bins, the minimum useful FFT size is 2×550 =
1100 → N_fft = 2048 (power of 2), giving FFT length 4096. The complex FFT
cost (~150k real ops) plus interpolation overhead is comparable to the brute-force
(~275k SAXPY ops that auto-vectorize very well). With 8× padding (N_fft = 8192,
FFT length 16384), the FFT is significantly slower.

Additionally, the FFT approach introduces interpolation error when mapping from
the DST Q-grid to the target Q-grid. With only 2× padding, this error accumulated
over 500k+ moves, degrading final chi² from 1666 to 1706.

### Multi-pair simultaneous accumulation

Process all 6 pairs for each bin simultaneously, reading each sin_table row once
instead of 6 times (reducing memory bandwidth from 13.2 MB to 2.2 MB).

**Why it didn't work:** The conditional `if pair_contrib_buf[p] != 0.0` in the
innermost loop prevented LLVM auto-vectorization. The resulting scalar code was
~1.7× slower despite better cache behavior.

### Rayon parallelization of ΔS(Q) pairs

Dispatch the 6 independent pair computations to rayon's thread pool.

**Why it didn't work:** Rayon's work-stealing overhead (~5–10 μs per dispatch)
is too large relative to per-pair cost (~45 μs). Over 650k moves, the thread
synchronization overhead dominated. User CPU time increased from 251s to 440s.

### Rayon parallelization of g(r) matvec

Parallelize the 1001 independent dot products across threads.

**Not attempted** after the ΔS(Q) parallelization showed that rayon overhead
is prohibitive for sub-millisecond tasks dispatched hundreds of thousands of times.

## Remaining optimization opportunities

### Reduce problem size (if scientifically acceptable)

- **Fewer RDF bins:** 275 bins (0.04 Å resolution) instead of 550 → ~2× faster
  ΔS(Q). Check if S(Q) quality degrades.
- **Fewer Q points:** 250 instead of 500 → ~2× faster ΔS(Q) and g(r) matvec.
- **Smaller RDF cutoff:** Reduces neighbor count and makes cell list more effective.

### Multiple independent replicas

Run N instances with different seeds on different cores (trivially parallel at
the process level). Pick the result with the lowest chi². This scales perfectly
with core count and doesn't change the algorithm.

### SIMD intrinsics for the SAXPY inner loop

The compiler auto-vectorizes reasonably well, but hand-tuned SIMD (NEON on
Apple Silicon, AVX2 on x86) could potentially squeeze out another 20–30% by
using fused multiply-add and optimal register allocation.

### Sparse histogram tracking

Track which bins actually change (indices of non-zero `dh`) to skip the
`if dh == 0.0 { continue; }` branch entirely. Whether this helps depends on
what fraction of bins change per move (~95% for large systems, so marginal).
