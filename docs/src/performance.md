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

## Tips

- Use `cargo build --release` -- debug builds are ~10x slower.
- Set `rdf_cutoff` and potential `cutoff` no larger than needed.
- Run multiple independent replicas with different seeds (trivially parallel) and pick the best result.
- The `--compute-sq-only` mode is useful for quickly checking S(Q) quality before committing to a long refinement run.
