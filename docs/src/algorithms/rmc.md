# RMC Algorithm

The core algorithm follows McGreevy & Pusztai (1988). Each iteration:

1. **Select** a random atom and propose a displacement drawn uniformly from [-max_step, +max_step] in each dimension.
2. **Check constraints**: reject immediately if minimum distance or coordination constraints are violated.
3. **Compute cost**: incrementally update the structure factor and evaluate chi2. If potentials are active, compute the energy change.
4. **Accept/reject** via the Metropolis criterion:
   - If cost decreases: always accept
   - If cost increases by delta: accept with probability exp(-delta / 2T)

## Cost function

```
cost = chi2 + weight * E
```

where chi2 is the total weighted chi-squared across all datasets:

```
chi2 = sum_datasets [ w_d * sum_i ((y_calc(x_i) - y_exp(x_i))^2 / sigma_i^2) ]
```

The sum runs only over data points within the specified fit range. When no potentials are active, cost = chi2.

## Incremental S(Q) updates

Recomputing S(Q) from scratch after each atom move would be O(N^2). Instead, reversesmith uses incremental updates:

1. **Histogram delta**: Compute the RDF histogram contribution of the moved atom at its old and new positions using the cell list. Only bins that change contribute to the S(Q) update.

2. **Partial S(Q) delta**: For each changed bin, compute the contribution to each partial S_ab(Q) via a precomputed sin(Q*r) lookup table:

   ```
   delta_S_ab(Q_k) = (4*pi*rho0*dr / Q_k) * sum_i [r_i * W(r_i) * delta_g_ab(r_i) * sin(Q_k * r_i)]
   ```

3. **Total S(Q)**: Reweight and sum: `S_X(Q) = sum_{ab} w_ab(Q) * S_ab(Q)`

4. **g(r) update** (if fitting g(r)): Apply the delta to the model g(r) via a precomputed FT matrix.

This reduces the per-move cost from O(N^2 * N_Q) to O(N_neighbors * N_bins * N_Q).

## Cell list

A linked-cell list partitions atoms into spatial cells of size >= rdf_cutoff. Neighbour searches only visit 27 adjacent cells, reducing iteration from O(N) to O(N_neighbors). The cell list is updated on accepted moves.

The same cell list is used for:
- RDF histogram computation
- Minimum distance constraint checking
- Coordination constraint checking
- Pair potential energy computation

## Metropolis criterion

At temperature T, the acceptance probability for a move with delta_cost > 0 is:

```
P = exp(-delta_cost / (2 * T))
```

The factor of 2 in the denominator follows the RMC convention (McGreevy & Pusztai 1988). At T=1 (standard RMC), this samples configurations consistent with the data within experimental uncertainty.

For optimization, use T < 1 (simulated annealing) or a cooling schedule. See [RMC Parameters](../config/rmc.md) for details.
