# Pair Potentials (Hybrid RMC)

## Why use potentials?

Standard RMC is underdetermined -- many atomic configurations can produce the same scattering pattern. Some of these configurations are physically unreasonable: distorted coordination polyhedra, unphysical bond lengths, or energetically unfavorable arrangements. Hard constraints (minimum distances, coordination bounds) help, but they are binary and cannot express the smooth energy landscape that governs real materials.

Adding pair potentials to the acceptance criterion biases the refinement toward energetically favorable configurations. This is the simplest form of Empirical Potential Structure Refinement (EPSR): the potential is fixed (not iteratively refined), and its contribution is weighted against chi2.

## How it works

### Acceptance criterion

Without potentials (standard RMC):
```
accept if delta_chi2 < 0, else with prob exp(-delta_chi2 / 2T)
```

With potentials (hybrid RMC):
```
delta_cost = delta_chi2 + weight * delta_E
accept if delta_cost < 0, else with prob exp(-delta_cost / 2T)
```

The `weight` parameter controls the balance:
- `weight = 0` -- pure RMC, potentials computed but don't affect acceptance
- Small `weight` -- potentials act as a gentle bias, chi2 dominates
- Large `weight` -- potential dominates, structure stays near energy minimum but may not fit the data well

### Energy computation

The energy change for a single atom move is computed efficiently:
```
delta_E = E_new(atom_i) - E_old(atom_i)
```
where `E(atom_i) = sum_j V(r_ij)` sums over all neighbours j within the potential cutoff. This uses the same cell list as the RDF computation, so the cost is O(N_neighbours) per move -- negligible compared to the S(Q) update.

### Best-structure tracking

The best structure is tracked on **chi2 alone** (not the combined cost). The potential guides the path, but the goal is still the best fit to the experimental data.

## Choosing the weight

This is the most important practical decision. There is no universal formula -- the right weight depends on the relative magnitudes of delta_chi2 and delta_E for your system.

### Automatic calibration

Reversesmith prints a calibration diagnostic after the first 1000 trial moves:

```
Calibration (1000 moves): avg |delta_chi2| = 0.0023, avg |delta_E| = 0.34 eV
  Current weight = 0.001000, suggested weight for equal balance = 0.006765
  Ratio current/suggested = 0.1478
```

The "suggested weight" equalizes the average contributions of chi2 and energy to the cost function. This is a starting point, not necessarily the optimal value.

### Practical guidelines

1. **Start with the suggested weight from calibration.** This gives roughly equal influence to data and potential.

2. **For an MD-derived starting structure**, the energy is already near its minimum. The goal is to keep it there while fitting the data. Start with a weight that makes the energy contribution ~0.5--1x the chi2 contribution.

3. **If energy increases significantly during refinement**, the weight is too low. The RMC is distorting the structure to fit scattering noise. Increase weight by 2--5x.

4. **If chi2 doesn't decrease**, the weight is too high. The potential is preventing the structure from adjusting to fit the data. Decrease weight by 2--5x.

5. **Compare final chi2 with and without potentials.** A good weight gives chi2 perhaps 10--30% higher than pure RMC, but with much better local structure (check with `--analyze`).

6. **Two to three trial runs** with different weights usually bracket the right range.

### Interpreting the status output

With potentials active, the status line shows all components:

```
Move 10000/500000: cost = -82.3 (chi2: 3.2 [sq: 2.1, gr: 1.1], w*E: -85.5 [E: -85500.00]), ...
```

- **cost** -- the quantity being minimized (chi2 + weight * E)
- **chi2** -- data fit only (sq + gr contributions)
- **w*E** -- weighted energy contribution (what the acceptance criterion sees)
- **E** -- raw potential energy in eV

The energy of an MD-equilibrated structure is large and negative (e.g., -85000 eV for 10000 atoms). During refinement:
- **E staying flat** -- good, the potential is acting as a guardrail
- **E increasing (less negative)** -- the structure is being distorted, consider increasing weight
- **E decreasing significantly** -- unusual, may indicate the structure is relaxing toward a lower-energy basin

## Potential types

### Which potential to use

Use the **same potential that generated your starting structure**. If you ran MD with Pedone + Coulomb in LAMMPS, use Pedone + Coulomb in rsmith. This ensures the energy surface the RMC explores is consistent with the initial equilibrium structure.

If you don't have analytical parameters, you can export tabulated potentials from LAMMPS using `pair_write` and use the `[[potential.tabulated]]` format.

### Pedone + Coulomb (recommended for oxide glasses)

The Pedone potential (Pedone et al., J. Phys. Chem. B, 2006, 110, 11780) is widely used for silicate and oxide glasses. It combines a Morse-type short-range interaction with r^-12 repulsion:

```
V_pedone(r) = D0 * [1 - exp(-alpha * (r - r0))]^2 - D0 + C0 / r^12
```

The Coulomb part uses Damped Shifted Force (DSF) electrostatics:

```
V_coul(r) = K * qi * qj * [erfc(alpha*r)/r - erfc(alpha*rc)/rc
            + (erfc(alpha*rc)/rc^2 + 2*alpha/sqrt(pi) * exp(-alpha^2*rc^2)/rc) * (r - rc)]
```

where K = 14.3997 eV A / e^2.

The Pedone and Coulomb contributions are summed for each pair. This exactly matches a LAMMPS `hybrid/overlay pedone coul/dsf` setup.

### Buckingham

The Buckingham potential is an alternative short-range form:

```
V(r) = A * exp(-r/rho) - C / r^6
```

Used by some older force fields (e.g., van Beest, Kramer, van Santen 1990 for silicates).

### Tabulated

For potentials not covered by the built-in forms, provide a two-column file with r (A) and V(r) (eV). The potential is interpolated onto a fine grid (dr = 0.001 A) and shifted so V(cutoff) = 0.

Generate tables from LAMMPS:
```
pair_write 1 3 5000 r 0.5 15.0 CaO_potential.dat Ca-O
```

## Cutoff considerations

- All potentials are shifted so V(r_cutoff) = 0, avoiding energy discontinuities.
- The potential cutoff must not exceed `rdf_cutoff` in `[sq]` (the same cell list is used for both).
- For Pedone + Coulomb, a cutoff of 15 A (matching the LAMMPS simulation) is typical. This requires `rdf_cutoff >= 15.0`.
- Larger cutoffs are more accurate but increase computation cost and require a larger cell list.
- Tabulation uses dr = 0.001 A (1000 points per A), which is fine enough for smooth interpolation.
- At very short range (r < 0.5 A for analytical, r < 0.3 A for Coulomb), potentials are capped to avoid numerical divergence. This region is normally excluded by minimum distance constraints.

## Example workflow

1. Run MD with Pedone + Coulomb in LAMMPS to generate starting structure
2. Set up rsmith with the same potential parameters
3. Run with `weight = 0` first to establish baseline chi2
4. Check calibration output for suggested weight
5. Run with suggested weight, monitor chi2 and E
6. Adjust weight if needed (2--3 iterations)
7. Compare structures with `--analyze` to verify improved local order

## References

- Pedone, A., Malavasi, G., Menziani, M.C., Cormack, A.N. & Segre, U. (2006). A New Self-Consistent Empirical Interatomic Potential Model for Oxides, Silicates, and Silica-Based Glasses. *J. Phys. Chem. B*, 110, 11780--11795.
- Fennell, C.J. & Gezelter, J.D. (2006). Is the Ewald summation still necessary? Pairwise alternatives to the accepted standard for long-range electrostatics. *J. Chem. Phys.*, 124, 234104.
- McGreevy, R.L. (2001). Reverse Monte Carlo modelling. *J. Phys.: Condens. Matter*, 13, R877--R913.
