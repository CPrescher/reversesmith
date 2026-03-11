# Pair Potentials Theory

## Potential energy in RMC

The total pair potential energy of a configuration is:

```
U = sum_{i<j} V_ab(r_ij)
```

where V_ab is the pair potential for species types a and b, and the sum runs over all pairs within the cutoff distance.

Moving a single atom i changes the energy by:

```
delta_U = sum_j [V(r_ij_new) - V(r_ij_old)]
```

where the sum runs over all neighbours j of atom i. This is O(N_neighbors), not O(N^2), using the cell list.

## Tabulation

All potentials (analytical and tabulated) are stored on a uniform grid with spacing dr = 0.001 A. Lookup uses linear interpolation between grid points. This ensures:
- Consistent, fast evaluation regardless of functional form
- No branching in the inner loop (same code path for all potential types)
- Smooth interpolation with negligible discretization error

## Cutoff and shifting

All potentials are shifted so V(r_cutoff) = 0:

```
V_shifted(r) = V(r) - V(r_cutoff)
```

This avoids energy discontinuities when atoms cross the cutoff boundary. The shift is a constant offset that does not affect forces or relative energy differences between configurations. It does change the absolute energy, but since only delta_U matters for acceptance, this is immaterial.

Short-range capping:
- Analytical potentials (Buckingham, Pedone): capped at r = 0.5 A
- Coulomb DSF: capped at r = 0.3 A
- Tabulated: extrapolated from the first data point

These regions are normally excluded by minimum distance constraints.

## Buckingham potential

```
V(r) = A * exp(-r / rho) - C / r^6
```

The first term is Born-Mayer repulsion, the second is van der Waals attraction. Used by some classical force fields for oxides (e.g., van Beest-Kramer-van Santen / BKS for silica).

## Pedone potential

```
V(r) = D0 * [1 - exp(-alpha * (r - r0))]^2 - D0 + C0 / r^12
```

Combines a Morse potential (well depth D0, width alpha, equilibrium distance r0) with a short-range r^-12 repulsion (coefficient C0). Developed for oxide glasses by Pedone et al. (2006).

The Morse part provides:
- A well at r = r0 with depth D0
- Asymmetric shape: steeper repulsion at r < r0, gradual attraction at r > r0

The r^-12 term prevents atomic overlap at very short distances.

## Coulomb DSF

The Damped Shifted Force (DSF) method provides a pairwise-additive approximation to Coulomb interactions that goes smoothly to zero at the cutoff:

```
V(r) = K * qi * qj * [erfc(alpha * r) / r - erfc(alpha * rc) / rc
       + (erfc(alpha * rc) / rc^2 + 2*alpha/sqrt(pi) * exp(-alpha^2 * rc^2) / rc) * (r - rc)]
```

where:
- K = 14.3997 eV A / e^2 (Coulomb constant in metal units)
- alpha = damping parameter (1/A)
- rc = cutoff distance
- erfc = complementary error function

The DSF form has three components:
1. `erfc(alpha*r)/r` -- damped Coulomb interaction (decays faster than 1/r)
2. `- erfc(alpha*rc)/rc` -- energy shift so V(rc) approaches 0
3. Linear correction term -- ensures V(rc) = 0 exactly and the force is continuous

This matches the LAMMPS `pair_style coul/dsf` implementation. The alpha parameter controls the damping rate -- larger alpha means faster decay, better convergence at the cutoff, but potentially missing long-range contributions.

## Combination of potentials

For each atom pair type, the total effective potential is the sum of all applicable analytical forms:

```
V_total(r) = V_short_range(r) + V_coulomb(r)
```

For example, the Ca-O interaction with Pedone + Coulomb DSF:

```
V_Ca-O(r) = [D0*(1-exp(-alpha*(r-r0)))^2 - D0 + C0/r^12]  +  [K * q_Ca * q_O * DSF(r)]
```

These are summed at the tabulation stage, producing a single lookup table per pair. During RMC, only one table lookup per pair per neighbour is needed.

## References

- Pedone, A. et al. (2006). *J. Phys. Chem. B*, 110, 11780--11795.
- Fennell, C.J. & Gezelter, J.D. (2006). *J. Chem. Phys.*, 124, 234104.
- Wolf, D. et al. (1999). *J. Chem. Phys.*, 110, 8254--8282.
- van Beest, B.W.H., Kramer, G.J. & van Santen, R.A. (1990). *Phys. Rev. Lett.*, 64, 1955.
