# `[potential]` -- Pair Potentials (Hybrid RMC)

Adding pair potentials biases the RMC refinement toward energetically favorable configurations. The acceptance criterion becomes:

```
delta_cost = delta_chi2 + weight * delta_E
```

where `delta_E` is the change in pair potential energy for the proposed move. See the [Pair Potentials](../potentials.md) chapter for theory and best practices.

## Top-level settings

```toml
[potential]
weight = 0.001               # Weight of energy term relative to chi2
cutoff = 8.0                 # Pair potential cutoff (A), default: rdf_cutoff
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `weight` | Float | 0.001 | Scales the energy contribution in the cost function |
| `cutoff` | Float | rdf_cutoff | Global cutoff for all pair potentials (A) |

The cutoff must not exceed `rdf_cutoff` in `[sq]` (the same cell list is used).

## Potential types

Multiple potential types can be combined. For each pair, the contributions are **summed** (e.g., Pedone short-range + Coulomb). If a tabulated potential is also defined for the same pair, it **replaces** (not adds to) the analytical forms.

### Buckingham

V(r) = A exp(-r/rho) - C/r^6

```toml
[[potential.buckingham]]
pair = "Si-O"
A = 18003.7572
rho = 0.205205
C = 133.5381
```

| Field | Type | Description |
|-------|------|-------------|
| `pair` | String | Species pair, e.g., `"Si-O"` |
| `A` | Float | Repulsive pre-exponential (eV) |
| `rho` | Float | Repulsive length scale (A) |
| `C` | Float | Attractive dispersion coefficient (eV A^6) |

### Pedone

V(r) = D0 [1 - exp(-alpha (r - r0))]^2 - D0 + C0/r^12

The Pedone potential (Pedone et al. 2006) combines a Morse-type interaction with a short-range r^-12 repulsion. Commonly used for oxide glasses.

```toml
[[potential.pedone]]
pair = "Ca-O"
D0 = 0.030211
alpha = 2.241334
r0 = 2.923245
C0 = 5.0

[[potential.pedone]]
pair = "Si-O"
D0 = 0.340554
alpha = 2.006700
r0 = 2.100000
C0 = 1.0

[[potential.pedone]]
pair = "O-O"
D0 = 0.042395
alpha = 1.379316
r0 = 3.618701
C0 = 22.0
```

| Field | Type | Description |
|-------|------|-------------|
| `pair` | String | Species pair |
| `D0` | Float | Well depth (eV) |
| `alpha` | Float | Morse width parameter (1/A) |
| `r0` | Float | Equilibrium distance (A) |
| `C0` | Float | r^-12 repulsion coefficient (eV A^12) |

### Coulomb DSF

Damped Shifted Force electrostatics, matching LAMMPS `pair_style coul/dsf`. Automatically generates pair potentials for all species combinations with nonzero charge products.

```toml
[potential.coulomb]
alpha = 0.25                          # DSF damping parameter (1/A)
charges = { Ca = 1.2, Si = 2.4, O = -1.2 }
```

| Field | Type | Description |
|-------|------|-------------|
| `alpha` | Float | DSF damping parameter (1/A), same as LAMMPS `coul/dsf` first argument |
| `charges` | Map | Per-species charges in electron units |

The DSF method provides a smooth, pairwise-additive approximation to Coulomb interactions that goes exactly to zero at the cutoff with zero force, matching the LAMMPS implementation. Each pair's Coulomb table is added to any existing short-range potential for that pair.

### Tabulated

Two-column file (r in A, V in eV), same format as experimental data files. The potential is interpolated onto a fine uniform grid and shifted so V(cutoff) = 0.

```toml
[[potential.tabulated]]
pair = "Ca-O"
file = "CaO_potential.dat"
```

| Field | Type | Description |
|-------|------|-------------|
| `pair` | String | Species pair |
| `file` | String | Path to two-column potential file (relative to config) |

Tabulated potentials **replace** any analytical form defined for the same pair (Buckingham, Pedone, Coulomb). Use this when you have effective potentials from LAMMPS `pair_write` or other sources.

## Combination rules

For each pair, potentials are processed in order:

1. Buckingham (if defined)
2. Pedone (if defined, added to Buckingham)
3. Coulomb DSF (added to whatever exists)
4. Tabulated (**replaces** everything above)

Within a single pair, analytical forms are summed. This means Pedone + Coulomb gives the total interatomic potential -- exactly matching a LAMMPS `hybrid/overlay pedone coul/dsf` setup.

## Complete example: CaSiO3 with Pedone + Coulomb

```toml
[potential]
weight = 0.001
cutoff = 15.0

[[potential.pedone]]
pair = "Ca-O"
D0 = 0.030211
alpha = 2.241334
r0 = 2.923245
C0 = 5.0

[[potential.pedone]]
pair = "Si-O"
D0 = 0.340554
alpha = 2.006700
r0 = 2.100000
C0 = 1.0

[[potential.pedone]]
pair = "O-O"
D0 = 0.042395
alpha = 1.379316
r0 = 3.618701
C0 = 22.0

[potential.coulomb]
alpha = 0.25
charges = { Ca = 1.2, Si = 2.4, O = -1.2 }
```

This matches the LAMMPS setup:
```
pair_style hybrid/overlay pedone 15.0 coul/dsf 0.25 15.0
```
