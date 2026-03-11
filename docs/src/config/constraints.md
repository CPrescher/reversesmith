# `[constraints]` -- Physical Constraints

Constraints are hard limits -- moves that violate them are rejected immediately, before computing chi2. This is much cheaper than evaluating the scattering function.

## Minimum distances

```toml
[constraints.min_distance]
"Ca-O" = 1.8
"Si-O" = 1.2
"O-O" = 2.0
"Ca-Si" = 2.2
"Ca-Ca" = 2.3
"Si-Si" = 2.3
```

No atom pair of the given type can be closer than the specified distance (in A). Pair order does not matter (`"Ca-O"` and `"O-Ca"` are equivalent).

These should be set to physically reasonable values -- slightly below the first peak of the corresponding partial g(r). Setting them too tight prevents the RMC from exploring; setting them too loose allows unphysical close contacts.

## Coordination constraints

```toml
[[constraints.coordination]]
pair = "Si-O"
min = 3          # Minimum coordination number
max = 6          # Maximum coordination number
cutoff = 2.2     # Neighbour search radius (A)
```

Multiple `[[constraints.coordination]]` blocks can be defined. After a proposed move, the coordination number of the moved atom **and its affected neighbours** is checked against the bounds.

| Field | Type | Description |
|-------|------|-------------|
| `pair` | String | Species pair, e.g., `"Si-O"` (counts O around Si) |
| `min` | Integer | Minimum allowed coordination number |
| `max` | Integer | Maximum allowed coordination number |
| `cutoff` | Float | Neighbour search radius (A) |

### Example: CaSiO3 glass

```toml
[[constraints.coordination]]
pair = "Si-O"
min = 3
max = 6
cutoff = 2.2
```

This ensures every Si has between 3 and 6 oxygen neighbours within 2.2 A. For CaSiO3 glass, Si is predominantly 4-fold coordinated (tetrahedral), so `[3, 6]` allows some flexibility while preventing extreme distortions.

## Constraints vs. potentials

Constraints are binary (accept/reject) while [pair potentials](potentials.md) are continuous (energy penalty). They serve complementary roles:

- **Constraints** enforce absolute limits: "no Si-O bond shorter than 1.2 A, ever"
- **Potentials** provide a smooth bias: "configurations near the potential energy minimum are preferred"

Using both together is recommended -- constraints prevent catastrophic distortions cheaply, while potentials guide the refinement toward physically reasonable structures.
