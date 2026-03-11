# `[system]` -- Structure Input

```toml
[system]
structure = "structure.data"   # Path to atomic structure file
format = "lammps"              # "lammps" or "xyz"
density = 3.27                 # Optional: target mass density in g/cm^3

[system.types]                 # Required for LAMMPS format: type ID -> element
1 = "Ca"
2 = "Si"
3 = "O"
```

## Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `structure` | String | Yes | Path to the input structure file |
| `format` | String | Yes | File format: `"lammps"` or `"xyz"` |
| `types` | Map | LAMMPS only | Maps LAMMPS type integers to element names |
| `density` | Float | No | Target mass density (g/cm^3). When set, the box and atom positions are rescaled isotropically before any computation begins. |

## Density rescaling

When `density` is set, reversesmith computes the current mass density of the input structure and derives an isotropic scale factor `(current/target)^(1/3)`. Both box dimensions and all atom positions are multiplied by this factor. This is useful for high-pressure simulations where the input structure's box size doesn't match the desired density.

The rescaling happens before S(Q)/g(r) computation and before RMC refinement, so all output reflects the target density.

## LAMMPS format

Reads the `Atoms # charge` section with columns: `atom_id type charge x y z [ix iy iz]`. Box bounds are parsed from `xlo xhi` / `ylo yhi` / `zlo zhi` lines.

The `[system.types]` mapping is required to associate LAMMPS numeric types with element names. The order determines the internal type indexing.

## XYZ format

Standard XYZ with box dimensions in the comment line (line 2):

- Extended XYZ: `Lattice="Lx 0 0 0 Ly 0 0 0 Lz"`
- Simple: `Lx Ly Lz`

Species are determined automatically from the element names in the coordinate lines. The `[system.types]` section is not needed.
