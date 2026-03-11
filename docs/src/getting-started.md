# Getting Started

## Installation

```bash
# From source
git clone https://github.com/CPrescher/reversesmith.git
cd reversesmith
cargo build --release
```

The binary is at `target/release/reversesmith`.

## Usage

```bash
# Run RMC refinement
reversesmith config.toml

# Compute S(Q) from a structure (no refinement)
reversesmith config.toml --compute-sq-only

# Analyze coordination and bond angles
reversesmith config.toml --analyze

# Analyze a specific structure
reversesmith config.toml --analyze structure.xyz
```

## Minimal Example

Create a `config.toml`:

```toml
[system]
structure = "glass.data"
format = "lammps"

[system.types]
1 = "Ca"
2 = "Si"
3 = "O"

[data.xray_sq]
file = "experimental.sq"
sigma = 0.02

[rmc]
max_moves = 500_000
max_step = 0.1

[sq]
qmax = 20.0
rdf_cutoff = 20.0     # Use 20+ A for accurate first peak in S(Q)

[constraints.min_distance]
"Si-O" = 1.2
"O-O" = 2.0
"Ca-O" = 1.8
"Ca-Si" = 2.2
"Ca-Ca" = 2.3
"Si-Si" = 2.3
```

Run it:

```bash
reversesmith config.toml
```

This refines `glass.data` against `experimental.sq` and writes:
- `refined.xyz` -- refined structure
- `refined_sq.dat` -- computed S(Q) of the refined structure

## Input File Formats

**LAMMPS data file** (`format = "lammps"`): Reads the `Atoms # charge` section. Expects columns: `atom_id type charge x y z [ix iy iz]`. Box bounds from `xlo xhi` / `ylo yhi` / `zlo zhi` lines. Requires `[system.types]` mapping type IDs to element names.

**XYZ file** (`format = "xyz"`): Standard XYZ with box dimensions in the comment line, either as `Lattice="Lx 0 0 0 Ly 0 0 0 Lz"` (extended XYZ) or `Lx Ly Lz`.

**Experimental data**: Two-column whitespace-separated files. Lines starting with `#` are skipped. First column is Q (1/A) or r (A), second column is S(Q) or g(r).
