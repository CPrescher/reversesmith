# Getting Started

## Installation

```bash
# From source
git clone https://github.com/CPrescher/rsmith.git
cd rsmith
cargo build --release
```

The binary is at `target/release/rsmith`.

## Usage

```bash
# Run RMC refinement
rsmith config.toml

# Analyze coordination and bond angles
rsmith config.toml --analyze

# Analyze a specific structure
rsmith config.toml --analyze structure.xyz

# Write output to a specific directory
rsmith config.toml --output-dir run01

# Suppress terminal output (log file only)
rsmith config.toml --quiet

# Override the RNG seed
rsmith config.toml --seed 123

# Resume an interrupted run from checkpoint
rsmith config.toml --resume

# Resume with a different output directory
rsmith config.toml --resume --output-dir run01
```

### Seed handling

The RNG seed controls the random atom displacements. Priority: `--seed` flag > `seed` in `[rmc]` > random.

For production runs, **omit `seed` from `[rmc]`** so each run gets a unique random seed. The seed is always logged in `rsmith.log`, so any run can be reproduced later by passing `--seed <value>`. The `--seed` flag and the config `seed` field are mainly useful for debugging and testing.

### Parallel ensemble runs

RMC refinement can converge to different local minima depending on the random trajectory. Running multiple independent fits and comparing results gives confidence in the solution. Omit `seed` from `[rmc]` to ensure each run explores a different path:

```bash
for i in $(seq 1 8); do
    rsmith config.toml --output-dir run$(printf '%02d' $i) --quiet &
done
wait
```

Each run gets a random seed (logged in its `rsmith.log`) and writes all output to its own directory. Analyze the results afterwards:

```bash
for i in $(seq 1 8); do
    rsmith config.toml --analyze run$(printf '%02d' $i)/refined.xyz --output-dir run$(printf '%02d' $i)
done
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
# sigma = 0.02              # Omit to auto-estimate from data noise
# sigma_alpha = 0.05     # Relax fit at high Q (optional)

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
rsmith config.toml
```

This refines `glass.data` against `experimental.sq` and writes:
- `start_xray_sq.dat` / `start_gr.dat` -- total X-ray S(Q) and partial g(r) of the starting structure
- `refined.xyz` -- refined structure
- `refined_xray_sq.dat` / `refined_gr.dat` -- total X-ray S(Q) and partial g(r) of the refined structure

## Input File Formats

**LAMMPS data file** (`format = "lammps"`): Reads the `Atoms # charge` section. Expects columns: `atom_id type charge x y z [ix iy iz]`. Box bounds from `xlo xhi` / `ylo yhi` / `zlo zhi` lines. Requires `[system.types]` mapping type IDs to element names.

**XYZ file** (`format = "xyz"`): Standard XYZ with box dimensions in the comment line, either as `Lattice="Lx 0 0 0 Ly 0 0 0 Lz"` (extended XYZ) or `Lx Ly Lz`.

**Experimental data**: Two-column whitespace-separated files. Lines starting with `#` are skipped. First column is Q (1/A) or r (A), second column is S(Q) or g(r).
