# rsmith

A Reverse Monte Carlo (RMC) structure refinement tool written in Rust. Refines atomic structures against experimental X-ray and neutron scattering data while enforcing physical constraints and pair potentials.

## Features

- **S(Q) and g(r) fitting** -- simultaneous refinement against structure factor and pair distribution function
- **Hybrid RMC** -- pair potentials (Buckingham, Pedone, Coulomb DSF, tabulated) bias refinement toward energetically favorable configurations
- **Hard constraints** -- minimum interatomic distances and coordination number bounds
- **Simulated annealing** -- exponential cooling schedule with adaptive step size
- **Structural analysis** -- coordination numbers and bond angle distributions for validation
- **Incremental updates** -- O(N_neighbors) per move via cell lists and precomputed lookup tables

## Quick Start

```bash
cargo build --release
./target/release/rsmith config.toml
```

Minimal `config.toml`:

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

[sq]
rdf_cutoff = 11.0

[constraints.min_distance]
"Si-O" = 1.2
"O-O" = 2.0
"Ca-O" = 1.8
```

## Modes

```bash
rsmith config.toml                      # RMC refinement
rsmith config.toml --compute-sq-only    # Compute S(Q) only
rsmith config.toml --analyze            # Structural analysis
```

## Documentation

Full documentation is available at [rsmith book](docs/src/SUMMARY.md) or build locally:

```bash
cd docs && mdbook serve
```

## References

- McGreevy, R.L. & Pusztai, L. (1988). Reverse Monte Carlo Simulation. *Mol. Simul.*, 1, 359-367.
- Pedone, A. et al. (2006). A New Self-Consistent Empirical Interatomic Potential Model for Oxides. *J. Phys. Chem. B*, 110, 11780-11795.
