# Reversesmith

Reversesmith is a Reverse Monte Carlo (RMC) structure refinement tool written in Rust. It refines atomic structures against experimental X-ray and neutron scattering data while enforcing physical constraints.

## What it does

Reversesmith iteratively displaces atoms in a model structure to minimize the difference between computed and experimental scattering functions. It supports:

- **S(Q)** fitting -- total X-ray or neutron structure factor in reciprocal space
- **g(r)** fitting -- total X-ray pair distribution function in real space
- **Pair potential constraints** -- Buckingham, Pedone, Coulomb DSF, or tabulated potentials to bias refinement toward energetically favorable configurations (hybrid RMC / EPSR-like)
- **Structural constraints** -- minimum interatomic distances and coordination number bounds
- **Simulated annealing** -- escape local minima and optimize systematically
- **Structural analysis** -- coordination numbers and bond angle distributions for validation

Both S(Q) and g(r) targets can be fitted simultaneously with independent weights and fit ranges.

## When to use it

Reversesmith is designed for refining atomic models of disordered materials (glasses, liquids, amorphous solids) against total scattering data. A typical workflow:

1. Generate an initial structure via molecular dynamics (e.g., melt-quench in LAMMPS)
2. Obtain experimental S(Q) and/or g(r) from X-ray or neutron diffraction
3. Refine the MD structure against the experimental data
4. Analyze the refined structure for local order (coordination, bond angles)

The pair potential feature (hybrid RMC) is particularly useful when pure RMC produces physically unreasonable structures -- it keeps the refinement close to the potential energy surface of the force field used to generate the initial structure.
