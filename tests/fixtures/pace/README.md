# Native PACE oracle fixtures

The `.yace` files in this directory are compact, repository-authored test
models. They are not production potentials. Their coefficients were chosen to
exercise radial bases, nonlinear embeddings, core energy, multiple elements,
directed cutoffs, and product ranks 1, 2, and 4.

The expected energies in `tests/pace_reference_fixture.rs` were generated with
LAMMPS 22 Jul 2025 Update 4, `pair_style pace product`, and ML-PACE
2023.11.25. The accompanying LAMMPS inputs make each reference reproducible.
The ignored live integration test uses `lammps_cu_cluster.in`.
