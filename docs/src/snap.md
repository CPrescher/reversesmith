# Native SNAP

rsmith can evaluate SNAP potentials fitted by FitSNAP directly in Rust. The
runtime consumes the same `.snapcoeff` and `.snapparam` pair used by LAMMPS,
but it does not launch or link against LAMMPS. Fit the potential with FitSNAP,
copy the two output files beside the rsmith configuration, and select the
`snap_native` backend.

## Quick start

Build rsmith normally; SNAP does not require a Cargo feature or an additional
runtime:

```console
cargo build --release
```

Given this layout:

```text
run/
├── config.toml
├── start.data
├── my_model.snapcoeff
└── my_model.snapparam
```

add the backend to `config.toml`:

```toml
[system]
structure = "start.data"
format = "lammps"

[system.types]
1 = "In"
2 = "P"

# Experimental data, [rmc], [sq], and constraints go here.

[ml_potential]
backend = "snap_native"
coefficient_file = "my_model.snapcoeff"
parameter_file = "my_model.snapparam"
weight = 0.001
```

The file paths are resolved relative to `config.toml`; absolute paths also
work. The element names in `[system.types]` must match names in the coefficient
file exactly. Multiple rsmith atom types may map to the same coefficient-file
element, but every configured type must have a matching element.

Start the run in the usual way:

```console
./target/release/rsmith run/config.toml
```

At startup rsmith prints the two resolved model paths, the energy weight, and
the maximum pair cutoff derived from the files. File, element, coefficient, and
cell compatibility is checked before the RMC loop begins.

## Choosing the weight

The ML potential is an energy regularizer in the RMC acceptance cost:

```text
delta_cost = delta_chi2 + weight * delta_E_SNAP
```

SNAP energies are interpreted in eV, so `weight` controls their influence
relative to the scattering-data residual. The default is `0.001`. A zero
weight is valid and is useful for checking model loading and reported energies;
production values should be chosen for the scale and purpose of the refinement.

`[ml_potential]` cannot currently be combined with `[potential]` or `[epsr]`.
The backend is available for normal RMC runs, not fitting: FitSNAP remains
responsible for creating the potential.

## Model and cell requirements

- Supply a matching FitSNAP/LAMMPS `.snapcoeff` and `.snapparam` pair.
- Model values are interpreted in LAMMPS `metal` units: distances in Angstrom
  and energies in eV. If the coefficient file declares `UNITS:`, it must be
  `metal`; including that metadata is recommended.
- The simulation cell must be orthorhombic and periodic. Every box length must
  be strictly greater than twice the largest active pair cutoff.
- Do not add `cutoff` to `[ml_potential]`. For atom types `a` and `b`, the
  backend derives `R_ab = rcutfac * (R_a + R_b)`,

  where the radii come from the coefficient file.
- The input structure should represent the same chemical domain and units as
  the structures used to fit and validate the potential. rsmith validates file
  shape and numerical invariants, but it cannot determine whether a potential
  is physically transferable to a new structure.

The box-size rule is a deliberate minimum-image restriction. If a model's
maximum cutoff is 4.8 Angstrom, each box length must be greater than 9.6
Angstrom. Create a larger supercell rather than changing the model cutoff.

## Supported SNAP models

The native evaluator supports standard linear and quadratic SNAP energies and
explicit multi-element chemical SNAP. The following model controls are
implemented:

- `rcutfac`, `twojmax`, `rfac0`, and `rmin0`;
- `switchflag` and inner switching through `switchinnerflag`, `sinner`, and
  `dinner`;
- `bzeroflag`, `bnormflag`, and both `wselfallflag` conventions;
- `chemflag` with ordered element-density channels;
- `quadraticflag`, including chemical plus quadratic models.

The LAMMPS execution-tuning keys `chunksize` and `parallelthresh` are accepted
but ignored because they do not change the mathematical model. Unknown model
parameters are rejected instead of being guessed.

Current limitations are intentional:

- rsmith evaluates energies and single-atom energy differences only; it does
  not expose SNAP forces, virials, or stresses;
- rsmith does not fit or refit coefficients;
- triclinic cells are not supported by the native minimum-image evaluator;
- only the SNAP contribution is loaded. Extra LAMMPS pair styles or overlays,
  such as a separate ZBL term, are not part of the `.snapcoeff` file and are
  therefore not evaluated.

## Checking a newly fitted potential

Before a long refinement:

1. Evaluate several representative cells with the pure SNAP contribution in
   LAMMPS and record the total energies.
2. Load the same files, type-to-element mapping, cell, and coordinates in
   rsmith. Confirm the startup-derived maximum cutoff and initial SNAP energy.
3. Include both equilibrium-like structures and displaced atoms. This checks
   energy differences rather than only a constant offset.
4. For a chemical model, include distinct central-element environments. For a
   quadratic model, include configurations that exercise both diagonal and
   cross terms.

LAMMPS is useful for this one-time validation but is not needed for subsequent
rsmith runs. The repository's `tests/data/snap` directory contains synthetic
models, diamond-Si and zincblende-InP cells, reproducible LAMMPS inputs, and
frozen descriptor and energy references. Run the native regression tests with:

```console
cargo test --test snap_model_files --test snap_reference_fixture
```

If LAMMPS example potentials are installed, the optional published-model checks
can be enabled with:

```console
RSMITH_LAMMPS_POTENTIALS=/path/to/lammps/potentials \
  cargo test --test snap_model_files
```

## Common startup errors

| Error | Resolution |
|-------|------------|
| Missing `coefficient_file` or `parameter_file` | Add both paths under `[ml_potential]`. |
| `cutoff must be omitted` | Remove `cutoff`; it is calculated from `rcutfac` and the element radii. |
| Element is absent from the coefficient file | Make `[system.types]` use exactly the coefficient-file element names. |
| Coefficient count does not match model flags | Use the `.snapcoeff` and `.snapparam` produced by the same FitSNAP fit. |
| Box length must exceed twice the maximum cutoff | Replicate the input structure into a larger orthorhombic supercell. |
| Unsupported SNAP parameter | Remove a LAMMPS-only option only after confirming it is execution-only, or open an issue with the model pair. |

For the internal evaluation and caching design, see [Native SNAP
Implementation](./algorithms/snap-native.md).
