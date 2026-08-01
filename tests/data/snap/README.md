# SNAP test data

`linear_two_element.snapcoeff` and `.snapparam` are small synthetic,
FitSNAP-compatible files authored for rsmith. They are test data, not a
scientifically meaningful potential.

The example cells exercise a single-element diamond-Si structure and a
two-element zincblende-InP structure. Later numerical reference fixtures are
generated with the corresponding example potentials distributed by LAMMPS:

- `Si_Zuo_JPCA2020.snapcoeff` / `.snapparam`
- `InP_JCPA2020.snapcoeff` / `.snapparam`

Those third-party potential files are not copied into rsmith. Tests that use
them discover a local LAMMPS potential directory or skip with an explicit
message.
