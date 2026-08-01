# SNAP test data

`linear_two_element.snapcoeff` and `.snapparam` are small synthetic,
FitSNAP-compatible files authored for rsmith. They are test data, not a
scientifically meaningful potential. `linear_two_element_si.reference.json`
records LAMMPS descriptors and energies for an equilibrium and a displaced
diamond-Si configuration using this model. The compact JSON fixture is what
rsmith's tests consume.

The fixture is reproducible with a LAMMPS build containing ML-SNAP:

```console
cd tests/data/snap
lmp_serial -log none -in linear_two_element_si.lammps.in
```

This writes `linear_two_element_si.dump`, whose first frame is the equilibrium
cell and whose second frame has atom 1 displaced by the vector documented in
the input. The checked JSON keeps one equilibrium environment and several
distinct displaced environments, plus both total energies. The 2x2x2 cell is
larger than twice the model cutoff so all references use the minimum-image
convention that the native evaluator will enforce.

The example XYZ cells exercise a single-element diamond-Si structure and a
two-element zincblende-InP structure. Additional numerical validation can use
the corresponding example potentials distributed by LAMMPS:

- `Si_Zuo_JPCA2020.snapcoeff` / `.snapparam`
- `InP_JCPA2020.snapcoeff` / `.snapparam`

Those third-party potential files are not copied into rsmith. Tests that use
them discover a local LAMMPS potential directory or skip with an explicit
message.
