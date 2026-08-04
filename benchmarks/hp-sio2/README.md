# High-pressure SiO2 recovery benchmark

This benchmark asks whether atomistic refinement can recover the structural
changes produced by compressing an ambient silica glass to 10 GPa. It begins
as a deliberately small synthetic pilot, not as another broad comparison
matrix.

The hidden target is the 3,000-atom Erhard-2024 ACE structure after sequential
0 -> 5 -> 10 GPa loading at 300 K and a 30 ps hold at 10 GPa. The common start
is the saved 0 GPa structure from the same trajectory, mapped by fractional
coordinates into the exact hidden-target box. Consequently the density is
known and fixed; the refinement must recover non-affine topology rather than
infer pressure or volume.

Synthetic neutron and X-ray totals are calculated from the hidden coordinates.
Hidden partial RDFs and Si/O coordination distributions are withheld from the
refinement and used only for scoring. The first gate compares pure rsmith RMC
with rsmith ACE-HRMC at 6,000, 18,000, and 30,000 attempted moves using one
seed. Native EPSR26 is added only if the pilot demonstrates both a measurable
0-to-10 GPa structural gap and useful recovery. This staged rule avoids the
large ambient-pressure matrix being repeated unnecessarily.

Because the hidden target was generated with the same ACE model available to
ACE-HRMC, this is an intentionally favorable synthetic recovery/oracle test.
It can validate algorithmic use of an informative MLIP prior, but cannot by
itself establish experimental accuracy or universal superiority over EPSR.

The protocol is frozen in `expected/ten-gpa-pilot.toml`, and the deterministically
generated fixture is pinned before refinement in
`expected/ten-gpa-pilot-inputs.toml`. Private source structures remain ignored;
hashes, repository commit, generation input, and the public ACE provenance are
committed.

The pilot passed its extension gate. Its observations are recorded in
`expected/ten-gpa-pilot-observed.toml`; the matched-budget EPSR26 comparison is
separately preregistered in `expected/ten-gpa-comparison.toml`. Only the 30,000
move endpoint is replicated over five seeds. The primary seed alone retains
6,000/18,000/30,000-move checkpoints.

## Reproduction

```bash
cd benchmarks/hp-sio2
python3 scripts/import_private_sources.py
python3 scripts/prepare_ten_gpa_pilot.py
python3 scripts/run_ten_gpa_pilot.py
python3 scripts/score_ten_gpa_pilot.py
python3 scripts/prepare_ten_gpa_comparison.py
python3 scripts/run_ten_gpa_comparison_rsmith.py
python3 scripts/run_ten_gpa_comparison_epsr.py
python3 scripts/score_ten_gpa_pilot.py
```
