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

## Result

The EPSR-native hidden-coordinate replay passed at `6.48e-9` neutron and
`1.94e-8` X-ray RMS. At the nominally matched 30,000-move endpoint, weight-30
rsmith PACE had lower held-out partial-RDF error than stock EPSR26 in all five
paired seeds. The ensemble medians are:

| Method | Common i(Q) RMS | Hidden partial-RDF RMS | Safety passes | Local wall time |
|---|---:|---:|---:|---:|
| rsmith PACE, weight 30 | 0.06758 | 0.43715 | 0/5 | 95.19 s |
| stock EPSR26 | 0.14807 | 0.94698 | 4/5 | 2.53 s |
| rsmith pure RMC | 0.05619 | 0.41053 | 0/5 | 1.79 s |

PACE therefore recovers 59.5% of the starting hidden-RDF gap, compared with
12.4% for EPSR, and its median error is 53.8% lower. This is a strong recovery
advantage, but the preregistered overall-superiority claim does **not** pass:
every PACE endpoint contains at least one distance below the frozen floors,
principally O-O minima of 1.871--1.966 A. EPSR fails one of five endpoints by a
small Si-Si violation. Pure RMC is still more pathological.

The ACE energy decreases by 0.194--0.198 eV/atom in all five PACE runs, so a
favorable global energy does not guarantee safe local minima. Stock EPSR is
also about 37.6 times faster in these local endpoint timings. The timing is a
diagnostic, not a publication speed ratio, because the implementations and
move semantics differ. The result supports superior hidden-structure recovery
for rsmith on this favorable same-ACE synthetic task, while leaving local
safety unresolved and forbidding a general claim that rsmith is superior to
EPSR.

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
python3 scripts/verify_ten_gpa_comparison.py
```

## Incremental pressure ladder

The next frozen pilot separates the available loading path into 0 -> 5 and
5 -> 10 GPa recovery steps. It uses one seed and one 30,000-move endpoint for
pure RMC, PACE weight 30, and stock EPSR26. This is pressure mapping, not a new
ensemble campaign. The protocol is in `expected/pressure-steps.toml`.
The generated inputs and pre-refinement gap audit are pinned in
`expected/pressure-steps-inputs.toml`.

Following the 10 GPa diagnostic review, no ambient absolute minimum-distance
gate is used. Short-range behavior is compared with the hidden target at the
same pressure through the 0.1% and 1% lower quantiles of the Si-Si, Si-O, and
O-O neighbor-distance distributions, together with counts below those target
quantiles. Absolute minima remain descriptive outlier diagnostics only.

```bash
python3 scripts/import_pressure_steps.py
python3 scripts/prepare_pressure_steps.py
python3 scripts/run_pressure_steps_rsmith.py
python3 scripts/run_pressure_steps_epsr.py
python3 scripts/score_pressure_steps.py
```

The upstream repository now contains the sequential loading endpoints through
200 GPa. The next frozen scan uses every 10 GPa increment through 70 GPa and
adds native RMCProfile to the three existing methods. The four arms are pure
rsmith RMC, rsmith PACE weight 30, the stock EPSR26 executable, and native
RMCProfile 6.7.9.5. Thus `native-epsr26` means the original EPSR26 program, not
rsmith's EPSR-inspired objective.

The 0 -> 10 step is rerun rather than copied from the earlier benchmark so all
seven steps share the same 13.5 A real-space scoring cutoff. RMCProfile uses
zero hard minimum distances: ambient-pressure cutoffs are not defensible at
70 GPa, and pressure-relative lower-tail distributions are instead withheld
for scoring. Native EPSR26 and RMCProfile both generate and replay their own
hidden-coordinate targets before refinement. The protocol, including source
commit and executable hashes, is frozen in
`expected/pressure-series-70.toml` before any outcomes are generated.

```bash
python3 scripts/import_pressure_series_70.py
python3 scripts/prepare_pressure_series_70.py
python3 scripts/run_pressure_series_70_rmcprofile.py --forward-only
python3 scripts/run_pressure_series_70_rsmith.py --only-missing
python3 scripts/run_pressure_series_70_epsr.py --only-missing
python3 scripts/run_pressure_series_70_rmcprofile.py --only-missing
python3 scripts/score_pressure_series_70.py
```

This is deliberately a single-seed pressure-trend scan. It can compare how the
four methods respond as pressure rises, but cannot establish universal
superiority, production speed, or stochastic uncertainty.

### Incremental result

Both pressure increments are informative and pass their EPSR-native replay
gates at approximately `1e-8` RMS. At 30,000 nominal moves:

| Step | Method | Common i(Q) RMS | Hidden RDF RMS | Mean lower-tail error |
|---|---|---:|---:|---:|
| 0 -> 5 | PACE weight 30 | 0.03200 | 0.26878 | 0.0712 A |
| 0 -> 5 | stock EPSR26 | 0.08156 | 0.74536 | **0.0451 A** |
| 5 -> 10 | PACE weight 30 | 0.02636 | 0.21446 | 0.0584 A |
| 5 -> 10 | stock EPSR26 | 0.07581 | 0.65301 | **0.0400 A** |

PACE recovers 68.6% and 71.3% of the starting hidden-RDF gaps, whereas EPSR
recovers 12.8% and 12.6%. PACE's hidden-RDF errors are 63.9% and 67.2% lower
than EPSR, so the recovery-advantage gate passes in both steps.

The pressure-relative lower-tail guard does not pass. PACE is closer to the
target for the Si-O lower tail, while EPSR is closer for O-O; the frozen
equal-pair mean therefore favors EPSR. Importantly, incremental PACE produces
minima of at least 2.109 A (Si-Si), 1.461 A (Si-O), and 2.117 A (O-O). The
severe short-contact outliers seen in the direct 0 -> 10 jump disappear. This
supports sequential pressure refinement as the appropriate workflow, while
retaining the lower-tail mismatch as a real structural diagnostic rather than
calling it an absolute safety failure.

Pure RMC obtains slightly lower scattering and hidden-RDF errors than PACE but
again develops very short local contacts and much larger lower-tail errors.
PACE energy decreases by 0.052 and 0.049 eV/atom. The result remains a
single-seed, same-ACE inverse-oracle pressure map; the frozen 0--70 GPa scan
above is the next computation.
