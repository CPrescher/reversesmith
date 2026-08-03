# Official EPSR amorphous-silica benchmark

Amorphous silica is one of the examples selected by the EPSR authors in the
EPSR manual. The comparison must use the official EPSR26 implementation and
its reference-potential plus empirical-potential workflow. rsmith's internal
experimental EPSR mode is not an acceptable stand-in for the comparator.

The official download page states that bundled raw data and demonstration files
may not be used in publicity or publication without express permission. They
may be used locally for execution testing, but they are not vendored here and
must not enter a paper archive until STFC/ISIS grants permission. The official
page, rather than this repository, is authoritative for the current terms.

On 2026-08-02 the user accepted those conditions for local execution testing.
The historical official ZIP targets returned HTTP 404 from the migrated live
site. On 2026-08-03, however, the complete `DTBsilicaNX` worked example was
identified inside the user's existing local EPSR26 installation. The local-only
importer copies it below an ignored directory and records file hashes; no
upstream data are committed or redistributed. See `ACQUISITION_LOG.md`.

The EPSR26 manual's silica exercise fixes several important starting points:

- number density from the worked `.ato`: 0.0663999705 atoms/A^3;
- equilibrate the initially random network at 10,000 K before cooling to 300 K;
- neutron normalization: the worked input's `nrtype = 5` convention;
- X-ray normalization: single-atom scattering;
- use both neutron and X-ray total-scattering data.

Required sequence:

1. obtain written publication/reuse permission before putting upstream data in
   a paper archive or supplement;
2. import an authorized local EPSR26 installation with
   `import_local_reference.sh --accept-local-testing-terms`;
3. record the file hashes, EPSR26 build, reference potentials, density, atom
   count, equilibration length, refinement length, and analysis commands;
4. reproduce native EPSR with ten independent seeds;
5. reproduce that ensemble with rsmith's EPSR mode before making any
   superiority comparison;
6. repeat with rsmith pure RMC, pair-potential HRMC, and MLIP-HRMC at matched
   fit quality;
7. compare fitted totals, partial RDFs, coordination, O-Si-O and Si-O-Si
   angles, ring topology, wall time, and seed-to-seed variability.

The comparison design is preregistered in `protocol.toml`. Numeric stochastic
equivalence limits will be frozen only after the native EPSR replicate spread
is known and before inspecting any rsmith result. This avoids defining
"equivalent" from the answer we want.

## Deterministic two-contrast forward gate

The zero-move gate uses the exact 6,000-atom `DTBsilica.ato`, recalculates one
native EPSR configuration with potential feedback disabled, and evaluates the
same coordinates with rsmith. EPSR's neutron `nrtype=5` and X-ray
single-atom-scattering normalizations are reconstructed from the upstream
weight files before comparing totals; rsmith's default output remains the
Faber-Ziman normalization.

| Curve | RMS / native dynamic range |
|---|---:|
| Si-Si partial i(Q) | 0.0558% |
| Si-O partial i(Q) | 0.0397% |
| O-O partial i(Q) | 0.0441% |
| EPSR-normalized neutron total i(Q) | 0.0639% |
| EPSR-normalized X-ray total i(Q) | 0.0661% |
| Si-Si g(r), rebinned to 0.12 A | 0.858% |
| Si-O g(r), rebinned to 0.12 A | 0.626% |
| O-O g(r), rebinned to 0.12 A | 0.934% |

An independent NumPy minimum-image histogram agrees with each rsmith partial
RDF to better than `6e-8` RMS/range. The committed limits in
`expected/native-forward.toml` are software-regression guards selected after
this diagnostic result; they are not preregistered stochastic equivalence
margins.

## Reference-potential gate and first joint epoch

The EPSR reference potential is defined jointly by the species records in
`si.ato` and `o.ato` and the control values in `DTBsilica.pcof`; the latter is
not, by itself, a complete pair-potential file. Source inspection and an
independent reconstruction establish the tutorial's exact conventions:

- geometric mixing of epsilon and arithmetic mixing of sigma;
- EPSR's modified 12--6 expression (algebraically the usual 12--6 form for
  this case);
- charge-product electrostatics with the EPSR constant of 1390 kJ A/mol and
  zero charge-cloud radii;
- a cosine short-range taper from 9 to 12 A and the separate Hummer Coulomb
  taper at 12 A.

Over 0.5--12 A, the reconstructed Si-Si, Si-O, and O-O curves agree with a
fresh native EPSR `.o01` at `3.25e-8`, `6.97e-8`, and `7.00e-8` RMS of their
dynamic ranges. On the exact 6,000-atom structure, native EPSR reports
-3539 kJ/mol per EPSR molecule (an atom in this input), while rsmith's 0.001 A
tables give -3535.287 kJ/mol per atom: a 0.105% difference. This passes the
frozen 1% energy guard and is consistent with the programs' different
potential-table grids.

The same script converts the measured neutron `nrtype=5` and X-ray `<f^2>`
totals to rsmith's Faber-Ziman normalization and runs one energy-only EPSR
epoch using both contrasts. The 6,000-move, one-thread smoke run completed with
7.0% acceptance; all three empirical pair potentials changed, with a maximum
absolute update of 0.001846 eV. The post-MC X-ray and neutron RMS residuals
were 0.001386 and 0.002137. An output audit corrected an earlier smoke result
that had mistakenly treated EPSR's `.v01` data-minus-model residual as its
calculated total; the calculated total is stored in `.u01`. These values prove
that the full
two-contrast reference-plus-empirical path executes; a single epoch is not a
converged EPSR reproduction or a speed comparison.

The deterministic forward and reference-potential blockers are therefore
closed. The remaining reproduction gate is a matched, multi-epoch native and
rsmith ensemble with frozen seed-to-seed equivalence margins.

## Ambient GAP and Pedone model references

Two 3,000-atom 300 K endpoints from the private `CPrescher/SiO2_glass`
repository are now pinned as secondary structural references at commit
`b05590846c87eb58cf1ed2e09a6787c1d67f9e53`:

- a Pedone-pre-equilibrated structure quenched and held with the Erhard et al.
  silica GAP;
- a Pedone-only structure produced with the same fast schedule.

Both use a single deterministic seed and a `10^14 K/s` quench, so they are
model endpoints rather than equilibrium ensembles. They also adopt different
0-bar model densities. The GAP endpoint is 2.282 g/cm3 (3.33% above the EPSR
worked number density), whereas Pedone is 2.560 g/cm3 (15.94% above it).

An independent minimum-image analysis gives:

| Diagnostic | GAP | Pedone |
|---|---:|---:|
| Si coordination defects, 2.2 A cutoff | 0.40% | 0.40% |
| O coordination defects, 2.2 A cutoff | 0.60% | 0.85% |
| mean O-Si-O angle | 109.40 deg | 109.31 deg |
| mean Si-O-Si angle | 139.21 deg | 148.06 deg |
| Si-Si RDF peak | 2.99 A | 3.11 A |
| Si-O RDF peak | 1.61 A | 1.61 A |
| O-O RDF peak | 2.59 A | 2.57 A |

The independently regenerated RDF table agrees with the source repository's
preserved curves to `1.16e-13` maximum absolute difference.

The analyzer also freezes the complete coordination counts, angle summaries,
minimum distances, and an edge-weighted shortest-cycle diagnostic for the Si
network. These references can reveal whether a fitted ensemble has gross
network defects and provide Pedone/GAP recovery tests. Closeness to GAP is not
itself the paper's superiority criterion: held-out scattering and independently
defined chemical diagnostics remain primary, and RDF comparisons must state
the density mismatch.

The scientific claim is not "smaller chi-squared than EPSR". The meaningful
test is whether MLIP-HRMC gives better held-out or chemically diagnostic
structure at the same data agreement.

## Symmetric GAP/Pedone cross-recovery smoke

The first correctness comparison uses plausible cross-starts rather than a
random 3,000-atom network. For the hidden GAP target, the Pedone endpoint is
isotropically expanded by `1.03914913527` into the GAP box. For the hidden
Pedone target, the GAP endpoint is compressed by `0.962325777945`. No random
displacement, relaxation, or target-potential preprocessing is applied. Exact
noise-free neutron and X-ray targets are generated from the hidden coordinates
on the frozen common grid.

All atom-moving programs start from those same density-rescaled coordinates.
Because EPSR26 and RMCProfile use program-specific reciprocal conventions,
each also generates a hidden-coordinate native target and must reproduce it in
an exact zero-move replay. EPSR passes that replay at `2.7e-8` RMS or better and
RMCProfile at exactly zero. Every recovered coordinate set is subsequently
rescored by the same independent minimum-image analyzer. PDFgui/PDFfit2 is
included as a forward-only PDF control: its target/cross-start differences are
5.26--6.74% RMS of the target PDF range, but PDFgui does not perform atom-wise
RMC recovery.

At 6,000 attempted moves (two moves per atom), the common total-scattering RMS
values are:

| Hidden target / method | Neutron start -> smoke | X-ray start -> smoke |
|---|---:|---:|
| GAP / rsmith RMC | 0.10880 -> 0.09122 | 0.12332 -> 0.10339 |
| GAP / rsmith EPSR | 0.10880 -> 0.10286 | 0.12332 -> 0.12060 |
| GAP / native EPSR26 | 0.10880 -> 0.10689 | 0.12332 -> 0.12246 |
| GAP / RMCProfile | 0.10880 -> 0.10837 | 0.12332 -> 0.12256 |
| Pedone / rsmith RMC | 0.11170 -> 0.08541 | 0.12596 -> 0.10174 |
| Pedone / rsmith EPSR | 0.11170 -> 0.10785 | 0.12596 -> 0.12271 |
| Pedone / native EPSR26 | 0.11170 -> 0.11052 | 0.12596 -> 0.12482 |
| Pedone / RMCProfile | 0.11170 -> 0.10939 | 0.12596 -> 0.12339 |

The Pedone- and GAP-regularized rsmith trajectories are nearly identical to
pure RMC in this smoke. The diagnostic calibration explains why: the frozen
energy weight of `0.001` is roughly three to four orders of magnitude below
the suggested per-move balance (`~27--30` for Pedone and `~0.67` for GAP).
These results validate that the adapters and energy paths execute; they do not
show an MLIP benefit. Production must measure a data-fit/energy-weight Pareto
curve and compare methods at matched data fit.

### Bounded HRMC energy-weight pilot

A frozen one-thread pilot then varied only the regularizer weight for 1,000
attempted moves (one third of a move per atom), retaining the same seed, starts,
targets, fit ranges, uncertainties, and minimum-distance constraints. Progress
below is the fraction of the same-budget pure-RMC reduction in RMS; `1.0` means
the regularized run made the same fit progress as pure RMC.

| Regularizer / weight | GAP-target N/X progress | GAP-target dE/atom | Pedone-target N/X progress | Pedone-target dE/atom |
|---|---:|---:|---:|---:|
| Pedone / 3 | 1.030 / 1.001 | +0.01476 eV | 0.977 / 0.969 | +0.03062 eV |
| Pedone / 10 | 1.063 / 0.974 | +0.01125 eV | 0.917 / 0.888 | +0.02281 eV |
| Pedone / 30 | 0.858 / 0.748 | +0.00430 eV | 0.653 / 0.630 | +0.00730 eV |
| Pedone / 100 | 0.547 / 0.370 | -0.00101 eV | 0.254 / 0.280 | -0.00136 eV |
| GAP / 0.1 | 1.017 / 1.012 | +0.57810 eV | 0.926 / 0.929 | +0.48420 eV |
| GAP / 0.3 | 0.937 / 0.913 | +0.55227 eV | 0.908 / 0.889 | +0.45585 eV |
| GAP / 1.0 | 0.000 / 0.000 | ~0 eV | 0.047 / 0.040 | +0.04112 eV |

Pedone weights 3, 10, and 30 therefore form the low/knee/high bracket for a
matched 6,000-move follow-up. Relative to weight `0.001`, weight 30 suppresses
the positive energy drift by 75--79% while retaining at least 63% of pure-RMC
progress; weight 100 fails the frozen 50% progress guard. GAP shows an abrupt
acceptance cliff: weight 0.3 retains 89--94% of fit progress, whereas weight
1.0 accepts only 0 and 10 of 1,000 moves. The interval 0.3--1.0 must therefore
be refined before selecting a high GAP production weight.

The pilot does not yet show better hidden-target structure from either energy
model. At this short budget the independently scored RDF distance generally
tracks the amount of scattering-fit progress. The defensible result is weight
calibration and detection of the GAP acceptance cliff, not HRMC superiority.

The preregistered follow-up resolves that cliff with weights 0.4, 0.5, 0.6,
and 0.8. Weight 0.4 is the largest value that passes the 50% progress guard in
both directions: it retains 79--85% of pure-RMC neutron/X-ray progress and
reduces positive energy drift by 10--12% relative to weight `0.001`. Weight
0.5 retains 79% progress in the Pedone-target direction but only 37--39% in the
GAP-target direction, so it fails the symmetric rule. The selected GAP
low/knee/high production bracket is therefore 0.1, 0.3, and 0.4.

### Erhard-2024 ACE/PACE validation and calibration

The comparison now also uses the public Si-O ACE potential of Erhard et al.
(2024; article DOI `10.1038/s41467-024-45840-9`, dataset DOI
`10.5281/zenodo.10419194`). Its `SBessel` radial basis is evaluated by
rsmith's native PACE backend. On the 3,000-atom GAP-density cross-start,
official `python-ace` gives `-25257.654797728312` eV and native rsmith prints
`-25257.654798` eV. After moving atom 1 by `(0.03, -0.02, 0.01)` A, the two
energies are `-25257.560323690243` and `-25257.560324` eV. This nonzero-move
oracle closes the energy-evaluation gate; it does not validate forces or
training-domain coverage.

A frozen 1,000-move scan found that weights through 3 behave almost like pure
RMC, weight 10 introduces a measurable energy bias while retaining 84--98% of
pure-RMC fit progress, and weight 30 retains 50--75% in both symmetric
directions. Weight 100 fails the X-ray 50% progress guard. The preregistered
PACE low/knee/high production bracket is therefore 3, 10, and 30. Native PACE
took 2.37--3.08 s per 1,000 moves, about 20--24 times faster than the matched
external GAP/QUIP pilot.

### Single-seed joint-acceptance production smoke

The first 6,000-move production attempt used delayed acceptance. GAP weights
0.1 and 0.3 accepted 0 of 6,000 moves: separating the data and energy tests
discarded favorable cancellation between opposing increments. That failed
optimization gate is retained in `hrmc-production-bracket-failure.toml`; the
scientific comparison uses the intended single Metropolis test on
`delta_chi2 + weight*delta_energy`.

The superseding run evaluates all selected Pedone, GAP, and PACE brackets from
both cross-starts. Progress is normalized to the same-seed 6,000-move pure-RMC
control; values above one are possible from stochastic trajectory differences.

| Regularizer / weight | GAP-target N/X progress | GAP-target RDF RMS | Pedone-target N/X progress | Pedone-target RDF RMS |
|---|---:|---:|---:|---:|
| GAP / 0.1 | 1.004 / 0.999 | 0.61159 | 0.970 / 0.973 | 0.45951 |
| GAP / 0.3 | 0.961 / 0.971 | 0.60968 | 0.969 / 0.966 | 0.46444 |
| GAP / 0.4 | 0.919 / 0.920 | 0.61066 | 0.943 / 0.956 | 0.46448 |
| PACE / 3 | 1.030 / 1.002 | 0.59508 | 1.006 / 1.012 | 0.44652 |
| PACE / 10 | 0.967 / 0.903 | 0.58861 | 0.983 / 0.932 | 0.46371 |
| PACE / 30 | 0.765 / 0.648 | 0.63613 | 0.737 / 0.619 | 0.54503 |
| Pedone / 3 | 1.042 / 0.988 | 0.59909 | 0.996 / 1.009 | 0.44614 |
| Pedone / 10 | 0.992 / 0.932 | 0.58353 | 0.986 / 0.984 | 0.45205 |
| Pedone / 30 | 0.830 / 0.692 | 0.60258 | 0.744 / 0.758 | 0.51642 |
| Pure-RMC control | 1.000 / 1.000 | 0.60126 | 1.000 / 1.000 | 0.45858 |

PACE weight 3 is the balanced PACE setting: it preserves essentially all fit
progress and improves hidden-target mean partial-RDF RMS by 1.0% and 2.6% in
the two directions. Pedone weight 3 gives a comparably small improvement, so
this single-seed synthetic smoke does **not** establish a PACE-specific
structural advantage. It establishes correctness, a usable PACE weight, and
practical native throughput. The matched 6,000-move PACE arms took 12.3--16.6
s versus 307.7--354.8 s for GAP/QUIP, a 21--25 times speedup. Repeated timing
and multi-seed structural ensembles remain publication gates.

### Local timing diagnostic

End-to-end one-thread wall times on an Apple M4 Pro were:

| Method, joint fit | GAP target | Pedone target |
|---|---:|---:|
| rsmith RMC | 0.529 s | 0.522 s |
| rsmith EPSR mode | 0.385 s | 0.417 s |
| native EPSR26 | 1.085 s | 1.116 s |
| RMCProfile | 6.920 s | 6.709 s |
| rsmith GAP/QUIP | 321.9 s | 343.3 s |

rsmith pure RMC is about 13 times faster than this installed RMCProfile path
and about twice as fast as native EPSR in the smoke. Those are not yet paper
speed claims. The official RMCProfile macOS wrapper fell back to its serial
binary because its parallel executable references an unavailable x86-64
`libgomp`; it also checkpoints coordinates by wall-clock time, so the
independently scored snapshots represent 5,315 and 5,490 of the 6,000 completed
moves. Conversely, GAP/QUIP model loading and local SOAP evaluation dominate
and make the present GAP-HRMC path roughly 600 times slower than pure RMC.
Production timing needs repeated cold/warm runs, a working parallel
RMCProfile build, exact saved move counts, and calibrated energy weights.

## Files

- `protocol.toml`: frozen methods, stages, contrasts, seeds, and observables;
- `PERMISSION_REQUEST.md`: ready-to-send publication/reuse request;
- `fetch_reference.sh`: explicit, locally gated acquisition with SHA-256 log;
- `import_local_reference.sh`: local-only DTBsilicaNX importer and hasher;
- `import_ambient_models.sh`: local-only importer for the pinned private
  `SiO2_glass` GAP and Pedone endpoints;
- `scripts/run_native_zero_move.py`: clean one-configuration EPSR26 runner;
- `scripts/compare_forward.py`: exact-structure rsmith comparison and
  independent RDF oracle;
- `scripts/verify_forward.py`: committed forward-regression verifier;
- `scripts/run_reference_potential_smoke.py`: independent EPSR reference-
  potential reconstruction, native curve/energy gate, normalization
  conversion, and one-epoch joint neutron/X-ray rsmith smoke run;
- `scripts/fetch_public_gap.py`: checksum-gated importer for the public Erhard
  silica GAP and its sparse companions;
- `scripts/fetch_public_ace2024.py` and `verify_pace2024_oracle.py`:
  checksum-gated Erhard-2024 ACE acquisition and independent energy-oracle
  verification;
- `scripts/prepare_cross_recovery.py`: symmetric hidden-target and cross-start
  fixture generator for the common neutron/X-ray protocol;
- `scripts/run_rsmith_cross.py`, `run_native_epsr_cross.py`, and
  `run_rmcprofile_cross.py`: matched adapter smoke runners;
- `scripts/run_pdfgui_cross_forward.py`: PDFgui/PDFfit2 forward-only control;
- `scripts/score_cross_recovery.py` and `verify_cross_recovery.py`: common
  coordinate scoring and committed smoke guards;
- `scripts/prepare_hrmc_weight_sweep.py`, `run_hrmc_weight_sweep.py`, and
  `verify_hrmc_weight_sweep.py`: frozen Pedone/GAP weight-pilot preparation,
  execution, scoring, and regression guards;
- `scripts/prepare_hrmc_gap_cliff_followup.py` and
  `verify_hrmc_gap_cliff_followup.py`: additive GAP cliff-grid preparation and
  its symmetric selection guard;
- `scripts/prepare_hrmc_pace2024_weight_pilot.py` and its runner/verifier:
  frozen native-PACE weight calibration;
- `scripts/prepare_hrmc_production_unscreened.py` and its runner/verifier:
  matched joint-acceptance Pedone/GAP/PACE production smoke;
- `scripts/analyze_ambient_models.py`: independent RDF, coordination, angle,
  minimum-distance, and shortest-ring analysis of both model endpoints;
- `scripts/verify_ambient_models.py`: deterministic verifier for the pinned
  structure hashes and derived observations;
- `expected/native-forward.toml`: provenance hashes, observations, and
  regression guards;
- `expected/reference-potential-smoke.toml`: frozen potential and smoke-run
  observations plus regression guards;
- `expected/ambient-model-endpoints.toml`: pinned source hashes, model
  provenance, structural observations, and claim boundary;
- `expected/cross-recovery.toml`: protocol frozen before cross-refinement
  outcomes were inspected;
- `expected/cross-recovery-smoke.toml`: post-diagnostic observations and
  adapter regression guards, explicitly not publication equivalence limits;
- `expected/hrmc-weight-sweep.toml` and `hrmc-weight-sweep-smoke.toml`: frozen
  pilot design, observed Pareto brackets, timings, and claim boundary;
- `expected/hrmc-gap-cliff-followup.toml` and its `-smoke.toml` record: frozen
  fine-grid design and the selected GAP 0.1/0.3/0.4 production bracket;
- `expected/pace2024-oracle.toml`: pinned model and official-python-ace energy
  oracle;
- `expected/hrmc-pace2024-weight-*.toml`: frozen PACE calibration protocols and
  observed bracket;
- `expected/hrmc-production-*.toml`: the retained delayed-acceptance failure,
  superseding joint-acceptance protocol, and observed three-model comparison;
- `reference/README.md`: provenance and redistribution rules for upstream data.

## Local reproduction

```bash
cd benchmarks/epsr-silica
EPSR26_ROOT=/path/to/EPSR26/EPSR ./import_local_reference.sh --accept-local-testing-terms
python3 scripts/run_native_zero_move.py
cargo build --release --manifest-path ../../Cargo.toml
python3 scripts/compare_forward.py
python3 scripts/verify_forward.py
python3 scripts/run_reference_potential_smoke.py

SIO2_GLASS_ROOT=/path/to/SiO2_glass ./import_ambient_models.sh
python3 scripts/analyze_ambient_models.py
python3 scripts/verify_ambient_models.py

python3 scripts/fetch_public_gap.py
python3 scripts/prepare_cross_recovery.py --force \
  --gap-model reference/local/public-gap/silica_gap.xml
python3 scripts/run_rsmith_cross.py --no-gap
# Build rsmith with --features gap-quip, then pass that binary:
python3 scripts/run_rsmith_cross.py --gap-only --binary /path/to/gap-enabled/rsmith
python3 scripts/run_native_epsr_cross.py --force
python3 scripts/run_rmcprofile_cross.py --force
/path/to/pdfgui/python scripts/run_pdfgui_cross_forward.py
python3 scripts/score_cross_recovery.py
python3 scripts/verify_cross_recovery.py

python3 scripts/prepare_hrmc_weight_sweep.py --force
python3 scripts/run_hrmc_weight_sweep.py --model pure
python3 scripts/run_hrmc_weight_sweep.py --model pedone
python3 scripts/run_hrmc_weight_sweep.py --model gap \
  --binary /path/to/gap-enabled/rsmith
python3 scripts/score_cross_recovery.py
python3 scripts/verify_hrmc_weight_sweep.py

python3 scripts/prepare_hrmc_gap_cliff_followup.py
python3 scripts/run_hrmc_weight_sweep.py --model gap --only-missing \
  --binary /path/to/gap-enabled/rsmith
python3 scripts/score_cross_recovery.py
python3 scripts/verify_hrmc_gap_cliff_followup.py

python3 scripts/fetch_public_ace2024.py
python3 scripts/verify_pace2024_oracle.py
python3 scripts/prepare_hrmc_pace2024_weight_pilot.py --force
python3 scripts/run_hrmc_pace2024_weight_pilot.py
python3 scripts/prepare_hrmc_pace2024_weight_followup.py
python3 scripts/run_hrmc_pace2024_weight_pilot.py --only-missing
python3 scripts/score_cross_recovery.py
python3 scripts/verify_hrmc_pace2024_weight_pilot.py

python3 scripts/prepare_hrmc_production_unscreened.py --force
python3 scripts/run_hrmc_production_unscreened.py --model pedone
python3 scripts/run_hrmc_production_unscreened.py --model pace
python3 scripts/run_hrmc_production_unscreened.py --model gap \
  --binary /path/to/gap-enabled/rsmith
python3 scripts/score_cross_recovery.py
python3 scripts/verify_hrmc_production_unscreened.py
```
