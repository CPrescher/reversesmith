# EPSR26 liquid-Ga benchmark

This benchmark starts with the simplest complete EPSR26 workshop case:
monatomic liquid Ga at 323 K. It is a reproduction gate, not yet a claim that
`rsmith` is superior to EPSR.

The first comparison is deliberately deterministic. Native EPSR and `rsmith`
calculate the Ga-Ga RDF and neutron interference function from the identical
5,000-atom tutorial configuration with zero Monte Carlo moves. EPSR's existing
worked output averages 702 configurations and therefore is not used as the
pointwise forward reference. A fresh native run resets the accumulators, uses
all 5,000 atoms as RDF origins, disables resolution broadening, and executes
one output iteration with zero MC cycles.

EPSR output conventions were verified against its bundled Fortran source:

- `.EPSR.f01` is the simulated intermolecular partial, here `S(Q)-1`;
- `.EPSR.u01` is the weighted total-data fit;
- `.EPSR.g01` is the simulated Ga-Ga partial RDF;
- `.EPSR.t01` is the interpolated input data.

For this monatomic neutron case, `.u01 = 0.531149447 * .f01`. The native
partial can therefore be compared directly with rsmith's `convention = "iq"`
output without a scattering-length normalization ambiguity.

## Local reproduction

The raw workshop data are excluded from Git because their disclaimer limits
them to tutorial/example use. After accepting those terms locally:

```bash
EPSR26_ROOT=/path/to/EPSR ./import_local_reference.sh --accept-local-testing-terms
python3 scripts/run_native.py
cargo build --release
python3 scripts/run_comparison.py
```

The frozen pre-result tolerances are in `manifest.toml`. RDF scoring uses
0.12 A blocks so that different native histogram point labels do not dominate
the result. The Q-space comparison is secondary because the two programs use
different real-space quadrature representations.

## Result

The frozen deterministic gate passes. Across 590 Q points from 0.5 to
29.95 inverse angstrom, the rsmith/native-EPSR RMS difference is 0.0381% of
the native `S(Q)-1` dynamic range and the maximum difference is 0.141%. For
the 0.12 A-rebinned RDF, the RMS and maximum differences are 0.752% and 7.05%
of the native range. An independent NumPy minimum-image pair histogram agrees
with rsmith to `9.41e-8` RMS/range. The observed values and the pre-result
limits are recorded in `expected/native-forward.toml`.

## Fixed-budget stochastic comparison

The next gate uses ten seeds and the same supplied 5,000-atom endpoint as the
starting configuration for both methods. Each replicate runs 40 empirical-
potential epochs with 25,000 attempted moves per epoch: 1,000,000 attempted
moves in total. Both executables use one thread. Native EPSR accumulates 40
configurations internally; the 40 rsmith epoch configurations are averaged by
an independent NumPy implementation on the same grids. The latter
post-processing is excluded from executable timing.

This protocol is a matched-start stability/re-equilibration test. It does not
test basin independence or convergence from a random initial model. Its
settings were selected after a pilot and then held fixed for the ten-seed
ensemble, so the results are observations rather than preregistered acceptance
limits.

Run the complete sequential comparison with:

```bash
cargo build --release
python3 scripts/run_stochastic.py --force --jobs 1
```

Existing local runs can be rescored without running either program:

```bash
python3 scripts/run_stochastic.py --summarize-only
```

### Accuracy result

All entries below are median [minimum, maximum] over ten seeds. Residuals are
RMS error divided by the dynamic range of the reference curve.

| Metric | Native EPSR26 | rsmith |
|---|---:|---:|
| Fit to common measured total, Q = 1.3--29.95 inverse angstrom | 2.128% [2.081%, 2.215%] | **1.350% [1.321%, 1.390%]** |
| g(r) versus supplied 702-configuration worked ensemble, r = 1.5--10 A | **2.166% [2.032%, 2.321%]** | 2.478% [2.466%, 2.504%] |
| Final move acceptance | 25.30% [24.12%, 26.54%] | 24.15% [23.60%, 25.90%] |

At the equal move budget, rsmith consistently fits the common measured total
more closely. Native EPSR remains slightly closer to EPSR's supplied long-run
structural ensemble. The median native-versus-rsmith differences are 1.711%
for the interference function and 3.247% for g(r), normalized by the native
dynamic range. The result therefore supports competitive refinement and a
different sampled ensemble; it does not by itself prove that the rsmith
structure is more physically correct. Coordination, angular, topological,
energetic, and held-out-data tests are still required for that claim.

Native EPSR scores its `.u01` total against the interpolated `.t01` data,
whereas rsmith is driven by EPSR's supplied background-corrected `.q01`
partial and is converted to the same neutron total for scoring. The common
external score is directly comparable, but the programs do not perform
identical internal background handling.

### Speed result

The clean timing run was sequential on a MacBook Pro Mac16,8 with an Apple M4
Pro (10 performance plus 4 efficiency cores), 48 GB RAM, macOS 26.5.2, and
native arm64 binaries. Each process was restricted to one thread. Times are
executable wall time and include program initialization and normal output, but
exclude input preparation and the independent post-run curve scoring.

| Program | Median wall time [minimum, maximum] | Attempted moves | Relative result |
|---|---:|---:|---:|
| Native EPSR26 | 29.09 s [28.13, 29.25] | 1,000,000 | 1.00x |
| rsmith 1.4.0, pre-Verlet baseline | 44.64 s [43.27, 45.14] | 1,000,000 | 1.535x native time |
| rsmith 1.4.0, optimized | **16.71 s [16.52, 16.85]** | 1,000,000 | 0.575x native time |

The optimized rsmith run is 2.670x faster than its frozen pre-optimization
baseline and 1.740x faster than native EPSR26 on this single-thread benchmark.
The optimization does not change the 12 A physical cutoff or any Monte Carlo
acceptance rule. It replaces a 3x3x3 cell list whose stencil visited essentially
all 5,000 atoms per trial with a directed cutoff-plus-skin Verlet list. A 3 A
skin gives a 15 A list radius and approximately 738 candidates per atom. Each
epoch builds the list once; none of the ten runs required an intra-epoch
rebuild.

The ten-seed acceptance, diffraction, and structural distributions are
numerically identical to the frozen pre-Verlet observations. Correctness is
also checked against brute-force periodic energy deltas through accepted moves,
cell crossings, and forced neighbor-list rebuilds. This makes the timing gain
an algorithmic result rather than a shorter-cutoff or altered-sampling result.

The machine-readable observations are frozen in
`expected/stochastic-fixed-budget.toml`. Raw EPSR inputs, program outputs, and
per-seed timing files remain local and ignored because of the workshop-data
terms.

During this comparison, the stochastic path also exposed a scattering-
convention defect: rsmith previously compared internal S(Q), whose baseline is
one, directly with input declared as `convention = "iq"`, whose baseline is
zero. The input is now converted back to internal S(Q) before the EPSR update,
with a round-trip regression test. Without this correction the empirical
potential would contain a spurious constant residual.

The next stages are (1) add independent physical-structure validation to this
Ga comparison, then (2) repeat the protocol for the 6,000-atom `DTBsilicaNX`
neutron-plus-X-ray example.
