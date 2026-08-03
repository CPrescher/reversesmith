# Validation and Reference-Code Benchmarks

rsmith uses reference-code benchmarks as scientific gates, not only as speed
tests. Deterministic forward calculations must agree before stochastic
refinement results are compared. Refinements use repeated seeds, a common
external score, structural diagnostics, and explicit reporting of censored or
failed runs.

The complete machine-readable programme is maintained in the repository's
[benchmarks directory](https://github.com/CPrescher/rsmith/tree/main/benchmarks).
Raw upstream inputs are committed only when their licences permit
redistribution.

## EPSR26 LiquidGa50C

The EPSR26 liquid-Ga workshop example currently provides the most complete
external comparison.

### Forward-calculation gate

Starting from the identical 5,000-atom configuration, rsmith reproduces the
native EPSR partial neutron interference function to 0.0381% RMS of its dynamic
range. The 0.12 A-rebinned RDF differs by 0.752% RMS/range, while an independent
minimum-image histogram agrees with rsmith to `9.41e-8` RMS/range.

### Supplied-endpoint re-equilibration

Across ten one-thread, one-million-move replicates from the supplied EPSR
endpoint, rsmith takes a median 16.71 s versus 29.09 s for native EPSR26. Its
external measured-total residual is 1.350% versus 2.128%. Native EPSR remains
closer to the supplied long-run EPSR g(r), so this result establishes
competitive refinement and speed rather than greater physical accuracy.

### Independent-start convergence

The supplied endpoint already fits the measured total well, so a separate
protocol uses ten independently generated periodic hard-sphere liquids. Each
native/rsmith pair receives the identical starting coordinates. Time to target
uses the common externally evaluated residual and requires the target to remain
satisfied at every later frozen checkpoint.

| Sustained residual | Native EPSR26 | rsmith | Comparison |
|---:|---:|---:|---:|
| 4.0% | 0.888 s | 0.752 s | 1.18x faster |
| 3.0% | 22.40 s | 2.02 s | **11.11x faster** |
| 2.5% | 22.58 s | 3.72 s | **6.07x faster** |
| 2.0% | 0/10 reached | 10/10 in 7.11 s | rsmith only |
| 1.5% | 0/10 reached | 10/10 in 13.88 s | rsmith only |

At one million moves, the median residuals are 2.390% for native EPSR26 and
1.369% for rsmith. One native seed reproducibly develops a 0.786 A pair and a
large energy excursion; the repeat is bit-identical. This is reported as one
deterministic instability under the frozen protocol, not as a claim that every
EPSR run fails.

The result supports superior time-to-data-agreement for this protocol. It does
not yet prove superior physical structure: the programs do not optimize
identical internal background representations, the starts are deliberately
unstructured rather than equilibrium samples, and the EPSR worked g(r) is not
an independent ground truth.

## Reproduction

After accepting the upstream local-testing terms and importing a private EPSR26
copy, run:

```bash
cd benchmarks/epsr-ga
python3 scripts/run_independent_starts.py --jobs 1
python3 scripts/verify_independent_starts.py
```

The frozen result, per-seed convergence histories, provenance hashes, and
plot-ready median curve are committed under `benchmarks/epsr-ga/expected/`.

## Multi-contrast EPSR readiness

The empirical-potential update now consumes every configured X-ray and neutron
S(Q) dataset rather than selecting one preferred dataset. Each contrast uses
its own Q-dependent or neutron scattering weights, including isotope-specific
overrides, and the projected updates are combined using the dataset weight and
pointwise uncertainty over the active fit range.

This infrastructure is covered by a synthetic conflicting-contrast regression
test: both an X-ray-sensitive and a neutron-sensitive partial change sign as
required, and increasing the neutron precision shifts the combined update
toward the neutron solution. An executable smoke test also exercises a joint
X-ray/neutron pure-EPSR epoch. The historical single-neutron LiquidGa path was
checked against the preceding release binary: the refined coordinates,
empirical potential, and calculated neutron S(Q) are byte-identical. These are
implementation gates; they do not replace the native EPSR silica comparison.

### EPSR26 DTBsilicaNX forward calculation

The exact 6,000-atom worked-example configuration has now passed a deterministic
neutron-plus-X-ray forward gate. Relative to a fresh one-configuration EPSR26
calculation, the Si-Si, Si-O, and O-O partial i(Q) curves differ by 0.0397--
0.0558% RMS of their respective dynamic ranges. After applying EPSR's documented
worked-file normalizations (`nrtype=5` for neutron and single-atom scattering
for X-ray), the neutron and X-ray total differences are 0.0639% and 0.0661%.
The 0.12 A-rebinned partial RDF differences are 0.626--0.934% RMS/range, while
an independent minimum-image histogram agrees with rsmith at better than
`6e-8` RMS/range.

The raw rsmith total files use Faber-Ziman normalization. Comparing those
directly with EPSR's differently normalized totals would create an artificial
31% neutron discrepancy; the benchmark therefore reconstructs the EPSR
normalizations from its weight files before scoring. Licensed upstream files
remain local and ignored, with committed hashes providing provenance.

### EPSR26 DTBsilicaNX reference potential

The second deterministic gate independently reconstructs the reference
potential from the local tutorial's species and control records. It reproduces
the native Si-Si, Si-O, and O-O potential curves to `3.25e-8`--`7.00e-8` RMS
of their dynamic ranges and the exact-configuration energy to 0.105%. The
implementation includes EPSR's geometric epsilon and arithmetic sigma mixing,
modified 12--6 form, charge-product electrostatics, and its distinct
short-range and Coulomb cutoff functions.

A one-thread, 6,000-move smoke trajectory then exercised the verified
reference potential together with both measured neutron and X-ray contrasts.
All three empirical pair potentials received nonzero updates. This closes the
reference-potential plumbing gate, but the one-epoch wall time and residuals
are deliberately not presented as performance or convergence results. The
next scientific gate is the preregistered matched multi-seed, multi-epoch
native/rsmith reproduction.

## Remaining publication gates

1. Validate liquid-Ga structures against held-out data or an independent
   atomistic/physical oracle.
2. Test independent equilibrium starts to separate basin sensitivity from
   convergence out of unstructured liquids.
3. Run the matched multi-seed native/rsmith silica ensembles, estimate native
   seed-to-seed spread, and freeze the stochastic equivalence margins before
   comparing pair-potential and MLIP-regularized refinements.
4. Apply the preregistered comparison to a complex high-pressure glass, where
   chemically diagnostic coordination, angle, and topology observables can
   establish or reject a physical-superiority claim.
