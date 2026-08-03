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

## Remaining publication gates

1. Validate liquid-Ga structures against held-out data or an independent
   atomistic/physical oracle.
2. Test independent equilibrium starts to separate basin sensitivity from
   convergence out of unstructured liquids.
3. Reproduce the EPSR silica neutron-plus-X-ray example before comparing
   pair-potential and MLIP-regularized rsmith refinements.
4. Apply the preregistered comparison to a complex high-pressure glass, where
   chemically diagnostic coordination, angle, and topology observables can
   establish or reject a physical-superiority claim.
