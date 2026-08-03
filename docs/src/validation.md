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
All three empirical pair potentials received nonzero updates. An audit also
fixed the target reader to use EPSR's `.u01` calculated total rather than its
similarly shaped `.v01` data-minus-model residual. This closes the
reference-potential plumbing gate, but the one-epoch wall time and residuals
are deliberately not presented as performance or convergence results. The
next scientific gate is the preregistered matched multi-seed, multi-epoch
native/rsmith reproduction.

### Ambient GAP and Pedone structure controls

The silica protocol additionally pins two 3,000-atom, 300 K endpoints from the
private `CPrescher/SiO2_glass` repository: a Pedone-pre-equilibrated GAP quench
and a Pedone-only fast quench. Independent analysis records RDFs, coordination,
O-Si-O and Si-O-Si angles, minimum distances, and a shortest-cycle diagnostic
for the Si network.

The GAP endpoint has 0.40% Si and 0.60% O coordination defects at a frozen
2.2 A Si-O cutoff; the Pedone endpoint has 0.40% and 0.85%. Mean Si-O-Si angles
are 139.21 and 148.06 degrees, respectively. These are useful model controls,
but neither is experimental ground truth: each is a single deterministic
endpoint from a `10^14 K/s` quench, and their densities differ substantially
(2.282 versus 2.560 g/cm3). The preregistered claim rule therefore forbids
using agreement with GAP alone as evidence of physical superiority.

### Symmetric cross-program silica smoke

The pinned GAP and Pedone endpoints now form two symmetric synthetic recovery
cases. In each direction the opposite endpoint is isotropically rescaled to
the hidden target box without relaxation or random displacement. rsmith RMC,
Pedone-HRMC, GAP/QUIP-HRMC, rsmith EPSR mode, native EPSR26, and RMCProfile all
execute from the same cross-start; PDFgui/PDFfit2 supplies a forward-only PDF
control because it is not an atom-by-atom RMC engine.

Program-native hidden-coordinate target/replay gates pass before refinement,
and all recovered coordinates are rescored with one independent Faber-Ziman
and minimum-image analyzer. In the 6,000-move adapter smoke every joint method
improves the mean common neutron/X-ray RMS over its cross-start. rsmith pure
RMC reduces the common totals from 0.10880/0.12332 to 0.09122/0.10339 for the
GAP target and from 0.11170/0.12596 to 0.08541/0.10174 for the Pedone target.

For native EPSR26 and RMCProfile, “matched” currently means the same hidden
coordinates, density-rescaled cross-start, neutron/X-ray information, 6,000
attempted atom moves, and independent final-coordinate scorer. It does not yet
mean the same accepted-move count, achieved residual, or convergence schedule.
EPSR26 has one short epoch in this adapter smoke. RMCProfile used the official
macOS wrapper's serial fallback, and its timed checkpoints contain 5,315 and
5,490 generated moves although the logs reach 6,000. These are therefore
correctness and plumbing controls, not a scientific ranking.

This is deliberately not evidence that GAP-HRMC is superior. The smoke's
`0.001` Pedone/GAP energy weight is orders of magnitude below the calibration
diagnostic, so both energy-regularized paths remain almost indistinguishable
from pure RMC. End-to-end rsmith RMC is about 2x faster than native EPSR26 and
13x faster than the installed serial-fallback RMCProfile path in this local
smoke, whereas the current GAP/QUIP path takes 322--343 s and is roughly 600x
slower than pure RMC. These timings expose the next optimization target; they
are not publication speed ratios.

### HRMC energy-weight pilot

The next frozen pilot used the same symmetric starts and targets for 1,000
moves while sweeping Pedone weights from `0.001` to `100` and GAP weights from
`0.001` to `3`. Pedone weights 3, 10, and 30 retain at least 63% of the
same-budget pure-RMC progress in both neutron and X-ray fits. At weight 30 the
positive energy drift is 75--79% smaller than at weight `0.001`; weight 100
reduces the energy but fails the preregistered 50% fit-progress guard.

For GAP, weight 0.3 retains 89--94% of pure-RMC progress and reduces positive
energy drift by 6--8%. Weight 1.0 instead crosses an acceptance cliff: only 0
and 10 of 1,000 moves are accepted in the two directions. A finer 0.3--1.0
pilot is required before the matched 6,000-move comparison. These results
calibrate regularizer influence; they do not establish convergence or a
structural advantage over EPSR.

The frozen fine-grid follow-up tested 0.4, 0.5, 0.6, and 0.8. Weight 0.4 is the
largest value that passes the 50% guard in both directions, retaining 79--85%
of pure-RMC neutron/X-ray progress and reducing positive energy drift by
10--12% relative to `0.001`. Weight 0.5 passes from the GAP start toward the
Pedone target but retains only 37--39% progress in the reverse direction. The
selected GAP low/knee/high production bracket is therefore 0.1, 0.3, and 0.4.

### Native PACE silica control

The benchmark also pins the public Si-O ACE model of Erhard et al. (2024), DOI
`10.5281/zenodo.10419194`. The model uses the `SBessel` radial basis and is
evaluated directly by rsmith's native PACE backend. Against the official
`python-ace` implementation, total energies for the 3,000-atom cross-start and
for a nonzero single-atom displacement agree to the six decimals printed by
the rsmith CLI. The pinned source hashes, software commits, configurations,
and numeric oracle are recorded in `expected/pace2024-oracle.toml`.

A preregistered 1,000-move scan selected PACE weights 3, 10, and 30. Weight 30
is the highest point retaining at least 50% of pure-RMC neutron and X-ray
progress in both directions; weight 100 fails the symmetric guard. Native PACE
takes 2.37--3.08 s per 1,000 moves in this control, 20--24 times less wall time
than the matched external GAP/QUIP path.

### Joint-acceptance Pedone/GAP/PACE production smoke

An initial delayed-acceptance optimization was rejected after GAP weights 0.1
and 0.3 accepted no moves in 6,000 attempts. The separate Metropolis stages
lost favorable cancellation between opposing data and energy increments. The
failed gate is retained for audit, and the replacement protocol uses the
intended single test on `delta_chi2 + weight*delta_energy`.

The single-seed replacement runs the frozen Pedone 3/10/30, GAP 0.1/0.3/0.4,
and PACE 3/10/30 brackets for 6,000 moves in both symmetric directions. PACE
weight 3 retains 100--103% of same-seed pure-RMC neutron/X-ray progress and
improves hidden-target mean partial-RDF RMS by 1.0% and 2.6%. Pedone weight 3
gives a comparably small improvement, so no PACE-specific structural advantage
is established. The matched PACE arms take 12.3--16.6 s versus 307.7--354.8 s
for GAP/QUIP, a 21--25 times speedup. This is an implementation, calibration,
and throughput result; repeated timings and multi-seed structural uncertainty
are still required for a paper claim.

GAP/QUIP is now archived at this point rather than carried into the multi-seed
campaign. Its completed single-seed results remain useful for documenting the
external MLIP path and the 21--25 times native-PACE speed difference, but the
next compute budget is assigned to pure RMC, Pedone weight 3, native PACE
weight 3, rsmith EPSR, native EPSR26, and RMCProfile. This decision and the
remaining matching limitations are frozen in
`expected/next-campaign-scope.toml`.

## Remaining publication gates

1. Validate liquid-Ga structures against held-out data or an independent
   atomistic/physical oracle.
2. Test independent equilibrium starts to separate basin sensitivity from
   convergence out of unstructured liquids.
3. Run pure RMC, Pedone weight 3, and native PACE weight 3 across the frozen
   seed set and cold/warm timing replicates; do not add further GAP/QUIP runs.
   Repair the RMCProfile parallel runtime and freeze exact coordinate-save and
   timing rules.
4. Run the matched multi-seed native/rsmith silica ensembles, estimate native
   seed-to-seed spread, and freeze stochastic equivalence margins before
   comparing pair-potential and MLIP-regularized refinements.
5. Apply the preregistered comparison to a complex high-pressure glass, where
   chemically diagnostic coordination, angle, and topology observables can
   establish or reject a physical-superiority claim.
