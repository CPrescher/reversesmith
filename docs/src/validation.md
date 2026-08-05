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

### Matched-fit hidden-structure recovery

Two 3,000-atom ambient-silica models provide symmetric synthetic recovery
tests: the Pedone structure starts the GAP-target refinement, and the GAP
structure starts the Pedone-target refinement. Native EPSR26 and rsmith are
compared only at common externally measured neutron/X-ray fit quality; the
hidden partial g(r) curves are never used for refinement.

In a ten-seed dense-checkpoint ensemble, rsmith has lower hidden-RDF error in
9/10 GAP-target seeds and 10/10 Pedone-target seeds. Mean native-minus-rsmith
errors are 0.02954 (95% paired-bootstrap interval 0.01905--0.03863) and
0.22145 (0.21572--0.22622), respectively.

A separate five-seed control matrix varies feedback, reference-potential
strength, and two independently reference-equilibrated common starts. Rsmith
is strong under the preregistered rule in eleven of twelve target-control
combinations. It wins all twenty common-start comparisons. The exception is
GAP recovery at 1.5 times reference strength: rsmith wins 2/5 seeds and the
mean difference is 0.00034 with interval -0.00999--0.01071. Consequently the
all-controls superiority claim is not retained. The validated conclusion is
that rsmith usually recovers these hidden ambient-silica structures better at
matched scattering fit and is robust to starting state, but its advantage is
sensitive to reference-potential strength. This result motivates a focused
high-pressure discriminator rather than additional ambient sensitivity grids.

### Synthetic 0-to-10 GPa silica recovery

The high-pressure discriminator uses a 3,000-atom Erhard-2024 ACE trajectory
from the `SiO2_glass` repository. The ambient structure is mapped fractionally
into the exact 10 GPa box, and synthetic neutron/X-ray totals from the hidden
10 GPa endpoint drive refinement. Hidden partial RDF and coordination data are
used only by the common scorer. The same ACE generates the hidden target and
regularizes HRMC, making this an intentionally favorable inverse-oracle test.

The stock-EPSR hidden-coordinate replay agrees at `6.48e-9` neutron and
`1.94e-8` X-ray RMS. Across five nominally matched 30,000-move endpoints,
weight-30 rsmith PACE has lower hidden partial-RDF error in 5/5 paired seeds.
Its median hidden-RDF error is 0.43715 versus 0.94698 for EPSR, and its common
neutron/X-ray RMS is 0.06758 versus 0.14807. Thus PACE recovers 59.5% of the
starting structural gap while EPSR recovers 12.4%.

The preregistered overall-superiority gate nevertheless fails. All five PACE
endpoints violate at least one local-distance floor, predominantly O-O minima
of 1.871--1.966 A; EPSR passes four of five endpoints. Pure RMC fits still
better but is much more pathological. PACE's total ACE energy decreases in all
five runs, demonstrating that a favorable global MLIP energy is not itself a
hard local-safety guarantee. Local endpoint timings also favor EPSR by about
37.6 times over native PACE, but differing move semantics make this diagnostic
rather than a production speed claim. The defensible conclusion is superior
hidden-structure recovery on this favorable synthetic task, not general
superiority over EPSR. The next algorithmic requirement is to prevent rare
local overlaps without discarding the MLIP-guided recovery advantage.

The follow-up separates the loading path into 0-to-5 and 5-to-10 GPa steps and
replaces ambient hard cutoffs with pressure-relative 0.1% and 1% lower-tail
distance quantiles. At 30,000 nominal moves, PACE recovers 68.6% and 71.3% of
the two hidden-RDF gaps, compared with 12.8% and 12.6% for EPSR. Its hidden-RDF
errors are 63.9% and 67.2% lower than EPSR, so the preregistered recovery gate
passes in both increments.

Sequential refinement also removes the severe contacts from the direct
0-to-10 jump: PACE minima remain at or above 2.109 A for Si-Si, 1.461 A for
Si-O, and 2.117 A for O-O. The target-relative lower-tail result is mixed,
however. PACE better reproduces Si-O while EPSR better reproduces O-O, so the
equal-pair mean lower-tail gate still favors EPSR. This is appropriately
reported as a distribution-level structural discrepancy, not an absolute
high-pressure safety failure.

The completed extension uses every 10 GPa loading endpoint through 70 GPa and
adds the installed native RMCProfile 6.7.9.5 executable. All four arms receive
synthetic neutron and X-ray reciprocal-space totals; partial RDFs and
coordination distributions remain hidden. EPSR26 and RMCProfile each generate
and replay their own native hidden-coordinate targets. The maximum EPSR replay
errors are `1.38e-8` neutron and `4.66e-8` X-ray RMS, and RMCProfile replay is
exact at stored precision. RMCProfile uses zero hard minimum distances so no
ambient cutoff is imposed as high-pressure truth.

| Target | Start RDF | rsmith RMC | PACE w30 | EPSR26 | RMCProfile |
|---:|---:|---:|---:|---:|---:|
| 10 GPa | 1.17607 | **0.44012** | 0.47914 | 1.03925 | 0.77136 |
| 20 GPa | 0.81465 | **0.24057** | 0.27772 | 0.69210 | 0.48452 |
| 30 GPa | 0.48126 | **0.16091** | 0.21673 | 0.45373 | 0.30509 |
| 40 GPa | 0.25542 | **0.11161** | 0.14761 | 0.27896 | 0.20233 |
| 50 GPa | 0.12928 | 0.09161 | **0.09061** | 0.18702 | 0.15194 |
| 60 GPa | 0.11140 | **0.08360** | 0.08682 | 0.17201 | 0.14273 |
| 70 GPa | 0.08772 | 0.07826 | **0.07672** | 0.15531 | 0.14108 |

PACE has lower hidden partial-RDF, common neutron/X-ray `i(Q)`, and
Si-coordination errors than both native programs in all seven single-seed
steps. Its recovered hidden-RDF fraction decreases from 59.3% at 10 GPa to
12.5% at 70 GPa as the mapped starting gap itself becomes smaller. EPSR26
worsens the mapped-start hidden RDF from the 30-to-40 GPa step onward;
RMCProfile does so from 40-to-50 GPa onward. These statements concern the
common independent coordinate scorer, not each program's differently
normalized native objective.

Pure RMC has the best hidden-RDF score in five of seven steps but produces
sub-angstrom Si-O contacts and loses to PACE on the pressure-relative
lower-tail score in every step. PACE lowers the ACE energy throughout and
keeps its observed minima within 2.050--2.278 A Si-Si, 1.371--1.472 A Si-O,
and 1.929--2.003 A O-O. Its lower-tail score is worse than EPSR through 30 GPa
but better at 40--70 GPa, while it beats RMCProfile at every pressure.

Median local endpoint walls are 1.81 s pure rsmith RMC, 2.84 s EPSR26, 77.45 s
RMCProfile, and 208.42 s PACE. They are adapter diagnostics rather than a
production-speed ranking because move semantics and native objectives differ.
The defensible conclusion is that MLIP-regularized rsmith provides the best
combination of held-out recovery and local structure in this favorable
same-ACE synthetic trend scan. It is not yet a universal superiority or
experimental high-pressure claim.

The next preregistered extension will compare rsmith `S(Q)`-only, local
`G(r)`-only, and joint `S(Q)+G(r)` refinement. The real-space X-ray total will
be transformed from the same input at `Qmax = 17 inverse Angstrom`; transform
window, local `r` range, uncertainties, and domain weights must be fixed before
outcomes because the two residuals are complementary weightings of the same
scattering information, not independent experiments.

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

### Ten-seed fixed-budget silica comparison

The completed ensemble contains 120 independently scored endpoints: six
methods, two symmetric 3,000-atom cross-recovery cases, and ten seeds after
6,000 attempted moves. RMCProfile now uses the intended neutron uncertainty,
an information-matched scalar X-ray uncertainty, explicit seeds, and audited
exact 6,000-move configurations. All programs use one thread.

| Hidden target / method | Combined fit RMS, median | Hidden partial-RDF RMS, median | Wall time, median |
|---|---:|---:|---:|
| GAP / native EPSR26 | 0.11489 | 0.76711 | 1.045 s |
| GAP / RMCProfile | 0.11422 | 0.72140 | 6.606 s |
| GAP / rsmith EPSR | 0.11177 | 0.71034 | 0.390 s |
| GAP / PACE weight 3 | 0.09747 | 0.59747 | 12.986 s |
| GAP / Pedone weight 3 | 0.09750 | 0.59841 | 0.675 s |
| GAP / pure RMC | 0.09757 | 0.60507 | 0.510 s |
| Pedone / native EPSR26 | 0.11796 | 0.76684 | 1.075 s |
| Pedone / RMCProfile | 0.11576 | 0.59904 | 6.404 s |
| Pedone / rsmith EPSR | 0.11559 | 0.72070 | 0.407 s |
| Pedone / PACE weight 3 | 0.09378 | 0.45780 | 15.690 s |
| Pedone / Pedone weight 3 | 0.09404 | 0.45839 | 0.668 s |
| Pedone / pure RMC | 0.09431 | 0.46323 | 0.515 s |

At identical seeds PACE lowers mean partial-RDF RMS relative to pure RMC by
`0.00761` and `0.00802`; the paired 95% bootstrap intervals exclude zero in
both directions. Its difference from Pedone weight 3 is unresolved in both
directions. This supports a modest benefit from energy-regularized hybrid RMC,
but not a PACE-specific benefit over Pedone.

The fixed move count is not a full EPSR comparison. Native EPSR26 resets the
empirical potential and performs one potential-refinement iteration; rsmith
also uses one iteration, so its first updated empirical potential is produced
only after the scored MC endpoint. The imported finished tutorial records 495
accumulated iterations. These EPSR endpoints therefore test initialization and
short-epoch throughput, not the converged method. This explains why the direct
data-fitting hybrid-RMC arms reach a different residual range after only two
attempted moves per atom.

No PACE endpoint overlaps native EPSR26 or RMCProfile within the preregistered
`0.002` achieved-fit tolerance. The ensemble consequently does not establish
rsmith superiority over either external code. Where support does overlap,
RMCProfile recovers a lower hidden RDF error than native EPSR for the GAP
target and than rsmith EPSR for the Pedone target. A proper EPSR test must run
the potential-refinement schedule toward convergence, save trajectory
checkpoints, and compare methods only at common achieved residuals.

### Iterative EPSR pilot

The follow-up uses five attempted moves per atom between empirical-potential
updates and independent seeded prefixes through iterations 1, 2, 5, 10, 25,
50, and 100. Removing rsmith's explicit minimum-distance constraints gives
byte-identical coordinates through iteration 50, so constraint asymmetry does
not explain the result.

A native EPSR control set pair `rminex` values to the same nominal Si-Si,
Si-O, and O-O minima. EPSR retained those inputs, but every checkpoint for
both targets was byte-identical to the default trajectory and still crossed
the requested Si-Si threshold. Therefore `rminex` is not an effective hard
move-rejection control in this EPSR26 path. The parity comparison uses stock
EPSR26 and unconstrained rsmith; adding a hard rejection to EPSR source would
define a separate sensitivity method.

Four native/rsmith checkpoint pairs fall within the frozen `0.002` combined-
fit tolerance. In all four, rsmith has lower hidden partial-RDF RMS and reaches
the endpoint 3.12--3.95 times faster. For example, toward GAP native iteration
10 gives fit/RDF `0.09597/0.41730` in 5.20 s, while rsmith iteration 5 gives
`0.09661/0.32782` in 1.64 s. Toward Pedone the corresponding values are
`0.09598/0.55230` and `0.09591/0.40366` in 5.55 and 1.74 s.

At iteration 100, native EPSR26 has Si-Si minima of 1.925 A toward GAP and
1.579 A toward Pedone, well below the hidden-target values of 2.257 and 2.445
A. Unguarded rsmith remains at 2.706 and 2.311 A and has lower common fit and
hidden-RDF errors. The planned 250/500 extension is stopped on structural-
stability grounds. This is promising evidence for rsmith and a reproducible
failure of this native setup, but it is a single-seed synthetic pilot. A paper
claim requires seed replication, denser matched-fit checkpoints, sensitivity
to EPSR controls and starting state, and held-out or experimental validation.

### Dense matched-fit EPSR ensemble

The preregistered follow-up completes ten seeds for both cross-recovery cases,
with 260 native EPSR26 and 280 unconstrained-rsmith endpoints. A monotone
one-to-one match gives 12 fit-matched pairs per GAP seed and 11 per Pedone seed;
every primary pair is within the frozen `0.002` combined-fit tolerance.

| Target | Rsmith lower hidden RDF | Mean native - rsmith RDF (95% CI) | Native Si-Si < 2.0 A at iteration 100 |
|---|---:|---:|---:|
| GAP from Pedone | 9/10 | 0.02954 (0.01905--0.03863) | 4/10 |
| Pedone from GAP | 10/10 | 0.22145 (0.21572--0.22622) | 10/10 |

At the best-fit and worst-fit ends of common support, rsmith has lower hidden-
RDF error in all ten seeds for both cases, with all four paired-bootstrap
intervals excluding zero. The hidden partial-S(Q) difference is unresolved for
GAP but favors rsmith for Pedone. Independently resampling native and rsmith
RDF values also leaves positive confidence intervals for GAP
(`0.02093--0.03785`) and Pedone (`0.21664--0.22590`), so the result does not
depend on pairing unlike RNG streams by seed label. Concurrent ensemble wall
times are excluded from speed claims. This establishes a reproducible advantage over stock
EPSR26 for hidden-structure recovery in these synthetic tests; experimental
validation and control/start sensitivity remain necessary before generalizing.

## Remaining publication gates

1. Validate liquid-Ga structures against held-out data or an independent
   atomistic/physical oracle.
2. Test independent equilibrium starts to separate basin sensitivity from
   convergence out of unstructured liquids.
3. Test the replicated native/rsmith result against reasonable EPSR feedback,
   reference-potential, and equilibrium-start sensitivity. The ten-seed dense
   matched-fit replication is complete; start/control generality is not.
4. Complete the separately preregistered cold/warm timing-only campaign. The
   ten fresh-process scientific timings are repeatable warm-cache diagnostics,
   not final publication speed ratios.
5. Apply the matched-fit comparison to a complex high-pressure glass, where
   chemically diagnostic coordination, angle, and topology observables can
   establish or reject a physical-superiority claim.
