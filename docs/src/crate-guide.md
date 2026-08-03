# Reverse Monte Carlo refinement with `rsmith`

`rsmith` refines periodic atomistic models against X-ray and neutron total-scattering
data. It is aimed at disordered systems—glasses, liquids, and amorphous solids—for
which a crystallographic unit-cell model is usually not an adequate description.

The command-line program is the usual entry point. The public Rust modules expose
the same structure readers, scattering calculations, constraints, potentials, and
RMC engine for custom workflows.

## The problem RMC solves

A total-scattering experiment measures a one-dimensional function such as the
structure factor \(S(Q)\) or a real-space transform such as \(g(r)\). The structure
we want is three-dimensional and may contain thousands of atomic coordinates.
For a multicomponent material the measured total function is also a weighted sum
of several partial pair correlations.

Reverse Monte Carlo (RMC) approaches this inverse problem directly:

1. Start from a periodic configuration with the desired composition and density.
2. Displace one randomly selected atom.
3. Reject the move if it violates a hard structural constraint.
4. Update the calculated scattering signal and compare it with experiment.
5. Accept improvements and sometimes accept worse moves according to a Metropolis
   probability.
6. Repeat until the fit and the structural diagnostics are satisfactory.

For data sets \(d\), the data-misfit term is a weighted chi-squared,

```text
chi2(X) = sum_d w_d sum_i ((y_calc,d(X, i) - y_exp,d(i)) / sigma_d(i))^2
```

where `X` is the atomic configuration. `rsmith` can fit several compatible data
sets simultaneously and can estimate pointwise uncertainties from the local data
noise when `sigma` is omitted.

RMC is therefore not molecular dynamics run backwards. It does not reconstruct a
trajectory, determine a unique structure, or by itself prove that every fitted
local environment is chemically realistic. It samples configurations that agree
with the supplied data and assumptions.

## Why the inverse problem needs regularization

Many different configurations can produce nearly the same total-scattering
signal. The finite experimental Q range, noise, the loss of directional
information in an isotropic measurement, and the mixing of partial correlations
all make the problem underdetermined. A pure data fit can consequently improve
\(S(Q)\) while creating implausible short contacts, coordination defects, or
distorted bond-angle distributions.

`rsmith` offers progressively richer ways to restrict the solution space:

| Method | Meaning | Best use | Main limitation |
|---|---|---|---|
| Minimum-distance constraint | An absolute forbidden region | Prevent atomic overlap cheaply | Gives no preference among allowed structures |
| Coordination constraint | A hard bound on neighbor counts | Preserve well-established local chemistry | Can encode the desired answer too rigidly |
| Pair-potential regularizer | A smooth energetic preference based on pair distances | Cheap hybrid RMC with a known force field | Pair terms cannot represent all angular or many-body chemistry |
| ML-potential regularizer | A learned many-body energy prior | Preserve richer local environments when a validated model exists | More expensive and only trustworthy inside the model's domain |
| EPSR | Iteratively learns an empirical pair correction from the scattering residual | Move the equilibrium ensemble toward experiment | A distinct iterative method, not a fixed transferable ML model |

Hard constraints answer “is this configuration allowed?” An energy regularizer
instead answers “among the allowed configurations, which ones are more plausible?”
They are complementary.

### Potentials as regulators

With a fixed pair or machine-learning potential, a trial move is judged by

```text
delta_cost = delta_chi2 + lambda * delta_E
```

and an uphill move is accepted with the configured Metropolis probability. The
experimental term pulls the model toward the measurement; the energy term resists
configurations that the potential considers unfavorable. The weight `lambda`
controls the trade-off.

This is most useful when the starting structure came from molecular dynamics or
another atomistic simulation. Using the same potential during RMC lets the
experimental data adjust that model without freely leaving its familiar energy
landscape. It is a guardrail, not a guarantee: a poor potential, an out-of-domain
ML model, or an excessive weight can bias the result away from the real material.

Pair potentials are attractive because the energy change for a one-atom move only
requires its neighbors. `rsmith` supports Lennard-Jones, Buckingham, Pedone, damped shifted-force
Coulomb, and tabulated pair terms. Analytical contributions for the same species
pair can be added, while a tabulated entry replaces them for that pair.

Machine-learning interatomic potentials can encode angular and higher-body
environment information that a simple pair potential cannot. `rsmith` supports:

- native linear or quadratic SNAP models from FitSNAP/LAMMPS model files;
- native PACE product-basis energies from pacemaker C-tilde `.yace` files;
- GAP models through an optional QUIP bridge; and
- MACE checkpoints through a persistent Python worker, with full, local, and
  incremental energy-delta strategies for the models that support them.

An ML potential must cover every element and the thermodynamic/structural regime
visited by the refinement. Validate accelerated local energy deltas against full
energy differences before a production run. Explicit long-range terms, global
charge or spin coupling, and dispersion corrections generally require the full
evaluation path unless they have been separately validated.

### Choosing the regularization strength

There is no universal energy weight because chi-squared depends on the data
uncertainties, number of fitted points, and data-set weights. A practical workflow
is:

1. Run a short data-only or weakly regularized calculation to establish the scale
   of `delta_chi2`.
2. Use the calibration printed by `rsmith`, which compares representative
   `|delta_chi2|` and `|delta_E|` values.
3. Run several weights around the suggested value, usually changing by factors of
   two to five rather than tiny increments.
4. Compare both the scattering residual and independent structural diagnostics.

If the energy rises steadily and local structure deteriorates, increase the
weight. If the scattering fit barely moves, decrease it. The best weight is the
weakest one that prevents clearly unphysical distortion; it is not necessarily
the weight that gives the smallest chi-squared.

## What `rsmith` supports

### Experimental targets

- X-ray total \(S(Q)\), \(i(Q)=S(Q)-1\), or \(F(Q)=Q[S(Q)-1]\)
- neutron total \(S(Q)\), \(i(Q)\), or \(F(Q)\)
- X-ray total \(g(r)\)
- X-ray reduced pair-distribution function \(f(r)\)
- simultaneous fitting of X-ray and neutron reciprocal-space data
- simultaneous reciprocal- and real-space fitting with independent ranges,
  uncertainties, and weights

`g(r)` and `f(r)` are alternative real-space targets and cannot both be selected
in one run. Because a real-space target may be derived from the same \(S(Q)\), do
not give both full independent weight without considering the double-counting.

### Structures and cells

The program reads LAMMPS data files, extended or simple XYZ, and VASP 5+
POSCAR/CONTCAR files. The simulation cell must be orthorhombic. LAMMPS atom types
are mapped to element symbols in the configuration; XYZ and VASP inputs contain
their species names directly. Optional isotropic density rescaling is available.

### Refinement and analysis

- adaptive single-atom moves and exponential simulated annealing;
- optional delayed acceptance that screens proposals against experimental data
  before invoking an expensive potential;
- checkpoint/resume and reproducible random seeds;
- minimum-distance and coordination-number constraints;
- pair-potential, SNAP, GAP/QUIP, or MACE energy regularization;
- hybrid or energy-only EPSR outer iterations;
- coordination-number and bond-angle analysis of starting and refined structures;
- incremental RDF and scattering updates backed by cell lists.

The built-in X-ray form factors and neutron scattering lengths cover the elements
listed in the full guide. Check [`xray`] and [`neutron`] or the
[supported-elements table](https://github.com/CPrescher/rsmith/blob/main/docs/src/elements.md)
before starting an unusual composition.

## Choosing a workflow

Use this as a starting decision guide:

- **Only want the scattering signal of an existing model?** Set
  `[rmc] max_moves = 0` and inspect the generated starting `*_sq.dat` and
  `start_gr.dat` files. The library functions in [`rdf`], [`sq`], [`xray`], and
  [`neutron`] are better for a custom calculation pipeline.
- **Have data but no trustworthy potential?** Start with minimum distances and
  deliberately broad coordination bounds. Run an ensemble of independent RMC
  fits and report the variation.
- **Have an MD structure and its force field?** Use matching pair or tabulated
  potentials for hybrid RMC, plus hard minimum distances.
- **Have a validated ML potential for the material and state point?** Use it as a
  many-body regularizer, first validating energy deltas and runtime on a short
  run.
- **Want an empirical potential whose equilibrium ensemble approaches the data?**
  Use EPSR with a reference pair potential when one is available.

## Tutorial 1: a minimal constrained RMC refinement

Build the command-line program:

```bash
cargo install rsmith
# or, from a source checkout:
cargo build --release
```

Prepare an orthorhombic structure and a two-column experimental file. For
LAMMPS input, a minimal `config.toml` is:

```toml
[system]
structure = "glass.data"
format = "lammps"

[system.types]
1 = "Ca"
2 = "Si"
3 = "O"

[data.xray_sq]
file = "experimental.sq"
weight = 1.0
# sigma = 0.02       # omit to estimate pointwise noise
sigma_alpha = 0.05   # optionally relax the high-Q tail
fit_min = 0.5
fit_max = 18.0

[rmc]
max_moves = 500_000
max_step = 0.1
checkpoint_every = 50_000

[sq]
qmin = 0.3
qmax = 20.0
nq = 500
rdf_cutoff = 11.0
rdf_nbins = 1000
lorch = true

[constraints.min_distance]
"Si-O" = 1.2
"O-O" = 2.0
"Ca-O" = 1.8
"Si-Si" = 2.3
"Ca-Si" = 2.2
"Ca-Ca" = 2.3

[[constraints.coordination]]
pair = "Si-O"
min = 3
max = 6
cutoff = 2.2

[analysis]
angle_triplets = ["O-Si-O", "Si-O-Si"]

[analysis.cutoffs]
"Si-O" = 2.2
"Ca-O" = 3.0
```

Paths are resolved relative to the TOML file. Run the refinement and then compare
the local structure before and after it:

```bash
rsmith config.toml --output-dir run01
rsmith config.toml --analyze --output-dir run01
```

Inspect at least:

- `refined.xyz` and the final calculated scattering files;
- the fitted residual over the requested range, not just the aggregate chi-squared;
- acceptance rate and step-size history in `rsmith.log`; and
- starting versus refined coordination and bond-angle distributions.

Do not tune constraints from the final RMC structure alone. Base them on chemical
knowledge, independent measurements, or a validated simulation, and record them as
part of the inference assumptions.

## Tutorial 2: add a classical-potential guardrail

Suppose the starting model was equilibrated with a Pedone short-range potential
plus damped shifted-force Coulomb interactions. Add the same parameters to the
previous configuration and ensure `rdf_cutoff` is at least the potential cutoff:

```toml
[sq]
qmin = 0.3
qmax = 20.0
nq = 500
rdf_cutoff = 15.0
rdf_nbins = 1500
lorch = true

[potential]
weight = 0.001
cutoff = 15.0

[[potential.pedone]]
pair = "Ca-O"
D0 = 0.030211
alpha = 2.241334
r0 = 2.923245
C0 = 5.0

[[potential.pedone]]
pair = "Si-O"
D0 = 0.340554
alpha = 2.006700
r0 = 2.100000
C0 = 1.0

[[potential.pedone]]
pair = "O-O"
D0 = 0.042395
alpha = 1.379316
r0 = 3.618701
C0 = 22.0

[potential.coulomb]
alpha = 0.25
charges = { Ca = 1.2, Si = 2.4, O = -1.2 }
```

The values above are an example for one model, not universal Ca-Si-O parameters.
Use the exact units, charges, cutoff, and parameters from the simulation that
created the starting structure. For an unsupported analytical form, export the
complete pair energy as a two-column table and configure it with
`[[potential.tabulated]]`.

Run a short calibration, then compare an ensemble at a few energy weights. A good
regularizer will usually allow some degradation in the numerical data fit in
exchange for a much more credible local structure.

## Tutorial 3: use a machine-learning potential

An ML potential is configured instead of `[potential]`; the two cannot be combined
in the current RMC implementation. It also cannot currently be combined with
`[epsr]`.

Native SNAP needs no external runtime:

```toml
[ml_potential]
backend = "snap_native"
coefficient_file = "potential.snapcoeff"
parameter_file = "potential.snapparam"
weight = 0.001
```

The cutoff is derived from the SNAP files and must not be set manually. The cell
must be larger than twice the largest model cutoff in each direction.

Native PACE likewise needs no external runtime:

```toml
[ml_potential]
backend = "pace_native"
model = "potential.yace"
weight = 0.001
```

The cutoff is derived from the directed bonds in the C-tilde model. Do not set
`cutoff` or `delta` for this backend.

MACE uses a Python environment containing MACE, PyTorch, and ASE:

```toml
[ml_potential]
backend = "mace_python"
model = "mace.model"
python = ".venv/bin/python"
weight = 0.001
cutoff = 5.0
delta = "full"
device = "cpu"
torch_threads = 8
```

Start with `delta = "full"` as the reference. For an ordinary finite-range MACE
checkpoint, compare representative trial deltas against `local` or `incremental`
mode before switching a production run. There is intentionally no silent fallback
when an incremental checkpoint feature is unsupported.

GAP uses the optional `gap-quip` Cargo feature and an external QUIP installation:

```toml
[ml_potential]
backend = "gap_quip"
model = "gap.xml"
init_args = "Potential xml_label=GAP_LABEL"
weight = 0.001
cutoff = 5.0
```

For GAP and MACE, the configured cutoff must cover the model's actual environment
range. A value that is too small can produce an incorrect local energy change
without making the scattering calculation itself fail.

See the [complete ML-potential configuration guide](https://github.com/CPrescher/rsmith/blob/main/docs/src/config/ml-potential.md)
for backend setup, compatibility, performance, and validation details.

## Tutorial 4: run and compare an ensemble

The inverse problem remains non-unique after regularization. Independent random
trajectories reveal whether a reported feature is stable or merely one local
solution. Omit `seed` from `[rmc]` so each run receives a random seed; every seed is
written to its log.

```bash
for i in 01 02 03 04 05 06 07 08; do
  rsmith config.toml --output-dir "run${i}" --quiet &
done
wait

for i in 01 02 03 04 05 06 07 08; do
  rsmith config.toml --analyze "run${i}/refined.xyz" \
    --output-dir "run${i}" --quiet
done
```

Compare the distributions of fit quality, partial pair correlations, coordination
numbers, bond angles, and regularizer energy. Report ensemble variation for
structural conclusions that are not uniquely fixed by the data.

## Reading the results critically

A completed run is not automatically a validated structure. Before interpreting
it, ask:

1. Does the calculated signal agree over the full scientifically relevant range,
   without fitting obvious noise or normalization errors?
2. Are the density, composition, scattering convention, and Fourier-transform
   settings consistent with the experiment?
3. Did coordination and bond-angle distributions remain plausible?
4. Does the conclusion survive changes in random seed, starting configuration,
   data weights, and reasonable regularization strength?
5. Is the potential valid for the elements, bonding environments, pressure, and
   temperature represented in the run?
6. Which structural features are actually constrained by the measured weighted
   totals, and which are inherited mainly from the prior?

The most defensible use of `rsmith` is therefore model refinement and hypothesis
testing: combine experimental information with explicit, documented physical
assumptions, then test how sensitive the conclusions are to those assumptions.

## Command-line reference

```text
rsmith <config.toml> [OPTIONS]

--analyze [structure]  Analyze the configured start and refined structures,
                       or only the explicitly supplied structure
--output-dir DIR       Write all outputs to a separate directory
--seed N               Override the random seed
--quiet                Write progress only to the log
--resume               Continue from checkpoint.dat in the output directory
--density RHO          Isotropically rescale to a mass density in g/cm^3
```

For all TOML fields and output files, continue with the
[full user guide](https://github.com/CPrescher/rsmith/blob/main/docs/src/SUMMARY.md),
the [complete configuration example](https://github.com/CPrescher/rsmith/blob/main/docs/src/full-example.md),
and the module documentation below. The main library entry points are [`config`],
[`rmc`], [`potential`], [`ml_potential`], [`constraints`], and [`analyze`].

## References

- R. L. McGreevy and L. Pusztai, “Reverse Monte Carlo Simulation: A New
  Technique for the Determination of Disordered Structures,” *Molecular
  Simulation* **1**, 359–367 (1988),
  [doi:10.1080/08927028808080958](https://doi.org/10.1080/08927028808080958).
- A. K. Soper, “Empirical potential Monte Carlo simulation of fluid structure,”
  *Chemical Physics* **202**, 295–306 (1996),
  [doi:10.1016/0301-0104(95)00357-6](https://doi.org/10.1016/0301-0104(95)00357-6).
- A. P. Bartók, M. C. Payne, R. Kondor, and G. Csányi, “Gaussian Approximation
  Potentials: The Accuracy of Quantum Mechanics, without the Electrons,”
  *Physical Review Letters* **104**, 136403 (2010),
  [doi:10.1103/PhysRevLett.104.136403](https://doi.org/10.1103/PhysRevLett.104.136403).
- I. Batatia *et al.*, “MACE: Higher Order Equivariant Message Passing Neural
  Networks for Fast and Accurate Force Fields,” *NeurIPS* (2022),
  [arXiv:2206.07697](https://arxiv.org/abs/2206.07697).
