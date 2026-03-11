# Full Example

A complete configuration file for CaSiO3 glass refinement with Pedone + Coulomb potentials, fitting both S(Q) and g(r).

```toml
# =============================================================
# System: CaSiO3 glass, 10000 atoms from MD melt-quench
# =============================================================

[system]
structure = "CaSiO3_glass.data"
format = "lammps"

[system.types]
1 = "Ca"
2 = "Si"
3 = "O"

# =============================================================
# Experimental data
# =============================================================

[data.xray_sq]
file = "CaSiO3_ambient_sample.sq"
weight = 1.0
sigma = 0.02
fit_min = 0.5
fit_max = 18.0

[data.xray_gr]
file = "CaSiO3_ambient_sample.gr"
weight = 0.3
sigma = 0.02
fit_min = 1.0
fit_max = 7.0
qmax = 17.97
lorch = true

# =============================================================
# RMC parameters
# =============================================================

[rmc]
max_moves = 1_000_000
max_step = 0.1
checkpoint_every = 50_000
seed = 42
print_every = 1000
target_acceptance = 0.3
adjust_step_every = 5000

# Simulated annealing
anneal_start = 2.0
anneal_end = 0.01
anneal_steps = 200_000

# Early stopping
convergence_threshold = 1e-4
convergence_window = 50_000

# =============================================================
# S(Q) computation
# =============================================================

[sq]
qmin = 0.0
qmax = 20.0
nq = 500
lorch = true
rdf_cutoff = 20.0     # Must be >= potential cutoff; 20 A for accurate FSDP
rdf_nbins = 1000      # dr = 0.02 A

# =============================================================
# Hard constraints
# =============================================================

[constraints.min_distance]
"Ca-O" = 1.8
"Si-O" = 1.2
"O-O" = 2.0
"Ca-Si" = 2.2
"Ca-Ca" = 2.3
"Si-Si" = 2.3

[[constraints.coordination]]
pair = "Si-O"
min = 3
max = 6
cutoff = 2.2

# =============================================================
# Pair potentials (Pedone 2006 + Coulomb DSF)
# Matches LAMMPS: hybrid/overlay pedone 15.0 coul/dsf 0.25 15.0
# =============================================================

[potential]
weight = 0.001          # Start here, adjust based on calibration output
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

# =============================================================
# Analysis settings (for --analyze mode)
# =============================================================

[analysis]
angle_triplets = ["O-Si-O", "Si-O-Si", "O-Ca-O", "Ca-O-Si"]

[analysis.cutoffs]
"Si-O" = 2.2
"Ca-O" = 3.0
```

## Workflow

```bash
# 1. Check S(Q) quality -- compare initial S(Q) with experiment
#    If the first peak (FSDP) is underestimated, increase rdf_cutoff
reversesmith config.toml --compute-sq-only

# 2. Run refinement
reversesmith config.toml

# 3. Validate the result
reversesmith config.toml --analyze
```

## Expected output

During refinement:

```
Loading structure from "CaSiO3_glass.data" ...
  10000 atoms, 3 species: ["Ca", "Si", "O"]
  Box: 51.5000 x 51.5000 x 51.5000 A

Pair potentials (weight = 0.001000, cutoff = 15.0 A):
  Ca-O: 15001 bins, dr = 0.0010 A
  Si-O: 15001 bins, dr = 0.0010 A
  O-O: 15001 bins, dr = 0.0010 A
  Ca-Ca: 15001 bins, dr = 0.0010 A
  ...

Initial potential energy = -85432.123456 eV (weight = 0.001000)

Calibration (1000 moves): avg |delta_chi2| = 0.0023, avg |delta_E| = 0.34 eV
  Current weight = 0.001000, suggested weight for equal balance = 0.006765
  Ratio current/suggested = 0.1478

Move 1000/1000000: cost = -82.1 (chi2: 3.3 [sq: 2.2, gr: 1.1], w*E: -85.4 [E: -85432.12]), ...
```

Use the calibration output to adjust the weight for subsequent runs.
