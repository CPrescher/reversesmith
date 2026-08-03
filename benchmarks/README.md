# rsmith scientific benchmark programme

This directory tracks the validation and comparison work required for a
publication-quality description of `rsmith`.  The governing rule is to reuse
the examples and success criteria chosen by the authors of the reference
programs wherever their inputs can be obtained and redistributed.

Benchmark results are not accepted from a single selected stochastic run.
Forward calculations use fixed numerical tolerances; refinements use repeated
seeds and compare distributions of observables.

## Benchmark ladder

| ID | Reference | Published purpose | rsmith reproduction target | Current status |
|---|---|---|---|---|
| `pdfgui-ni` | PDFgui Ni tutorial | Demonstrate PDF calculation and small-box refinement | Reproduce the calculated PDF on the published grid, including finite-Q and resolution settings where available | Forward convergence validated; 0.55% RMS/range representation difference |
| `rmcprofile-caf2` | RMCProfile7, section 3.4 | Compare histogram and back-Fourier-transform PDFs at realistic and reduced `Qmax` | Reproduce the X-ray and neutron CaF2 curves and their `Qmax` dependence from the same configuration | Published Qmax protocol reproduced; native output comparison awaits licensed inputs |
| `rmcprofile-srtio3` | Official RMCProfile Exercise 4 | Demonstrate neutron F(Q), G(r), Bragg, and big-box modelling of the 105 K transition | Pointwise zero-move comparison of neutron F(Q), total G(r), and partial RDFs from the identical tutorial configuration | Native forward parity validated: 1.03% F(Q) RMS/range |
| `epsr-ga` | EPSR26 workshop LiquidGa50C example | Demonstrate empirical-potential refinement for a simple liquid | First match the zero-move partial RDF and neutron interference function, then reproduce the fitted ensemble | Forward parity validated; from ten independent paired starts rsmith reaches 3.0% residual 11.11x faster and reaches 2.0% in every run while native EPSR does not reach it within the budget |
| `epsr-silica` | Published EPSR silica-glass example | Demonstrate a chemically plausible network-glass model | First reproduce native EPSR totals and structural distributions, then compare MLIP-HRMC at matched fit quality | Forward/reference-potential gates and symmetric cross-program smoke passed for rsmith, EPSR26, RMCProfile, PDFgui, Pedone, and GAP; calibrated multi-seed production pending; EPSR publication permission required |
| `rmcprofile-ceo2-mlip` | Cuillier, Tucker & Zhang (2024) | Recover oxygen-vacancy ordering with a LAMMPS/MLIP constraint | Reproduce the vacancy-pair convergence experiment, then compare full, local, and incremental energy deltas and wall time | Requires vacancy/swap moves and upstream inputs |
| `hp-geo2` | High-pressure neutron/isotope-substitution literature | Resolve the pressure-driven GeO4/GeO5/GeO6 network transformation | Compare pure RMC, EPSR/Dissolve, pair-HRMC, and MLIP-HRMC with held-out contrast validation | Design complete; data and pressure-valid MLIP pending |

## Primary references

- C. L. Farrow *et al.*, "PDFfit2 and PDFgui: computer programs for studying
  nanostructure in crystals", *J. Phys.: Condens. Matter* **19**, 335219
  (2007), <https://doi.org/10.1088/0953-8984/19/33/335219>.
- W. A. Slawinski *et al.*, "RMCProfile7: reverse Monte Carlo for multiphase
  systems", *J. Appl. Cryst.* **57**, 1251-1262 (2024),
  <https://doi.org/10.1107/S1600576724004175>.
- P. Cuillier, M. G. Tucker and Y. Zhang, "Integrating machine learning
  interatomic potentials with hybrid reverse Monte Carlo structure refinements
  in RMCProfile", *J. Appl. Cryst.* **57**, 1780-1788 (2024),
  <https://doi.org/10.1107/S1600576724009282>.
- D. T. Bowron, "Experimentally consistent atomistic modeling of bulk and
  local structure in liquids and disordered materials by empirical potential
  structure refinement", *Pure Appl. Chem.* **80**, 1211-1227 (2008),
  <https://doi.org/10.1351/pac200880061211>.
- T. G. A. Youngs, "Dissolve: next generation software for the interrogation
  of total scattering data by empirical potential generation", *Molecular
  Physics* **117**, 3464-3477 (2019),
  <https://doi.org/10.1080/00268976.2019.1651918>.

## Acceptance rules

### Deterministic forward calculations

1. Use the same atomic configuration, composition, density, grids, scattering
   convention, form-factor or scattering-length table, termination function,
   and resolution correction.
2. Archive both the native reference output and the converted comparison file.
3. Report maximum absolute error, RMS error, and the error normalized by the
   dynamic range.  A tolerance is fixed in the case manifest before inspecting
   the result.
4. Compare against an independent direct-summation oracle in addition to the
   external program, so agreement cannot be caused by a shared approximation.

### Stochastic refinements

1. Use identical starting configurations and fitted data wherever the move
   sets permit this.
2. Run at least ten independent seeds for scientific comparisons.
3. Report the full distribution of fit residual, energy, coordination,
   angular, topological, and held-out observables; do not select only the best
   run.
4. Compare methods at matched data-fit quality, preferably as a Pareto curve of
   data residual versus independently evaluated structural error.
5. Include initialization and model-loading cost in end-to-end timing.

## Per-case artifact contract

Each benchmark directory must eventually contain:

```text
README.md              provenance, licence, and reproduction instructions
manifest.toml          versions, grids, tolerances, seeds, and hardware fields
reference/             unmodified upstream inputs and outputs when redistributable
rsmith/                rsmith configuration and converted inputs
scripts/               deterministic conversion and comparison scripts
expected/              frozen machine-readable metrics
```

If upstream files cannot legally be redistributed, `reference/README.md` must
give a stable DOI or download location and checksums for the files used.

## Explicit initial scope

The first paper does not need to reproduce RMCProfile's multiphase Bragg-profile
or rigid-molecule demonstrations.  Those features solve different problems and
are not currently implemented in `rsmith`.  The initial comparison is limited
to the overlapping domain: single-phase atomic liquids, glasses, total
scattering/PDF forward models, and potential-regularized big-box refinement.
