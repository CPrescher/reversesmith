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
7.0% acceptance; all three empirical pair potentials changed (maximum changes
0.0671, 0.2012, and 0.2025 eV for Si-Si, Si-O, and O-O). The post-MC X-ray and
neutron RMS residuals were 0.2391 and 0.2083. These values prove that the full
two-contrast reference-plus-empirical path executes; a single epoch is not a
converged EPSR reproduction or a speed comparison.

The deterministic forward and reference-potential blockers are therefore
closed. The remaining reproduction gate is a matched, multi-epoch native and
rsmith ensemble with frozen seed-to-seed equivalence margins.

The scientific claim is not "smaller chi-squared than EPSR". The meaningful
test is whether MLIP-HRMC gives better held-out or chemically diagnostic
structure at the same data agreement.

## Files

- `protocol.toml`: frozen methods, stages, contrasts, seeds, and observables;
- `PERMISSION_REQUEST.md`: ready-to-send publication/reuse request;
- `fetch_reference.sh`: explicit, locally gated acquisition with SHA-256 log;
- `import_local_reference.sh`: local-only DTBsilicaNX importer and hasher;
- `scripts/run_native_zero_move.py`: clean one-configuration EPSR26 runner;
- `scripts/compare_forward.py`: exact-structure rsmith comparison and
  independent RDF oracle;
- `scripts/verify_forward.py`: committed forward-regression verifier;
- `scripts/run_reference_potential_smoke.py`: independent EPSR reference-
  potential reconstruction, native curve/energy gate, normalization
  conversion, and one-epoch joint neutron/X-ray rsmith smoke run;
- `expected/native-forward.toml`: provenance hashes, observations, and
  regression guards;
- `expected/reference-potential-smoke.toml`: frozen potential and smoke-run
  observations plus regression guards;
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
```
