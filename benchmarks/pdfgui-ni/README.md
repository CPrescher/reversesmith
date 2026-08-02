# PDFgui / DiffPy Ni PDF benchmark

This is the first cross-code benchmark because Ni is the introductory PDFgui
example and also appears in the official DiffPy-CMI examples.  It is a
real-space forward-model benchmark, not evidence that a big-box RMC method is
better than a small-box crystallographic refinement.

## Provenance

- PDFgui tutorial: <https://www.diffpy.org/diffpy.pdfgui/tutorial.html>
- DiffPy-CMI example repository: <https://github.com/diffpy/cmi_exchange/tree/main/cmi_scripts/fitNiPDF>
- Pinned example commit: `b76dda7dd899b0a391009225d11cf5651295e887`
- `ni.cif` SHA-256: `ee4d128747da0e0ab8cb0deca8170d9f39d612f5756a08007348d1d96f31b1b3`
- `ni-q27r60-xray.gr` SHA-256: `87c484a36e2a6f987381b586fb7e124139830a897c2b7d4c7c7857c7128e93f2`

The CIF states that its Crystallography Open Database content is public domain.
The example repository did not expose a repository-level licence in the pinned
checkout, so the experimental PDF is not vendored here.  Review its reuse terms
before distributing it with a publication archive.

## Reproduction

```bash
./fetch_reference.sh
python3 scripts/prepare_reference.py
python3 scripts/calculate_pdfgui_reference.py
python3 scripts/generate_ni_supercell.py --cells 12 --lattice 3.52 --uiso 0.00126651 --seed 20260802
cargo run --release --manifest-path ../../Cargo.toml -- rsmith/forward.toml --output-dir results
python3 ../scripts/compare_curves.py \
  reference/pdfgui-calculated.dat results/start_total_fr.dat \
  --fit-min 1.7 --fit-max 20 --output results/forward-metrics.json

# Reduce finite sampling noise using independent static thermal replicas.
python3 scripts/run_ensemble.py --replicas 10
python3 ../scripts/compare_curves.py \
  reference/pdfgui-calculated.dat results/ensemble_mean_fr.dat \
  --fit-min 1.7 --fit-max 20 --output results/ensemble-forward-metrics.json

# Reproduce the box-size, RDF-bin, and Q-grid convergence table.
python3 scripts/run_convergence.py
```

The PDFgui tutorial specifies an initial `Qdamp` of 0.08 1/A and a lower fit
bound of 1.7 A. The bundled `Ni.stru` at PDFgui commit
`9496a699c130efaa9ce620c29e7fe95b88fcad6f` supplies `a = 3.52 A` and
`Uiso = 0.00126651 A^2`; those values define the unrefined forward comparison.

The primary forward-model comparison is against `pdfgui-calculated.dat`, not
against the experiment. This isolates numerical agreement at identical model
parameters. The experimental curve is reserved for the subsequent refinement
comparison.

## Convergence result

The seven-case, three-replica convergence study varies the box from 6,912 to
16,384 atoms, RDF spacing from 0.020 to 0.005 A, and nominal Q spacing from
0.020 to 0.005 1/A. Across those cases, RMS error relative to the PDFgui curve
is 0.537--0.580% of its dynamic range and the maximum absolute difference is
2.120--2.131%. The stability of the residual under all three refinements shows
that it is not dominated by finite-box sampling or either numerical grid. It
is instead the expected representation difference between PDFfit2's analytic
crystal peak model and rsmith's explicitly sampled histogram/back-transform
model.

The frozen regression record is `expected/convergence.toml`. Its 0.60% RMS and
2.2% maximum-difference bounds were selected after this study and are therefore
labelled as post-result engineering regression limits, not preregistered
scientific acceptance criteria.
