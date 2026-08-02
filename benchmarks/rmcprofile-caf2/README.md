# RMCProfile CaF2 forward-PDF benchmark

RMCProfile7 Figure 11 compares its histogram PDF with the back-Fourier-transform
PDF for neutron CaF2 at Qmax = 30 1/A and at reduced Qmax values from 12 to
20 1/A. This directly probes the systematic difference now visible in the Ni
benchmark and is therefore the next forward-model case.

Primary reference: W. A. Slawinski et al., *J. Appl. Cryst.* **57** (2024),
1251–1262, DOI `10.1107/S1600576724004175`.

## What is reproduced now

The executable benchmark reconstructs the Figure 11 protocol without copying
or digitizing its curves. It builds a 15 x 15 x 15 conventional-cell fluorite
configuration (40,500 atoms, 81.93 A per side), computes the neutron-weighted
histogram PDF, and back-transforms the same calculated S(Q) at Qmax = 12, 14,
16, 18, 20, and 30 1/A. Natural-element coherent scattering lengths are 4.70
fm for Ca and 5.654 fm for F.

The paper does not state the numerical lattice parameter or ADPs used for its
plotted configuration. This reconstruction therefore records its assumptions:
`a = 5.462 A`, `Uiso(Ca) = 0.005 A^2`, and `Uiso(F) = 0.007 A^2`. They are not
presented as recovered RMCProfile inputs.

Run from the repository root (Python 3.11+ with NumPy is required):

```bash
cargo build --release
python3 benchmarks/rmcprofile-caf2/scripts/run_qmax_series.py
```

The deterministic result over 1.5--15 A is:

| Qmax (1/A) | RMS / PDF range | maximum / PDF range |
|---:|---:|---:|
| 12 | 5.890% | 16.050% |
| 14 | 4.776% | 13.414% |
| 16 | 3.726% | 12.531% |
| 18 | 2.473% | 7.471% |
| 20 | 1.372% | 4.259% |
| 30 | 0.275% | 1.375% |

This reproduces the published conclusion: reduced Qmax produces progressively
larger termination ripples and peak broadening, whereas the Qmax = 30 1/A
back-transform is nearly indistinguishable from the histogram representation.
The script enforces the monotonic trend and the frozen regression envelopes in
`expected/qmax-series.toml`.

## What remains for a native cross-code comparison

Required upstream artifacts are:

- the exact CaF2 configuration used for Figure 11;
- RMCProfile7 input files and executable version;
- histogram and back-transform output curves at every reported Qmax;
- resolution parameters and neutron scattering-length convention.

These files were not linked from the article text. The public RMCProfile disk
image also presents a licence/usage notice before mounting, so it has not been
accepted automatically. Request the files from the authors/RMCProfile team or
accept the distribution terms explicitly, rather than treating digitized plot
pixels as native reference output. Once acquired, run the same configuration
through rsmith and compare all curves using `../scripts/compare_curves.py`. No
native cross-code acceptance tolerance is set until those curves are available.
