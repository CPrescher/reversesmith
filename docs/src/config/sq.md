# `[sq]` -- S(Q) Computation Parameters

Controls the model S(Q) grid and the underlying RDF computation.

```toml
[sq]
qmin = 0.0                  # Minimum Q for model S(Q) grid (1/A)
qmax = 20.0                 # Maximum Q (1/A)
nq = 500                    # Number of Q grid points
lorch = true                # Apply Lorch window in r-space for S(Q) computation
rdf_cutoff = 11.0           # Maximum r for RDF histogram (A)
rdf_nbins = 550             # Number of RDF histogram bins
```

## Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `qmin` | Float | 0.3 | Minimum Q for S(Q) grid (1/A) |
| `qmax` | Float | 20.0 | Maximum Q (1/A) |
| `nq` | Integer | 500 | Number of uniformly spaced Q points |
| `lorch` | Bool | true | Lorch modification function in r-space |
| `rdf_cutoff` | Float | 10.0 | Maximum r for RDF histograms (A) |
| `rdf_nbins` | Integer | 500 | Number of RDF bins |

## Guidelines

- The Q grid should extend beyond the experimental Q range to avoid edge effects during interpolation.
- `rdf_cutoff` determines the maximum real-space range used in the Fourier transform. Larger values improve low-Q accuracy but increase computation cost. Should be less than half the box length.
- `rdf_nbins` controls the resolution of the RDF histogram. The bin width is `rdf_cutoff / rdf_nbins`. A resolution of ~0.02 A (e.g., 550 bins for 11 A cutoff) is typical.
- The Lorch window suppresses truncation ripples in S(Q) caused by the finite `rdf_cutoff`. It is generally recommended.
- When using [pair potentials](potentials.md), the potential cutoff must not exceed `rdf_cutoff` (the same cell list is used for both).
