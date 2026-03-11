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

## Choosing `rdf_cutoff`

The `rdf_cutoff` is a **critical** parameter that directly affects the quality of the computed S(Q). It determines how much of the real-space pair distribution function g(r) is used in the Fourier transform to S(Q).

### Impact on S(Q)

S(Q) is computed as:

```
S(Q) = 1 + (4*pi*rho0 / Q) * integral_0^rmax r * [g(r) - 1] * sin(Qr) dr
```

Truncating the integral at a finite `rmax = rdf_cutoff` is equivalent to multiplying g(r) by a step function, which introduces artifacts in S(Q):

- **Truncation ripples**: Oscillatory artifacts, especially at low Q. The Lorch window mitigates these but cannot fully compensate.
- **First Sharp Diffraction Peak (FSDP)**: The FSDP (typically Q ~ 1.5--3 A^-1 in glasses) encodes medium-range order extending to 15--25 A in real space. A short cutoff systematically reduces the FSDP intensity and broadens it.
- **Low-Q region**: Features below ~5 A^-1 are most affected by the cutoff.

### Practical recommendations

| `rdf_cutoff` | Effect | Use case |
|---|---|---|
| 10--12 A | Fast, but first S(Q) peak may be inaccurate | Quick tests, high-Q-only fitting |
| 15--20 A | Good balance of accuracy and speed | Most production runs |
| 20--25 A | Excellent S(Q) quality including FSDP | When low-Q accuracy is critical |
| > L/2 | Invalid -- exceeds minimum image convention | Never |

where L is the shortest box dimension. Use `--compute-sq-only` to check S(Q) quality before committing to a refinement run.

### Comparison with LAMMPS

If you compare the initial S(Q) from rsmith with one computed from the same structure in LAMMPS (e.g., via `compute rdf`), differences typically arise from:

1. **RDF cutoff**: LAMMPS often uses 20--25 A cutoffs for `compute rdf` while rsmith defaults to 10 A. This is usually the dominant source of discrepancy, particularly for the first peak.
2. **Trajectory averaging vs. single snapshot**: LAMMPS `fix ave/time` averages g(r) over many frames, producing smoother curves. Reversesmith uses a single snapshot, so the g(r) has more statistical noise.
3. **Form factor source**: rsmith uses Cromer-Mann parameterization; other tools may use xraylib. These are very similar but not identical at high Q.

To achieve agreement with LAMMPS-derived S(Q), increase `rdf_cutoff` to match the LAMMPS `compute rdf` cutoff.

### Performance trade-off

Increasing `rdf_cutoff` costs performance in two ways:

1. **More histogram bins**: Each atom move updates more bins in the RDF histogram, increasing the per-move cost of the S(Q) update.
2. **More neighbors**: The cell list covers a larger volume, so each atom has more neighbors to iterate over.

Roughly, doubling the cutoff increases neighbor count by ~8x and histogram bins by ~2x. In practice, going from 11 A to 20 A makes each RMC move ~2--3x slower. For a 10,000-atom system this is still fast (< 1 ms per move). See [Performance Notes](../performance.md) for details.

### Example

```toml
[sq]
qmax = 20.0
lorch = true
rdf_cutoff = 20.0    # Good for accurate FSDP
rdf_nbins = 1000     # dr = 0.02 A
```
