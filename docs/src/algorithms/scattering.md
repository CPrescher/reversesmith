# Scattering Functions

## S(Q) computation (forward transform)

Partial structure factors S_ab(Q) are computed from partial pair distribution functions g_ab(r):

```
S_ab(Q) = 1 + (4*pi*rho0 / Q) * integral[ r * (g_ab(r) - 1) * W(r) * sin(Q*r) dr ]
```

where W(r) is the optional Lorch modification function:

```
W(r) = sin(pi*r/r_max) / (pi*r/r_max)
```

This suppresses truncation ripples in S(Q) caused by the finite RDF cutoff.

For the initial computation, a DST-II (Discrete Sine Transform Type II) via FFT is used for efficiency, with 8x zero-padding for accurate interpolation onto the requested Q grid.

During RMC, the same formula is evaluated by direct summation using precomputed sin(Q*r) tables, enabling incremental updates.

## X-ray weighting (Faber-Ziman convention)

The total X-ray structure factor combines partial S_ab(Q) with Q-dependent weights:

```
S_X(Q) = sum_{a<=b} (2 - delta_ab) * c_a * c_b * f_a(Q) * f_b(Q) / <f(Q)>^2 * S_ab(Q)
```

where:
- c_a = concentration (atomic fraction) of species a
- f_a(Q) = X-ray atomic form factor (Cromer-Mann parameterization)
- \<f(Q)\> = concentration-weighted average form factor
- delta_ab = Kronecker delta (1 if a=b, 0 otherwise)

## Neutron weighting

Neutron weighting uses Q-independent coherent scattering lengths b_a instead of form factors:

```
S_N(Q) = sum_{a<=b} (2 - delta_ab) * c_a * c_b * b_a * b_b / <b>^2 * S_ab(Q)
```

Since the scattering lengths are Q-independent, the neutron weights are constant vectors (unlike X-ray weights which vary with Q through the form factors).

## Simultaneous X-ray and neutron fitting

X-ray and neutron S(Q) can be fitted simultaneously. Each dataset gets its own weight vector and total S(Q), computed from the same partial structure factors but with different weighting:

- X-ray: Q-dependent form factor weights (see above)
- Neutron: Q-independent scattering length weights

Both contribute independently to the total chi2, with their own `weight`, `sigma`, and `fit_min`/`fit_max`. This is particularly powerful because X-ray and neutron contrast differently for different atom pairs -- for example, in oxide glasses, oxygen is nearly invisible to X-rays but has a large neutron cross-section.

## g(r) from S(Q) (inverse transform)

The model g(r) for comparison with experimental data is computed by inverse Fourier transform of the corresponding total S(Q) (X-ray or neutron weighted):

```
g_X(r) = 1 + (dQ / (2*pi^2 * rho0 * r)) * sum_k Q_k * W(Q_k) * (S_X(Q_k) - 1) * sin(Q_k * r)
```

The sum runs over Q points up to Q_max (matching the experimental data range). An optional Q-space Lorch window damps high-Q contributions:

```
W(Q) = sin(pi*Q/Q_max) / (pi*Q/Q_max)
```

This is precomputed as a matrix-vector product for efficient incremental updates:

```
g(r_i) = 1 + sum_k M[i][k] * (S(Q_k) - 1)
```

When S(Q) changes by delta_S(Q), the g(r) update is simply:

```
delta_g(r_i) = sum_k M[i][k] * delta_S(Q_k)
```

## Effect of RDF cutoff on S(Q)

The finite `rdf_cutoff` (r_max) truncates the Fourier integral:

```
S(Q) = 1 + (4*pi*rho0/Q) * integral_0^r_max r * [g(r) - 1] * W(r) * sin(Qr) dr
```

In a glass, g(r) oscillates around 1 with decreasing amplitude, typically reaching 1 +/- 0.01 by ~15 A and 1 +/- 0.005 by ~25 A. Cutting off these oscillations removes information:

- **Q > 5 A^-1**: Dominated by short-range correlations (first and second coordination shells at 1.5--4 A). Insensitive to the cutoff as long as `rdf_cutoff` > ~8 A.
- **Q ~ 1.5--5 A^-1**: The First Sharp Diffraction Peak (FSDP) encodes intermediate-range order (ring statistics, network connectivity). These real-space correlations extend to 15--25 A. Insufficient cutoff reduces the FSDP intensity and broadens it.
- **Q < 1 A^-1**: Very sensitive to long-range density fluctuations. Difficult to converge; the Lorch window is essential here.

The Lorch window `W(r) = sinc(pi*r/r_max)` smoothly damps g(r) to zero at r_max, converting the sharp step into a gradual taper. This eliminates truncation ripples but attenuates the signal near r_max, effectively reducing the useful range by ~20%. A system with `rdf_cutoff = 20 A` and Lorch damping captures roughly the same information as `rdf_cutoff = 16 A` without damping, but with much cleaner S(Q).

**Recommendation**: For production refinement of glasses, use `rdf_cutoff >= 15 A` (with Lorch). For careful comparison with experiment, especially the FSDP, use 20--25 A.

## DST-II implementation

The forward Fourier transform is computed using a Discrete Sine Transform Type II (DST-II) via FFT, with 8x zero-padding for accurate interpolation:

1. Compute `f(r_i) = r_i * [g(r_i) - 1] * W(r_i)` on the RDF grid
2. Zero-pad to `n_fft = next_pow2(8 * n_bins)` points
3. Compute forward FFT, extract DST-II coefficients via twiddle factors
4. Evaluate `S(Q_k) = 1 + 4*pi*rho0 * DST[k] / Q_k`
5. Linearly interpolate from the fine FFT Q-grid onto the requested Q points

This gives S(Q) on a fine grid (dQ ~ 0.02 A^-1 for typical parameters) from which the user-requested Q points are interpolated.

During RMC, the same formula is evaluated by **direct summation** using precomputed `sin(Q*r)` tables, which enables efficient incremental updates when a single atom is moved.

## g(r) consistency

The model and experimental g(r) must be derived using identical transform parameters:

- **Q_max**: the maximum Q used in the FT
- **Lorch window**: whether and how it is applied
- **Number density rho0**: computed from the model structure

Mismatched parameters cause systematic differences in peak heights, biasing the refinement. Set `qmax` and `lorch` on `[data.xray_gr]` (or `[data.neutron_gr]` when available) to match the values used to produce the experimental g(r).
