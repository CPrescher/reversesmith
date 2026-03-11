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

For the initial computation and `--compute-sq-only` mode, a DST-II (Discrete Sine Transform Type II) via FFT is used for efficiency, with 8x zero-padding for accurate interpolation onto the requested Q grid.

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

## g(r) from S(Q) (inverse transform)

The model g(r) for comparison with experimental data is computed by inverse Fourier transform of the total X-ray S(Q):

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

## g(r) consistency

The model and experimental g(r) must be derived using identical transform parameters:

- **Q_max**: the maximum Q used in the FT
- **Lorch window**: whether and how it is applied
- **Number density rho0**: computed from the model structure

Mismatched parameters cause systematic differences in peak heights, biasing the refinement. Set `qmax` and `lorch` on `[data.xray_gr]` to match the values used to produce the experimental g(r).
