# `[data]` -- Experimental Data

Each dataset is a two-column whitespace-separated file. Lines starting with `#` are skipped. At least one dataset must be specified.

## `[data.xray_sq]` -- X-ray Structure Factor

```toml
[data.xray_sq]
file = "experimental.sq"    # Two columns: Q (1/A), S(Q)
weight = 1.0                # Relative weight in total chi2 (default: 1.0)
sigma = 0.02                # Uncertainty in S(Q) units (omit to auto-estimate from data)
# sigma_column = 3          # Or read a per-point uncertainty from this 1-based column
sigma_alpha = 0.05          # Linear Q-scaling: sigma *= 1 + alpha*Q (default: 0)
fit_min = 0.5               # Only fit Q > fit_min (default: 0, fit all)
fit_max = 18.0              # Only fit Q < fit_max (default: inf, fit all)
convention = "sq"           # "sq", "iq", or "fq" (default: "sq")
```

## `[data.neutron_sq]` -- Neutron Structure Factor

Same format as `xray_sq`. Uses coherent neutron scattering lengths instead of X-ray form factors.

```toml
[data.neutron_sq]
file = "neutron.sq"
weight = 1.0
sigma = 0.02
sigma_alpha = 0.05
scattering_lengths = { Ge = 10.0, O = 5.803 } # fm; optional isotope contrast
```

Both X-ray and neutron S(Q) can be fitted simultaneously.

## `[data.xray_gr]` -- X-ray Pair Distribution Function

```toml
[data.xray_gr]
file = "experimental.gr"    # Two columns: r (A), g(r)
weight = 0.3                # Relative weight in total chi2 (default: 1.0)
sigma = 0.02                # Uncertainty in g(r) units (omit to auto-estimate from data)
fit_min = 0.0               # Only fit r > fit_min (default: 0)
fit_max = 7.0               # Only fit r < fit_max (default: inf)
qmax = 17.97                # Q_max for inverse FT (default: experimental S(Q) Qmax)
lorch = true                # Apply Lorch window in Q-space (default: true)
qdamp = 0.04                # Gaussian PDF resolution damping in 1/A (default: 0)
```

## `[data.xray_fr]` -- X-ray Reduced Pair Distribution Function

A density-independent alternative to `[data.xray_gr]`. The reduced pair distribution function f(r) is defined as:

```
f(r) = (2/π) ∫ Q·[S(Q)-1]·W(Q)·sin(Qr) dQ
```

Unlike g(r), f(r) does not depend on the number density ρ₀, making it useful when comparing ensemble runs at different densities (`--density`). The relationship is `f(r) = 4πρ₀r·[g(r) - 1]`.

```toml
[data.xray_fr]
file = "experimental.fr"    # Two columns: r (A), f(r)
weight = 0.3                # Relative weight in total chi2 (default: 1.0)
sigma = 0.02                # Uncertainty in f(r) units (omit to auto-estimate from data)
fit_min = 0.0               # Only fit r > fit_min (default: 0)
fit_max = 7.0               # Only fit r < fit_max (default: inf)
qmax = 17.97                # Q_max for inverse FT (default: experimental S(Q) Qmax)
lorch = true                # Apply Lorch window in Q-space (default: true)
```

**Note**: `[data.xray_gr]` and `[data.xray_fr]` are mutually exclusive -- you cannot use both in the same config.

## Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `file` | String | Required | Path to data file (relative to config) |
| `weight` | Float | 1.0 | Relative weight in total chi2 |
| `sigma` | Float | auto | Per-point uncertainty in absolute units (same units as S(Q) or g(r)). If omitted, estimated from data noise using windowed second finite differences |
| `sigma_column` | Integer | none | One-based column containing a positive per-point uncertainty (3 or greater). Mutually exclusive with `sigma` |
| `sigma_alpha` | Float | 0.0 | Linear Q-scaling factor alpha: `sigma(Q) *= 1 + alpha * Q`. Only applies to S(Q) datasets |
| `fit_min` | Float | 0.0 | Lower bound of fitting range |
| `fit_max` | Float | inf | Upper bound of fitting range |
| `qmax` | Float | auto | Q_max for g(r) inverse FT (g(r) only) |
| `lorch` | Bool | true | Lorch window in Q-space (g(r) only) |
| `qdamp` | Float | 0.0 | PDF resolution envelope in 1/A: the real-space signal is multiplied by `exp[-0.5*(qdamp*r)^2]`. Valid for `xray_gr` and `xray_fr` |
| `convention` | String | `"sq"` | Reciprocal-space input convention: `"sq"` for S(Q), `"iq"` for S(Q)-1, or `"fq"` for Q[S(Q)-1] |
| `scattering_lengths` | Table | natural elements | Coherent neutron scattering lengths in fm, keyed by structure species. Valid only for `neutron_sq`; unspecified species use built-in natural-element values |

`convention` applies to X-ray and neutron reciprocal-space datasets. It lets
`rsmith` transform its calculated Faber-Ziman S(Q) before evaluating the residual,
so the experimental file does not need to be pre-converted. It is ignored for
real-space datasets.

For isotope-substitution measurements, define `scattering_lengths` separately
for the neutron dataset. Values are effective coherent scattering lengths for
that sample, in fm. Every key must match a species in the loaded structure;
species omitted from the table use the natural-element value.

## Sigma estimation

Sigma is in **absolute units** -- the same units as the data it applies to (dimensionless for S(Q), dimensionless for g(r)). A value of `sigma = 0.02` means you expect the data to be accurate to ±0.02 in S(Q). It enters the chi2 cost function as:

```
chi2 += (S_calc(Q) - S_exp(Q))^2 / sigma(Q)^2
```

When `sigma_column` is set, rsmith reads one positive uncertainty for each data
point from that column. This is preferred for quantitative comparisons when the
experimental reduction provides uncertainties. If neither `sigma` nor
`sigma_column` is given, rsmith estimates per-point uncertainties directly from
the data using windowed second finite differences. The second difference `d2[i]
= y[i+1] - 2y[i] + y[i-1]` removes the smooth signal and isolates noise. The
local RMS of d2 in a sliding window, scaled by 1/sqrt(6), gives the noise estimate
at each point.

This automatically produces larger sigma where the data is noisy (typically at high Q) and smaller sigma where it is smooth.

## Q-dependent sigma scaling

S(Q) data is typically noisier at high Q. The `sigma_alpha` parameter provides additional Q-dependent relaxation on top of the base sigma (whether constant or auto-estimated):

```
sigma_effective(Q) = sigma_base(Q) * (1 + alpha * Q)
```

For example, with `sigma = 0.02`, `sigma_alpha = 0.05`, and Q_max = 20 A^-1:
- At Q = 0: sigma = 0.02 (unchanged)
- At Q = 10: sigma = 0.02 * 1.5 = 0.03
- At Q = 20: sigma = 0.02 * 2.0 = 0.04

This prevents the fit from chasing noise in the high-Q tail at the expense of the physically important low-Q features.

This option only applies to S(Q) datasets (`xray_sq`, `neutron_sq`) and is ignored for g(r).

## g(r) consistency

The model g(r) is computed by inverse Fourier transform of the corresponding total S(Q) (X-ray weighted for `xray_gr`). The `qmax` and `lorch` parameters **must match** those used to derive the experimental g(r) from experimental S(Q). Otherwise peak heights will be inconsistent -- the model g(r) may have sharper or broader peaks than the experimental g(r), creating a systematic bias.

## Balancing S(Q) and g(r)

Since g(r) is derived from S(Q), they contain correlated information. When fitting both simultaneously, the g(r) `weight` should typically be smaller (0.01--0.3) to avoid double-counting. Monitor both components in the status output to ensure neither dominates.
