# `[data]` -- Experimental Data

Each dataset is a two-column whitespace-separated file. Lines starting with `#` are skipped. At least one dataset must be specified.

## `[data.xray_sq]` -- X-ray Structure Factor

```toml
[data.xray_sq]
file = "experimental.sq"    # Two columns: Q (1/A), S(Q)
weight = 1.0                # Relative weight in total chi2 (default: 1.0)
sigma = 0.02                # Experimental uncertainty (default: 0.01)
fit_min = 0.5               # Only fit Q > fit_min (default: 0, fit all)
fit_max = 18.0              # Only fit Q < fit_max (default: inf, fit all)
```

## `[data.neutron_sq]` -- Neutron Structure Factor

Same format as `xray_sq`. Uses coherent neutron scattering lengths instead of X-ray form factors.

```toml
[data.neutron_sq]
file = "neutron.sq"
weight = 1.0
sigma = 0.02
```

Both X-ray and neutron S(Q) can be fitted simultaneously.

## `[data.xray_gr]` -- X-ray Pair Distribution Function

```toml
[data.xray_gr]
file = "experimental.gr"    # Two columns: r (A), g(r)
weight = 0.3                # Relative weight in total chi2 (default: 1.0)
sigma = 0.02                # Experimental uncertainty (default: 0.01)
fit_min = 0.0               # Only fit r > fit_min (default: 0)
fit_max = 7.0               # Only fit r < fit_max (default: inf)
qmax = 17.97                # Q_max for inverse FT (default: experimental S(Q) Qmax)
lorch = true                # Apply Lorch window in Q-space (default: true)
```

## Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `file` | String | Required | Path to data file (relative to config) |
| `weight` | Float | 1.0 | Relative weight in total chi2 |
| `sigma` | Float | 0.01 | Per-point experimental uncertainty |
| `fit_min` | Float | 0.0 | Lower bound of fitting range |
| `fit_max` | Float | inf | Upper bound of fitting range |
| `qmax` | Float | auto | Q_max for g(r) inverse FT (g(r) only) |
| `lorch` | Bool | true | Lorch window in Q-space (g(r) only) |

## g(r) consistency

The model g(r) is computed by inverse Fourier transform of the model's total X-ray S(Q). The `qmax` and `lorch` parameters **must match** those used to derive the experimental g(r) from experimental S(Q). Otherwise peak heights will be inconsistent -- the model g(r) may have sharper or broader peaks than the experimental g(r), creating a systematic bias.

## Balancing S(Q) and g(r)

Since g(r) is derived from S(Q), they contain correlated information. When fitting both simultaneously, the g(r) `weight` should typically be smaller (0.01--0.3) to avoid double-counting. Monitor both components in the status output to ensure neither dominates.
