# `[rmc]` -- Refinement Parameters

```toml
[rmc]
max_moves = 1_000_000       # Total MC moves to attempt
max_step = 0.1              # Initial maximum atomic displacement (A)
checkpoint_every = 50_000   # Save checkpoint every N moves
seed = 42                   # RNG seed (optional; random if omitted)
print_every = 1000          # Print status every N moves
target_acceptance = 0.3     # Target acceptance rate for step adaptation
adjust_step_every = 5000    # Adjust step size every N moves

# Simulated annealing (omit for standard RMC at T=1)
anneal_start = 2.0          # Starting temperature
anneal_end = 0.01           # Final temperature
anneal_steps = 200_000      # Moves over which to anneal (rest at anneal_end)

# Early stopping
convergence_threshold = 1e-4  # Stop if chi2 improves less than this...
convergence_window = 50_000   # ...over this many moves (0 = disabled)
```

## Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `max_moves` | Integer | 1,000,000 | Total MC moves to attempt (not additional moves on `--resume`) |
| `max_step` | Float | 0.1 | Initial maximum displacement per dimension (A) |
| `checkpoint_every` | Integer | 50,000 | Checkpoint interval (moves) |
| `seed` | Integer | random | RNG seed for reproducibility. When omitted, a random seed is generated from system entropy and logged. Can also be overridden via `--seed N` on the command line. |
| `print_every` | Integer | 1,000 | Status output interval (moves) |
| `target_acceptance` | Float | 0.3 | Target acceptance ratio for step adaptation |
| `adjust_step_every` | Integer | 5,000 | Step size adjustment interval (moves) |
| `anneal_start` | Float | 1.0 | Starting temperature |
| `anneal_end` | Float | 1.0 | Final temperature |
| `anneal_steps` | Integer | max_moves | Moves over which to anneal |
| `convergence_threshold` | Float | 0.0 | Early stopping threshold (0 = disabled) |
| `convergence_window` | Integer | 50,000 | Window for convergence check (moves) |

## Simulated annealing

Standard RMC uses T=1, which samples configurations consistent with the data. For optimization, use simulated annealing to drive chi2 lower:

| Phase | Temperature | Behavior |
|-------|-------------|----------|
| Early (T > 1) | High | Explores broadly, escapes local minima |
| Middle (T ~ 0.1--1) | Moderate | Selective exploration |
| Late (T < 0.1) | Low | Greedy optimization, mostly downhill moves |

The temperature follows an exponential schedule:

```
T(n) = T_start * (T_end / T_start) ^ (n / anneal_steps)
```

For n > anneal_steps, T stays at anneal_end. The remaining moves optimize at constant low temperature.

**Guidelines:**
- `anneal_start = 1.0--2.0` -- higher values explore more aggressively
- `anneal_end = 0.01--0.1` -- lower values optimize more strictly
- `anneal_steps` -- typically 10--30% of `max_moves`

## Adaptive step size

- **During annealing:** step scales with temperature, capped at `max_step`
- **After annealing:** step adapts to maintain the target acceptance rate. If acceptance > target + 5%, step grows by 5%; if below target - 5%, step shrinks by 5%
- Step is always clamped to [0.001, 2.0] A

## Convergence detection

When `convergence_threshold > 0`, the run monitors chi2 improvement over sliding windows of `convergence_window` moves. If the improvement is less than the threshold, the run stops early. During annealing, convergence checking is deferred until T reaches `anneal_end`.

## Best-structure tracking

The lowest-chi2 configuration encountered during the entire run is tracked and restored at completion. This prevents the output from being worse than intermediate states (which can happen during annealing or random walks at high temperature).
