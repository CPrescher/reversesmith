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

# Simulated annealing (omit for constant T = 0.1)
anneal_start = 2.0          # Starting temperature
anneal_end = 0.1            # Final temperature (also used as constant T when annealing is omitted)
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
| `anneal_start` | Float | 0.1 | Starting temperature |
| `anneal_end` | Float | 0.1 | Final temperature (constant T when annealing is omitted) |
| `anneal_steps` | Integer | max_moves | Moves over which to anneal |
| `convergence_threshold` | Float | 0.0 | Early stopping threshold (0 = disabled) |
| `convergence_window` | Integer | 50,000 | Window for convergence check (moves) |

## Temperature

The temperature T controls the Metropolis acceptance criterion: moves that increase the cost by delta are accepted with probability `exp(-delta / (2T))`. Lower T means stricter optimization (mostly downhill moves); higher T allows more exploration.

When annealing is **not** configured (i.e., `anneal_start` and `anneal_end` are omitted or equal), the simulation runs at a constant temperature of `anneal_end` (default **T = 0.1**). This is suitable for most refinements — it accepts small cost increases to avoid getting trapped in local minima while still driving chi2 down effectively.

| Temperature | Behavior |
|-------------|----------|
| T > 1 | Very permissive — accepts most moves, chi2 will typically increase |
| T ~ 0.1--1 | Moderate exploration with steady optimization |
| T < 0.1 | Greedy — mostly downhill, risk of trapping in local minima |

## Simulated annealing

For enhanced exploration, use simulated annealing with a cooling schedule:

```toml
anneal_start = 2.0       # Start hot (broad exploration)
anneal_end = 0.1         # Cool to this temperature
anneal_steps = 200_000   # Moves over which to anneal
```

The temperature follows an exponential schedule:

```
T(n) = T_start * (T_end / T_start) ^ (n / anneal_steps)
```

For n > anneal_steps, T stays at `anneal_end`. The remaining moves optimize at constant low temperature.

**Guidelines:**
- `anneal_start = 1.0--2.0` -- higher values explore more aggressively
- `anneal_end = 0.01--0.1` -- lower values optimize more strictly
- `anneal_steps` -- typically 10--30% of `max_moves`

## Adaptive step size

- **During annealing:** step scales with temperature, capped at `max_step`
- **After annealing:** step adapts to maintain the target acceptance rate. If acceptance > target + 5%, step grows by 5%; if below target - 5%, step shrinks by 5%
- Step is always clamped to [0.001, 2.0] A

## Convergence detection

When `convergence_threshold > 0`, the run monitors chi2 improvement over sliding windows of `convergence_window` moves. If the improvement is less than the threshold, the run stops early.

When annealing is active, convergence checking is deferred until the annealing phase completes. The first post-annealing window establishes a fresh baseline (no comparison), and actual convergence testing begins one window later. This prevents premature convergence triggered by the temporary chi2 increase during high-temperature annealing.

## Best-structure tracking

The lowest-chi2 configuration encountered during the entire run is tracked and restored at completion. This prevents the output from being worse than intermediate states (which can happen during annealing or random walks at high temperature).
