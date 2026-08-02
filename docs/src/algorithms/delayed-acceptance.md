# Delayed Acceptance

Delayed acceptance reduces the average cost of hybrid RMC when evaluating the
configured pair or machine-learning potential is much more expensive than
updating the experimental residual. It screens a proposal with the data term
before calling the energy backend.

Enable it in the RMC section:

```toml
[rmc]
delayed_acceptance = true
energy_calibration_moves = 1000
```

The implementation is shared by all energy backends: analytical and tabulated
pair potentials, native SNAP, GAP/QUIP, and MACE/Python. It is normally most
useful for GAP and MACE because those evaluations dominate the cost of an RMC
move.

## Standard and delayed decisions

Standard hybrid RMC uses the cost

```text
C(X) = chi2(X) + weight * E(X)
```

and accepts a proposed configuration `Y` from the accepted configuration `X`
with probability

```text
min(1, exp(-(delta_chi2 + weight * delta_E) / (2T))).
```

Delayed acceptance replaces this one decision with two consecutive decisions:

1. Check hard structural constraints. Constraint failures are rejected as
   usual.
2. Calculate `delta_chi2` using the incremental scattering update.
3. Apply the data-stage probability
   `min(1, exp(-delta_chi2 / (2T)))`.
4. If the data stage rejects, restore the proposed coordinate and skip the
   energy backend completely.
5. If it passes, calculate the exact configured `delta_E` and apply
   `min(1, exp(-weight * delta_E / (2T)))`.
6. Accept the move only if both stages pass.

The same temperature is used for both stages, including during simulated
annealing. A data-improving move (`delta_chi2 <= 0`) always reaches the energy
stage. An energy-improving move (`weight * delta_E <= 0`) always passes the
second stage once evaluated.

## Why the target distribution is unchanged

For a Metropolis factor `a(delta) = min(1, exp(-delta / (2T)))`, the ratio of
the forward and reverse probabilities is

```text
a(delta) / a(-delta) = exp(-delta / (2T)).
```

The two delayed factors therefore give

```text
[a(delta_chi2) / a(-delta_chi2)]
* [a(weight * delta_E) / a(-weight * delta_E)]
= exp(-(delta_chi2 + weight * delta_E) / (2T)).
```

This is the same detailed-balance ratio as standard hybrid RMC. The potential
is not approximated, truncated, or replaced: every proposal that reaches stage
two uses the exact energy-delta mode configured for that backend.

Delayed acceptance does change the transition kernel. When `delta_chi2` and
`weight * delta_E` have opposite signs, the combined rule can allow one term to
compensate for the other before testing. The delayed rule tests them separately
and can reject more often. It consequently has the same target distribution
but may mix more slowly.

## Expected performance

Let `p_data` be the fraction of constraint-valid proposals that pass the data
stage. The approximate average cost is

```text
cost_per_attempt ~= data_cost + p_data * potential_cost.
```

For example, with a 0.4 ms data update and a 21 ms local GAP evaluation, a 30%
data-stage pass rate gives approximately

```text
0.4 ms + 0.30 * 21 ms = 6.7 ms per attempted move.
```

This is only a throughput estimate. Judge the feature using accepted moves per
second and, for production sampling, independent structural configurations per
unit time. A large skip percentage is not beneficial if the lower acceptance
rate produces substantially worse mixing.

The progress and final log messages report `potential skip`, the fraction of
data-stage trials that avoided an energy evaluation. Proposals rejected by hard
constraints are not included in that denominator.

## Backend guidance

| Backend | Recommendation |
|---------|----------------|
| Pair potentials | Usually leave disabled. Their local interpolation cost is often smaller than the possible loss in acceptance. Benchmark if the cutoff or neighbor count is unusually large. |
| SNAP | Supported. It can help when the data-stage rejection rate is high, but native local SNAP is already inexpensive, so measure accepted moves per second. |
| GAP `local` | Recommended candidate. The cached accepted energies remain untouched after a stage-one rejection, and only proposals passing the data test enter QUIP. |
| GAP `full` | Supported and potentially valuable because a skipped proposal avoids both full-system evaluations. Prefer `local` when the model is compatible. |
| MACE `incremental` | Recommended candidate. A stage-one rejection avoids graph construction and partial message-passing inference. |
| MACE `local` or `full` | Supported and potentially valuable; combine delayed acceptance with the fastest delta mode validated for the model. |

Delayed acceptance and a backend's `delta` option solve different problems.
`delta = "local"` or `"incremental"` reduces the cost of one energy evaluation;
delayed acceptance reduces how often that evaluation is needed. They can and
usually should be combined for expensive compatible models.

## Energy-weight calibration

The existing energy-weight calibration estimates representative
`abs(delta_chi2)` and `abs(delta_E)` values. Sampling `delta_E` only after a
data-stage pass would bias that estimate, so rsmith performs
`energy_calibration_moves` complete, unscreened evaluations first. Delayed
screening begins after the calibration window.

The default is 1000. If a production energy weight has already been established
with the same data normalization, model, structure class, and temperature
scheme, start screening immediately with:

```toml
[rmc]
delayed_acceptance = true
energy_calibration_moves = 0
```

With zero calibration moves, rsmith does not print a new weight-calibration
estimate for that run.

## Reproducibility and limitations

- Enabling delayed acceptance changes the random-number sequence and accepted
  trajectory, even when the seed is unchanged.
- Checkpoints store the accepted structure and RMC state, not a promise to
  reproduce the trajectory of a run using a different acceptance setting.
- The feature applies to hybrid RMC with an energy model. It has no effect when
  no pair or ML potential is configured.
- Pure energy-only EPSR/MC has no experimental first stage and therefore cannot
  use this screening strategy. Hybrid RMC portions can use it.
- rsmith currently implements the data-only first stage. A separate cheap
  surrogate-potential first stage is not implemented.

For the configuration fields, see [RMC Parameters](../config/rmc.md). For
backend-specific energy-delta modes, see [ML Potentials](../config/ml-potential.md).
