# `[analysis]` -- Structural Analysis Settings

Controls the `--analyze` mode for computing coordination numbers and bond angle distributions.

```toml
[analysis]
angle_bins = 180                    # Bins for 0-180 degree range (default: 180)
angle_triplets = ["O-Si-O"]        # Only compute these triplets (default: all)

[analysis.cutoffs]
"Si-O" = 2.2
"Ca-O" = 3.0
```

## Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `cutoffs` | Map | -- | Pair cutoffs for CN and angle analysis (A) |
| `angle_bins` | Integer | 180 | Number of bins for angle histograms (0--180 degrees) |
| `angle_triplets` | List | all | Filter: only compute these angle triplets |

## Cutoff resolution

If `[analysis.cutoffs]` is not specified, cutoffs from `[[constraints.coordination]]` are used as fallback. Both sources are merged, with explicit `[analysis.cutoffs]` taking precedence.

This means if you already have coordination constraints, `--analyze` works without any additional configuration:

```toml
[[constraints.coordination]]
pair = "Si-O"
min = 3
max = 6
cutoff = 2.2
# This cutoff is automatically used by --analyze
```

## Angle triplets

The `angle_triplets` filter uses canonical labels with **sorted end-species**: `"O-Si-O"`, not `"Si-O-O"`. The central atom is always in the middle.

See the [Structural Analysis](../analysis.md) chapter for full usage details and interpretation.
