# Configuration Reference

Reversesmith uses a single TOML configuration file. All file paths are resolved relative to the config file location.

The configuration has the following sections:

| Section | Required | Description |
|---------|----------|-------------|
| [`[system]`](config/system.md) | Yes | Input structure and format |
| [`[data]`](config/data.md) | Yes | Experimental datasets to fit |
| [`[rmc]`](config/rmc.md) | Yes | Refinement parameters |
| [`[sq]`](config/sq.md) | No | S(Q) computation grid and cutoffs |
| [`[constraints]`](config/constraints.md) | No | Minimum distances and coordination bounds |
| [`[potential]`](config/potentials.md) | No | Pair potentials for hybrid RMC |
| [`[analysis]`](config/analysis.md) | No | Settings for `--analyze` mode |

See [Full Example](full-example.md) for a complete configuration file with all options.
