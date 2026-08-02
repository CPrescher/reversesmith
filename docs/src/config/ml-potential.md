# `[ml_potential]` -- ML-Potential-Regularized RMC

`[ml_potential]` enables an interatomic machine-learning potential as an RMC
energy regularizer:

```text
delta_cost = delta_chi2 + weight * delta_E_ml
```

This is RMC regularized by an ML potential, not EPSR refinement of the ML model.
In v1, `[ml_potential]` is only supported for normal RMC runs and cannot be
combined with `[epsr]` or `[potential]`.

## GAP/QUIP

```toml
[ml_potential]
backend = "gap_quip"
model = "gap.xml"
init_args = "Potential xml_label=GAP_2026"
weight = 0.001
cutoff = 5.0
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `backend` | String | required | `gap_quip` |
| `model` | String | required | GAP XML model path, relative to the config file |
| `init_args` | String | none | Optional QUIP initialization string, e.g. an `xml_label` |
| `weight` | Float | 0.001 | Scales the ML energy contribution in the RMC cost |
| `cutoff` | Float | required | Local environment cutoff in A |

`cutoff` must be at least the maximum descriptor cutoff used by the GAP model.
If it is too small, local energy deltas will be wrong because some changed
atomic environments will be missed.

## Building With QUIP

The default build does not require QUIP. The `gap-quip` feature links to an
external QUIP/libAtoms build plus the small rsmith C ABI shim.

Recommended workflow:

```bash
sh scripts/build_gap_quip_backend.sh
. $HOME/Software/rsmith-gap-quip/env.sh
cargo build --features gap-quip
```

The helper writes both:

```text
$HOME/Software/rsmith-gap-quip/env.sh
$HOME/Software/rsmith-gap-quip/lib/pkgconfig/rsmith-gap-quip.pc
```

The env file sets `QUIP_INCLUDE_DIR`, `QUIP_LIB_DIR`, `QUIP_LIBS`, and
`PKG_CONFIG_PATH`. If the explicit `QUIP_*` variables are absent, rsmith's
build script tries `pkg-config rsmith-gap-quip`.

For manual or cluster builds:

```bash
QUIP_INCLUDE_DIR=/path/to/quip/include \
QUIP_LIB_DIR=/path/to/quip/lib \
QUIP_LIBS="rsmith_quip_gap_shim quip gap quip_core gapfit atoms quiputils ..." \
cargo build --features gap-quip
```

QUIP must be built separately on each platform with a consistent Fortran/C/C++
toolchain and GAP support enabled. Static QUIP builds often need BLAS/LAPACK
and Fortran runtime libraries in `QUIP_LIBS`.

The `gap-quip` feature expects a small C ABI shim that owns the QUIP calculator
and exposes per-atom GAP energies. The Rust side calls the following symbols:

```c
rsmith_gap_quip_create
rsmith_gap_quip_destroy
rsmith_gap_quip_set_structure
rsmith_gap_quip_move_atom
rsmith_gap_quip_per_atom_energy
```

`rsmith_gap_quip_create()` receives both the model path and the optional
`init_args` string. When `init_args` is omitted, the shim receives a null
pointer and should use its own safe default for the supplied GAP XML.

The backend fails at startup if GAP/QUIP support is requested but the binary was
not built with `--features gap-quip`.

## MACE/Python

MACE support uses a small Python worker process. The Rust binary stays pure
Rust, while the worker imports `mace-torch`, PyTorch, and ASE from the Python
environment selected in the config.

```toml
[ml_potential]
backend = "mace_python"
model = "mace.model"
weight = 0.001
cutoff = 5.0
delta = "full"
device = "cpu"
torch_threads = 8
dtype = "float32"
compile_mode = "reduce_overhead"
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `backend` | String | required | `mace_python` |
| `model` | String | required | MACE model path, relative to the config file |
| `weight` | Float | 0.001 | Scales the MACE energy contribution in the RMC cost |
| `cutoff` | Float | required | Model cutoff in A; should match the MACE model cutoff |
| `delta` | String | `full` | Energy delta mode: `full`, `local`, or `incremental` |
| `device` | String | `cpu` | PyTorch device, e.g. `cpu`, `cuda`, or `mps` |
| `torch_threads` | Integer | PyTorch default | Sets `torch.set_num_threads()` in the worker |
| `dtype` | String | Checkpoint dtype | Optional `float32` or `float64` model conversion |
| `compile_mode` | String | disabled | Optional PyTorch mode: `default`, `reduce_overhead`, or `max_autotune` |
| `python` | String | `python3` | Python executable used to launch the worker |
| `worker` | String | embedded worker | Optional custom worker script for testing or site integration |

Create the optional MACE Python environment with uv:

```bash
uv sync --group mace
cargo build --release
```

Then point rsmith at `.venv/bin/python`, or set `python = ".venv/bin/python"` in
`[ml_potential]`. No Rust feature flag is required for `mace_python` because
MACE is loaded by the Python worker at runtime.

The default MACE backend computes:

```text
delta_E = E_full(after move) - E_full(before move)
```

This is slower, but it is robust for message-passing models where a moved atom
can affect neighbors beyond one cutoff. Parallelism comes from PyTorch inside
each energy evaluation through `torch_threads` or the selected accelerator
device.

For short-range MACE models, `delta = "local"` evaluates a local non-periodic
cluster with explicit periodic image atoms as context. It sums MACE per-atom
energies only for atoms within `num_interactions * cutoff` of the moved atom
before or after the trial move. Context atoms are included out to the same
message-passing radius around those central atoms. At fixed density this makes
the cost effectively independent of total atom count once the simulation cell
is larger than the cluster. A dedicated Rust cell list generates compact
`(atom index, periodic image shift)` arrays and places the affected central
atoms in one coherent image around the moved atom; the Python worker only
assembles and evaluates the resulting cluster. The worker caches accepted
full-system per-atom energies, so a trial evaluates only the proposed cluster
once, without force/stress autograd. Accepting the trial commits its affected
per-atom energies to the cache; rejecting it restores the position and drops
the pending update.

For standard short-range `MACE` and `ScaleShiftMACE` checkpoints,
`delta = "incremental"` caches the accepted hidden node features at every
interaction layer. It recomputes only the moved atom's layer-by-layer causal
cone and uses cached features for unchanged source atoms. Its context cluster
therefore needs only one cutoff beyond the final affected atoms. The result is
an exact message-passing update for supported models, not a reduced-cutoff
approximation. Incremental mode currently rejects pair-repulsion terms, joint
embeddings, fused interaction kernels, and `compile_mode`.

`dtype = "float32"` converts the checkpoint for faster inference and must be
validated against the desired numerical accuracy. `compile_mode` maps to
PyTorch's `default`, `reduce-overhead`, or `max-autotune` modes; compilation can
increase initialization time and is not guaranteed to improve small local
graphs.

Use `delta = "local"` and `delta = "incremental"` only for ordinary short-range
MACE models. Do not use them for models with explicit long-range electrostatics,
global charge/spin coupling, or dispersion corrections unless those terms are
separately validated. The real MACE integration test compares local,
incremental, and full deltas when a test model is provided.

MACE code and MACE model files have separate licenses. Check the license of the
specific model file before redistribution or publication.

### Real MACE Integration Test

The normal test suite uses a mock worker and does not require PyTorch or MACE.
To run against an actual MACE checkpoint, create the uv environment and download
the small MIT-licensed MACE-MP model:

```bash
uv sync --group mace
MODEL=$(.venv/bin/python scripts/download_mace_test_model.py)
RSMITH_MACE_TEST_MODEL="$MODEL" \
RSMITH_MACE_TEST_PYTHON=.venv/bin/python \
RSMITH_MACE_TEST_DEVICE=cpu \
RSMITH_MACE_TEST_TORCH_THREADS=1 \
cargo test mace_python -- --nocapture
```

If `RSMITH_MACE_TEST_MODEL` is not set, the real-model tests skip themselves.
