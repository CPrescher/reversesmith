# `[ml_potential]` -- ML-Potential-Regularized RMC

`[ml_potential]` enables an interatomic machine-learning potential as an RMC
energy regularizer:

```text
delta_cost = delta_chi2 + weight * delta_E_ml
```

This is RMC regularized by an ML potential, not EPSR refinement of the ML model.
In v1, `[ml_potential]` is only supported for normal RMC runs and cannot be
combined with `[epsr]` or `[potential]`.

## Native SNAP

```toml
[ml_potential]
backend = "snap_native"
coefficient_file = "potential.snapcoeff"
parameter_file = "potential.snapparam"
weight = 0.001
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `backend` | String | required | `snap_native` |
| `coefficient_file` | String | required | FitSNAP/LAMMPS `.snapcoeff` file, relative to the config file |
| `parameter_file` | String | required | Matching `.snapparam` file, relative to the config file |
| `weight` | Float | 0.001 | Scales the SNAP energy contribution in the RMC cost |

Do not set `cutoff`. rsmith derives every species-pair cutoff from `rcutfac`
and the element radii in the model files. A manually configured SNAP cutoff is
rejected so that it cannot silently omit affected atomic environments.

Native SNAP is part of the default Rust build. FitSNAP is used to fit the model,
but neither FitSNAP nor LAMMPS is needed when rsmith evaluates it. See
[Native SNAP](../snap.md) for setup, compatibility, and troubleshooting, and
[Native SNAP Implementation](../algorithms/snap-native.md) for the descriptor
math, local cache, and LAMMPS validation strategy.

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
device = "cpu"
torch_threads = 8
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `backend` | String | required | `mace_python` |
| `model` | String | required | MACE model path, relative to the config file |
| `weight` | Float | 0.001 | Scales the MACE energy contribution in the RMC cost |
| `cutoff` | Float | required | Model cutoff in A; retained for logging and future local-delta support |
| `device` | String | `cpu` | PyTorch device, e.g. `cpu`, `cuda`, or `mps` |
| `torch_threads` | Integer | PyTorch default | Sets `torch.set_num_threads()` in the worker |
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

For correctness, the first MACE backend computes:

```text
delta_E = E_full(after move) - E_full(before move)
```

This is slower than a local delta, but it is robust for message-passing models
where a moved atom can affect neighbors beyond one cutoff. Parallelism comes
from PyTorch inside each energy evaluation through `torch_threads` or the
selected accelerator device.

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
cargo test mace_python_backend_runs_real_model_when_configured -- --nocapture
```

If `RSMITH_MACE_TEST_MODEL` is not set, the real-model test skips itself.
