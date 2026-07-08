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
