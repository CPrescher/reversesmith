# rsmith

A Reverse Monte Carlo (RMC) structure refinement tool written in Rust. Refines atomic structures against experimental X-ray and neutron scattering data while enforcing physical constraints, pair potentials, and optional machine-learning potential regularizers.

## Features

- **S(Q) and g(r) fitting** -- simultaneous refinement against structure factor and pair distribution function
- **Hybrid RMC** -- pair potentials (Buckingham, Pedone, Coulomb DSF, tabulated) bias refinement toward energetically favorable configurations
- **ML-potential-regularized RMC** -- optional native SNAP, GAP/QUIP, or MACE energy regularizers
- **Hard constraints** -- minimum interatomic distances and coordination number bounds
- **Simulated annealing** -- exponential cooling schedule with adaptive step size
- **Structural analysis** -- coordination numbers and bond angle distributions for validation
- **Incremental updates** -- O(N_neighbors) per move via cell lists and precomputed lookup tables

## Quick Start

```bash
cargo build --release
./target/release/rsmith config.toml
```

Minimal `config.toml`:

```toml
[system]
structure = "glass.data"
format = "lammps"

[system.types]
1 = "Ca"
2 = "Si"
3 = "O"

[data.xray_sq]
file = "experimental.sq"
sigma = 0.02

[rmc]
max_moves = 500_000

[sq]
rdf_cutoff = 11.0

[constraints.min_distance]
"Si-O" = 1.2
"O-O" = 2.0
"Ca-O" = 1.8
```

## Modes

```bash
rsmith config.toml                      # RMC refinement
rsmith config.toml --compute-sq-only    # Compute S(Q) only
rsmith config.toml --analyze            # Structural analysis
```

## Documentation

Full documentation is available at [rsmith book](docs/src/SUMMARY.md) or build locally:

```bash
cd docs && mdbook serve
```

## Building With GAP/QUIP

GAP support is intentionally built through a small C ABI shim:

```text
rsmith Rust code -> rsmith GAP/QUIP shim -> QUIP/libAtoms
```

Do not bind Rust directly to QUIP internals. QUIP/libAtoms is primarily a
Fortran library, and the stable boundary for rsmith is the C ABI declared in
[`include/rsmith_gap_quip_shim.h`](include/rsmith_gap_quip_shim.h). The shim is
responsible for owning the QUIP calculator, keeping an atom structure
synchronized with rsmith, and returning summed per-atom GAP energies for a list
of atom indices.

The shim must export these symbols:

```c
rsmith_gap_quip_create
rsmith_gap_quip_destroy
rsmith_gap_quip_set_structure
rsmith_gap_quip_move_atom
rsmith_gap_quip_per_atom_energy
```

`rsmith_gap_quip_per_atom_energy()` is required. A total-energy-only wrapper is
not sufficient for production RMC, because rsmith computes local GAP energy
deltas from the moved atom plus atoms whose local environments include it.

### External Build Layout

QUIP source and build outputs should be kept outside this repository as a
standalone dependency. The default layout used by the helper is:

```text
$HOME/Software/QUIP/
  build/lammps/      # QUIP build directory containing libquip.a

$HOME/Software/rsmith-gap-quip/
  lib/               # compiled rsmith shim library
  obj/               # shim object files
  env.sh             # environment variables for cargo
```

Use `QUIP_ROOT=/some/QUIP` for a different QUIP checkout and
`RSMITH_QUIP_PREFIX=/some/path` for a different shim/env location.
Do not vendor QUIP into the rsmith repository.

### 1. Build QUIP

Build QUIP with GAP support and produce a linkable QUIP library. For a standard
QUIP checkout this is typically:

```bash
git clone --recursive https://github.com/libAtoms/QUIP.git
cd QUIP
make config
make libquip
```

On the local macOS setup used for this repository, QUIP is built at
`$HOME/Software/QUIP` with `QUIP_ARCH=lammps`, Homebrew GCC/GFortran 15, and
GAP enabled. That produces:

```text
$HOME/Software/QUIP/build/lammps/libquip.a
```

Use the compiler/MPI/OpenMP/BLAS/LAPACK setup appropriate for your machine. Keep
the chosen QUIP architecture name and library directory; they are needed when
building the shim and rsmith.

### 2. Build the rsmith GAP/QUIP shim outside the repo

The helper script builds only the rsmith shim into the external prefix. It can
use either the default `$HOME/Software/QUIP` build or an explicit `QUIP_ROOT`.

For the default local layout:

```bash
RSMITH_QUIP_PREFIX=$HOME/Software/rsmith-gap-quip \
QUIP_ROOT=$HOME/Software/QUIP \
QUIP_BUILD_DIR=$HOME/Software/QUIP/build/lammps \
CXX=/opt/homebrew/bin/g++-15 \
EXTRA_QUIP_LIB_DIR=/opt/homebrew/opt/openblas/lib:/opt/homebrew/opt/gcc/lib/gcc/current \
EXTRA_QUIP_LIBS="openblas" \
sh scripts/build_gap_quip_backend.sh
```

`EXTRA_QUIP_LIB_DIR` and `EXTRA_QUIP_LIBS` are optional. Use them when your QUIP
archive depends on libraries outside the QUIP build directory, such as
Homebrew OpenBLAS or a LAMMPS-generated linear algebra archive.

The shim should be compiled with the same C++ toolchain used by the QUIP build.
The helper tries to read that from `build/<arch>/Makefile.inc` or `quip.config`;
if that is not available, set `CXX` explicitly, for example
`CXX=/opt/homebrew/bin/g++-15`.

For a fresh external QUIP checkout, omit `QUIP_ROOT`; the script downloads QUIP
into `$HOME/Software/QUIP`. You still need to configure/build QUIP if no
`libquip.a` exists yet:

```bash
RSMITH_QUIP_PREFIX=$HOME/Software/rsmith-gap-quip \
sh scripts/build_gap_quip_backend.sh
```

The script writes:

```text
$RSMITH_QUIP_PREFIX/lib/librsmith_quip_gap_shim.a
$RSMITH_QUIP_PREFIX/env.sh
$RSMITH_QUIP_PREFIX/lib/pkgconfig/rsmith-gap-quip.pc
```

The shim includes `include/rsmith_gap_quip_shim.h` and calls QUIP's
`quip_lammps_wrapper`, because that wrapper already returns per-atom
`local_energy`. Its implementation:

1. Load the GAP XML model and optional QUIP initialization string in `rsmith_gap_quip_create()`.
2. Build/update the QUIP atoms object in `rsmith_gap_quip_set_structure()`.
3. Update one atom position in `rsmith_gap_quip_move_atom()`.
4. Recompute per-atom energies for the current structure.
5. Sum only the requested atom indices in `rsmith_gap_quip_per_atom_energy()`.

Return `0` on success from mutating/evaluation functions and nonzero on error.
Return a null handle from `rsmith_gap_quip_create()` if initialization fails.

### 3. Build rsmith with the GAP backend

The recommended build path is to source the generated environment. This sets
both the explicit `QUIP_*` variables and `PKG_CONFIG_PATH`:

```bash
. $HOME/Software/rsmith-gap-quip/env.sh
cargo build --release --features gap-quip
```

Alternatively, use the generated pkg-config metadata directly:

```bash
export PKG_CONFIG_PATH=$HOME/Software/rsmith-gap-quip/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}
cargo build --release --features gap-quip
```

When `QUIP_INCLUDE_DIR` and `QUIP_LIB_DIR` are set, they take precedence over
pkg-config. This is useful for cluster module systems where link flags are
managed explicitly:

```bash
export QUIP_INCLUDE_DIR=/path/to/QUIP/build/<arch>
export QUIP_LIB_DIR=/path/to/rsmith-gap-quip/lib:/path/to/QUIP/build/<arch>:/path/to/extra/libs
export QUIP_LIBS="rsmith_quip_gap_shim quip gap quip_core gapfit atoms quiputils FoX_wxml FoX_sax FoX_fsys FoX_wcml FoX_utils FoX_common FoX_wkml gfortran quadmath stdc++ openblas"
cargo build --release --features gap-quip
```

`QUIP_LIB_DIR` uses the platform path separator (`:` on Linux/macOS, `;` on
Windows). `QUIP_LIBS` is split on spaces or commas and passed to the linker as
`-l` entries in order. Adjust it for your QUIP build. Static QUIP builds often
need Fortran runtime and BLAS/LAPACK libraries explicitly; shared-library QUIP
builds may need fewer entries.

### Platform Notes

The default rsmith build is pure Rust and does not require QUIP. The
`gap-quip` feature is a native build and must be configured per platform.

macOS:

```bash
brew install gcc openblas pkg-config
cd $HOME/Software/QUIP
QUIP_ARCH=lammps QUIP_ROOT=$PWD make libquip
cd /path/to/rsmith
EXTRA_QUIP_LIB_DIR=/opt/homebrew/opt/openblas/lib:/opt/homebrew/opt/gcc/lib/gcc/current \
EXTRA_QUIP_LIBS=openblas \
sh scripts/build_gap_quip_backend.sh
. $HOME/Software/rsmith-gap-quip/env.sh
cargo build --release --features gap-quip
```

Linux:

```bash
# Package names vary by distribution.
# Install gcc, g++, gfortran, make, pkg-config, and BLAS/LAPACK/OpenBLAS.
cd $HOME/Software/QUIP
make config
make libquip
cd /path/to/rsmith
QUIP_ROOT=$HOME/Software/QUIP \
QUIP_BUILD_DIR=$HOME/Software/QUIP/build/<arch> \
EXTRA_QUIP_LIBS="openblas" \
sh scripts/build_gap_quip_backend.sh
. $HOME/Software/rsmith-gap-quip/env.sh
cargo build --release --features gap-quip
```

Clusters/modules:

```bash
module load gcc openblas pkg-config
export QUIP_ROOT=/path/to/site/QUIP
export QUIP_BUILD_DIR=/path/to/site/QUIP/build/<arch>
export RSMITH_QUIP_PREFIX=$HOME/Software/rsmith-gap-quip
sh scripts/build_gap_quip_backend.sh
. $RSMITH_QUIP_PREFIX/env.sh
cargo build --release --features gap-quip
```

Windows is not currently a supported GAP/QUIP target. It may be possible with a
consistent Fortran/C++ toolchain, but it needs separate testing.

At runtime, the dynamic loader must be able to find both the shim and QUIP
libraries if any of them are shared libraries. For example:

```bash
# Linux
LD_LIBRARY_PATH=/path/to/quip-and-shim-libs ./target/release/rsmith config.toml

# macOS
DYLD_LIBRARY_PATH=/path/to/quip-and-shim-libs ./target/release/rsmith config.toml
```

Using an installed rpath is preferable for cluster/module deployments.

### 4. Configure an ML-potential RMC run

`[ml_potential]` is an alternative to `[potential]` for normal RMC. In the first
implementation it cannot be combined with `[potential]` or `[epsr]`.

For a linear SNAP potential fitted by FitSNAP, rsmith can evaluate the model
natively without installing or launching LAMMPS:

```toml
[ml_potential]
backend = "snap_native"
coefficient_file = "potential.snapcoeff"
parameter_file = "potential.snapparam"
weight = 0.001
```

The model files use the standard LAMMPS/FitSNAP `.snapcoeff` and `.snapparam`
formats. rsmith derives the species-pair cutoffs from `rcutfac` and the element
radii, maintains its own SNAP-sized cell list, and caches accepted per-atom
energies. A trial atom move therefore recomputes only the local environments
within the old or new cutoff shells. The configured cell must be larger than
twice the largest model cutoff in every direction so that the minimum-image
environment is unambiguous. Do not set `cutoff` for `snap_native`; unlike the
GAP and MACE backends, the native SNAP cutoff is derived from the model files,
and a manual value is rejected to avoid silently using the wrong range.

Native evaluation supports linear standard and explicit multi-element
(`chemflag`) SNAP models, including `switchflag`, `rmin0`, per-element `wj` and
`radelem`, `bzeroflag`, `bnormflag`, both `wselfallflag` conventions, and
quadratic (`quadraticflag`) energy models. FitSNAP remains responsible for
fitting; rsmith only loads and evaluates the resulting potential. LAMMPS is
used as a numerical test oracle and is not a runtime dependency.

Small reproducible model files, example cells, and frozen LAMMPS reference
energies live in `tests/data/snap`. Run their native checks with:

```bash
cargo test --test snap_model_files --test snap_reference_fixture
```

For GAP/QUIP, use:

```toml
[ml_potential]
backend = "gap_quip"
model = "gap.xml"
init_args = "Potential xml_label=GAP_2026"
weight = 0.001
cutoff = 5.0
```

`cutoff` must be at least the maximum descriptor cutoff used by the GAP model.
If it is too small, rsmith will miss affected local environments and the local
energy delta will be wrong.

`init_args` is optional. Use it when your GAP XML requires a specific QUIP
initialization string, such as an `xml_label`. If omitted, the shim receives a
null pointer and must choose an appropriate default for the model.

For MACE, use the Python worker backend:

```toml
[ml_potential]
backend = "mace_python"
model = "mace.model"
weight = 0.001
cutoff = 5.0
delta = "full"       # full, local, or incremental
device = "cpu"        # cpu, cuda, or mps
torch_threads = 8     # optional CPU threading control
dtype = "float32"     # optional; defaults to the checkpoint dtype
compile_mode = "reduce_overhead" # optional torch.compile mode
```

Create the optional MACE Python environment with uv:

```bash
uv sync --group mace
cargo build --release
```

`mace_python` does not require a Rust feature flag. rsmith launches an embedded
Python worker and keeps the MACE calculator alive across trial moves. By
default it computes full-system energy differences for correctness:

```text
delta_E = E_full(after move) - E_full(before move)
```

For ordinary short-range MACE models, `delta = "local"` evaluates a bounded
local cluster with explicit periodic image atoms as context and sums only the
affected per-atom energies. rsmith builds the cluster with a dedicated Rust
cell list, places affected atoms in one coherent periodic image around the
moved atom, and sends compact atom-index/image-shift arrays to the Python
worker. The worker caches the accepted full-system per-atom energies, so each
trial needs only one energy-only MACE forward pass for the proposed local
cluster; accepting a move updates the affected cache entries and rejecting it
discards them. The affected radius is `num_interactions * cutoff`, so at fixed
density the trial cost becomes effectively independent of total atom count once
the simulation cell is larger than the local cluster.

For standard short-range `MACE` and `ScaleShiftMACE` checkpoints,
`delta = "incremental"` additionally caches the accepted hidden features after
every interaction layer. A trial recomputes the first layer only for the moved
atom and its immediate neighbours, then expands the dirty set by one neighbour
shell per layer. Unchanged source features come from the accepted cache. This
is exact message-passing inference rather than a truncated approximation.
Incremental mode currently rejects models with pair-repulsion terms, joint
embeddings, fused interaction kernels, or `compile_mode`.

`dtype = "float32"` can materially improve CPU throughput, but changes the
checkpoint's numerical precision and should be validated for the intended
model. `compile_mode` accepts `default`, `reduce_overhead`, or `max_autotune`;
whether compilation helps is hardware- and workload-dependent.

Use local and incremental MACE deltas only for short-range models. Keep
`delta = "full"` for models with explicit long-range electrostatics, global
charge/spin coupling, or dispersion corrections unless the alternative delta
has been separately validated.
MACE code and MACE model files have separate licenses; check the license of the
specific model before redistribution.

To run the real MACE integration test with the small MIT-licensed MACE-MP model:

```bash
MODEL=$(.venv/bin/python scripts/download_mace_test_model.py)
RSMITH_MACE_TEST_MODEL="$MODEL" \
RSMITH_MACE_TEST_PYTHON=.venv/bin/python \
RSMITH_MACE_TEST_DEVICE=cpu \
RSMITH_MACE_TEST_TORCH_THREADS=1 \
cargo test mace_python -- --nocapture
```

To measure CPU thread scaling through the same rejected single-atom trial path
used by RMC, run the release-mode benchmark example:

```bash
MODEL=$(.venv/bin/python scripts/download_mace_test_model.py)
cargo run --release --example mace_scaling -- \
  --model "$MODEL" \
  --python .venv/bin/python \
  --delta local
```

By default this benchmarks 216-, 1,000-, and 1,728-atom diamond-Si cells with
1, 2, 4, 8, 10, and 14 PyTorch threads. It prints CSV containing per-trial
statistics, speedup relative to one thread, and parallel efficiency. Use
`--delta full`, `--delta local`, or `--delta incremental` to select the
energy-delta strategy, and use `--help` to change the system sizes, thread
counts, precision, compilation mode, warm-up, or sample count.

### 5. Smoke-test the binding

If you have a small GAP XML test model, the feature-gated integration test can
check that the backend initializes:

```bash
. $HOME/Software/rsmith-gap-quip/env.sh
RUST_MIN_STACK=134217728 \
RSMITH_GAP_TEST_MODEL=/path/to/gap.xml \
RSMITH_GAP_TEST_INIT_ARGS="Potential xml_label=GAP_LABEL" \
cargo test --features gap-quip gap_quip_backend_initializes_when_test_model_is_available
```

If the test cannot find `RSMITH_GAP_TEST_MODEL`, it skips itself. `RUST_MIN_STACK`
is recommended because QUIP/GAP initialization can use more stack than Rust's
default test-thread stack.

## References

- McGreevy, R.L. & Pusztai, L. (1988). Reverse Monte Carlo Simulation. *Mol. Simul.*, 1, 359-367.
- Pedone, A. et al. (2006). A New Self-Consistent Empirical Interatomic Potential Model for Oxides. *J. Phys. Chem. B*, 110, 11780-11795.
