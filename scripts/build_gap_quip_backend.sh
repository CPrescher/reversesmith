#!/usr/bin/env sh
set -eu

repo_dir="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
prefix="${RSMITH_QUIP_PREFIX:-$HOME/Software/rsmith-gap-quip}"
default_quip_root="$HOME/Software/QUIP"

mkdir -p "$prefix/src" "$prefix/lib" "$prefix/lib/pkgconfig" "$prefix/obj"

ar="${AR:-ar}"

quip_root="${QUIP_ROOT:-}"
quip_build_dir="${QUIP_BUILD_DIR:-}"

if [ -z "$quip_root" ]; then
  quip_root="$default_quip_root"
  if [ ! -d "$quip_root/.git" ]; then
    mkdir -p "$(dirname "$quip_root")"
    git clone --recursive https://github.com/libAtoms/QUIP.git "$quip_root"
  fi
fi

if [ -z "$quip_build_dir" ]; then
  quip_arch="${QUIP_ARCH:-lammps}"
  if [ -z "$quip_arch" ]; then
    quip_arch="$(find "$quip_root/build" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | head -n 1 | xargs basename 2>/dev/null || true)"
  fi
  if [ ! -d "$quip_root/build/$quip_arch" ]; then
    quip_arch="$(find "$quip_root/build" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | head -n 1 | xargs basename 2>/dev/null || true)"
  fi
  if [ -z "$quip_arch" ]; then
    cat >&2 <<EOF
No QUIP build was found.

Set QUIP_BUILD_DIR to an existing QUIP build directory, for example:
  QUIP_BUILD_DIR=/path/to/QUIP/build/<arch> sh scripts/build_gap_quip_backend.sh

Or configure/build QUIP first:
  cd $quip_root
  make config
  make libquip
EOF
    exit 1
  fi
  quip_build_dir="$quip_root/build/$quip_arch"
fi

if [ ! -f "$quip_build_dir/libquip.a" ]; then
  cat >&2 <<EOF
QUIP build directory does not contain libquip.a:
  $quip_build_dir

Build QUIP first, or set QUIP_BUILD_DIR to the directory containing libquip.a.
EOF
  exit 1
fi

quip_config=""
search_dir="$quip_build_dir"
for _ in 1 2 3 4 5 6 7 8; do
  if [ -f "$search_dir/Makefile.inc" ]; then
    quip_config="$search_dir/Makefile.inc"
    break
  elif [ -f "$search_dir/quip.config" ]; then
    quip_config="$search_dir/quip.config"
    break
  fi
  parent="$(dirname "$search_dir")"
  if [ "$parent" = "$search_dir" ]; then
    break
  fi
  search_dir="$parent"
done

if [ -z "${CXX:-}" ] && [ -f "$quip_config" ]; then
  cxx="$(awk -F= '/^CPLUSPLUS[[:space:]]*=/ { print $2; exit }' "$quip_config" | xargs)"
  if [ -z "$cxx" ]; then
    cxx="c++"
  fi
else
  cxx="${CXX:-c++}"
fi

if [ -z "${FC:-}" ] && [ -f "$quip_config" ]; then
  fc="$(awk -F= '/^F95[[:space:]]*=/ { print $2; exit }' "$quip_config" | xargs)"
  if [ -z "$fc" ]; then
    fc="gfortran"
  fi
else
  fc="${FC:-gfortran}"
fi

if ! command -v "$cxx" >/dev/null 2>&1; then
  fc_major="$("$fc" -dumpversion 2>/dev/null | awk -F. '{ print $1 }')"
  matching_cxx="$(dirname "$fc")/g++-$fc_major"
  if [ -n "$fc_major" ] && [ -x "$matching_cxx" ]; then
    cxx="$matching_cxx"
  else
    cat >&2 <<EOF
The C++ compiler recorded by QUIP is unavailable:
  $cxx

Set CXX to the GNU C++ compiler matching your QUIP/GFortran installation.
EOF
    exit 1
  fi
fi

"$cxx" -std=c++17 -O2 -fPIC \
  -I"$repo_dir/include" \
  -c "$repo_dir/c_src/rsmith_gap_quip_lammps_shim.cpp" \
  -o "$prefix/obj/rsmith_gap_quip_lammps_shim.o"

"$fc" -O2 -fPIC -cpp -ffree-line-length-none \
  -I"$quip_build_dir" \
  -c "$repo_dir/c_src/rsmith_gap_quip_local_energy_wrapper.F90" \
  -o "$prefix/obj/rsmith_gap_quip_local_energy_wrapper.o"

"$ar" rcs "$prefix/lib/librsmith_quip_gap_shim.a" \
  "$prefix/obj/rsmith_gap_quip_lammps_shim.o" \
  "$prefix/obj/rsmith_gap_quip_local_energy_wrapper.o"

fox_dirs="$(find "$quip_root/src/fox" -path '*/lib' -type d 2>/dev/null | paste -sd: - || true)"
lib_dirs="$prefix/lib:$quip_build_dir"
if [ -n "$fox_dirs" ]; then
  lib_dirs="$lib_dirs:$fox_dirs"
fi
if [ -n "${EXTRA_QUIP_LIB_DIR:-}" ]; then
  lib_dirs="$lib_dirs:$EXTRA_QUIP_LIB_DIR"
fi

default_libs="${QUIP_LIBS:-rsmith_quip_gap_shim quip gap quip_core gapfit atoms quiputils FoX_wxml FoX_sax FoX_fsys FoX_wcml FoX_utils FoX_common FoX_wkml gfortran quadmath stdc++}"
if [ -n "${EXTRA_QUIP_LIBS:-}" ]; then
  default_libs="$default_libs $EXTRA_QUIP_LIBS"
fi

cat > "$prefix/env.sh" <<EOF
# Source this file before building or testing rsmith with --features gap-quip.
export QUIP_ROOT="$quip_root"
export QUIP_BUILD_DIR="$quip_build_dir"
export QUIP_INCLUDE_DIR="$quip_build_dir"
export QUIP_LIB_DIR="$lib_dirs"
export QUIP_LIBS="\${QUIP_LIBS:-$default_libs}"
export PKG_CONFIG_PATH="$prefix/lib/pkgconfig\${PKG_CONFIG_PATH:+:\$PKG_CONFIG_PATH}"
EOF

lib_flags=""
old_ifs=$IFS
IFS=:
for lib_dir in $lib_dirs; do
  if [ -n "$lib_dir" ]; then
    lib_flags="$lib_flags -L$lib_dir"
  fi
done
IFS=$old_ifs
for lib in $default_libs; do
  lib_flags="$lib_flags -l$lib"
done

cat > "$prefix/lib/pkgconfig/rsmith-gap-quip.pc" <<EOF
prefix=$prefix
quip_root=$quip_root
quip_build_dir=$quip_build_dir

Name: rsmith-gap-quip
Description: rsmith GAP/QUIP shim and QUIP link metadata
Version: 1.0
Cflags: -I\${quip_build_dir}
Libs:$lib_flags
EOF

cat <<EOF
Built rsmith GAP/QUIP shim:
  $prefix/lib/librsmith_quip_gap_shim.a

Wrote environment:
  $prefix/env.sh

Wrote pkg-config metadata:
  $prefix/lib/pkgconfig/rsmith-gap-quip.pc

Use it with:
  . "$prefix/env.sh"
  cargo build --release --features gap-quip

Or with pkg-config:
  export PKG_CONFIG_PATH="$prefix/lib/pkgconfig\${PKG_CONFIG_PATH:+:\$PKG_CONFIG_PATH}"
  cargo build --release --features gap-quip

If linking fails on BLAS/LAPACK symbols, append your BLAS/LAPACK library to QUIP_LIBS
and add its directory to QUIP_LIB_DIR. For example, with Homebrew OpenBLAS:
  export QUIP_LIB_DIR="\$QUIP_LIB_DIR:/opt/homebrew/opt/openblas/lib"
  export QUIP_LIBS="rsmith_quip_gap_shim quip gap quip_core gapfit atoms quiputils FoX_wxml FoX_sax FoX_fsys FoX_wcml FoX_utils FoX_common FoX_wkml openblas gfortran quadmath stdc++"

You can also pass EXTRA_QUIP_LIB_DIR and EXTRA_QUIP_LIBS to this script to bake
site-specific library paths into env.sh.
EOF
