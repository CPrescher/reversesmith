#!/usr/bin/env bash
set -euo pipefail

case_dir="$(cd "$(dirname "$0")" && pwd)"
epsr_root="${EPSR26_ROOT:-}"

if [[ "${1:-}" != "--accept-local-testing-terms" ]]; then
    echo "Usage: $0 --accept-local-testing-terms" >&2
    echo "This imports EPSR tutorial files only for private local testing." >&2
    exit 2
fi
if [[ -z "$epsr_root" ]]; then
    echo "Set EPSR26_ROOT to the local EPSR installation directory." >&2
    exit 2
fi

source_dir="$epsr_root/workshop2022/EPSR_WorkedExamples_2022/DTBsilicaNX"
binary="$epsr_root/bin/epsr"
destination="$case_dir/reference/local"
upstream="$destination/upstream"

if [[ ! -d "$source_dir" ]]; then
    echo "DTBsilicaNX worked example not found at: $source_dir" >&2
    exit 1
fi
if [[ ! -x "$binary" ]]; then
    echo "EPSR executable not found or not executable at: $binary" >&2
    exit 1
fi

required=(
    DTBsilica.ato
    DTBsilica.EPSR.inp
    DTBsilica.pcof
    DTBsilica.NWTStot.wts
    DTBsilica.XWTS.wts
    SilicaGlassRT.mint01
    SiO2XRD.int01
)
for filename in "${required[@]}"; do
    if [[ ! -f "$source_dir/$filename" ]]; then
        echo "Required tutorial file is missing: $source_dir/$filename" >&2
        exit 1
    fi
done

mkdir -p "$destination"
if [[ -e "$upstream" ]]; then
    echo "Refusing to overwrite existing local import: $upstream" >&2
    echo "Move it aside first if a fresh import is intended." >&2
    exit 1
fi
cp -R "$source_dir" "$upstream"

(
    cd "$upstream"
    find . -type f -print0 | sort -z | xargs -0 shasum -a 256
) > "$destination/SHA256SUMS"

{
    echo "epsr_root=$epsr_root"
    echo "source_dir=$source_dir"
    echo "binary=$binary"
    echo "case=DTBsilicaNX"
    echo "local_testing_terms_accepted=true"
    echo "imported_utc=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
} > "$destination/IMPORT.txt"

echo "Imported private EPSR26 test files to $upstream"
echo "Checksums: $destination/SHA256SUMS"
echo "The imported files remain ignored and must not be redistributed."
