#!/usr/bin/env bash
set -euo pipefail

case_root="$(cd "$(dirname "$0")" && pwd)"
source_root="${SIO2_GLASS_ROOT:-}"
force=0

usage() {
    echo "Usage: SIO2_GLASS_ROOT=/path/to/SiO2_glass $0 [--force]" >&2
}

for argument in "$@"; do
    case "$argument" in
        --force) force=1 ;;
        -h|--help) usage; exit 0 ;;
        *) usage; exit 2 ;;
    esac
done

if [[ -z "$source_root" || ! -d "$source_root/.git" ]]; then
    echo "SIO2_GLASS_ROOT must name a local checkout of CPrescher/SiO2_glass" >&2
    exit 2
fi

gap_relative="results/production/staged/data/sio2_gap_from_pedone_glass_300K.data"
pedone_relative="results/production/pedone-reference/data/sio2_pedone_fast_glass_300K_charge.data"
rdf_relative="analysis/results/rdf_gap_vs_pedone.csv"
for relative in "$gap_relative" "$pedone_relative" "$rdf_relative"; do
    if [[ ! -f "$source_root/$relative" ]]; then
        echo "Missing required structure: $source_root/$relative" >&2
        exit 2
    fi
done

target="$case_root/reference/local/ambient-models"
if [[ -e "$target" && "$force" -ne 1 ]]; then
    echo "Target already exists: $target (pass --force to replace its files)" >&2
    exit 2
fi
mkdir -p "$target"
cp "$source_root/$gap_relative" "$target/gap-300K.data"
cp "$source_root/$pedone_relative" "$target/pedone-300K.data"
cp "$source_root/$rdf_relative" "$target/upstream-partial-rdf.csv"

commit="$(git -C "$source_root" rev-parse HEAD)"
remote="$(git -C "$source_root" remote get-url origin 2>/dev/null || true)"
{
    echo "repository=CPrescher/SiO2_glass"
    echo "remote=$remote"
    echo "commit=$commit"
    echo "source_root=$source_root"
    echo "gap_source=$gap_relative"
    echo "pedone_source=$pedone_relative"
    echo "rdf_source=$rdf_relative"
    echo "redistribution=local-only-until-repository-license-is-explicit"
} > "$target/IMPORT.txt"

shasum -a 256 "$target/gap-300K.data" "$target/pedone-300K.data" "$target/upstream-partial-rdf.csv" > "$target/SHA256SUMS"
echo "Imported ambient SiO2 model endpoints into $target"
