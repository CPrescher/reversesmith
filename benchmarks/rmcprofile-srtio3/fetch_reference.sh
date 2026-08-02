#!/usr/bin/env bash
set -euo pipefail

case_dir="$(cd "$(dirname "$0")" && pwd)"
archive="$case_dir/reference/RMCProfile_V7b35_Linux_ARM64.tgz"
url="https://sourceforge.net/projects/rmcprofile/files/Version_7/RMCProfile_V7b35_Linux_ARM64.tgz/download"
expected="828b9292167a6691e02d1a9073006a42be7300ce770309702a827e5886b105c3"

curl -L --fail --show-error -o "$archive" "$url"
actual="$(shasum -a 256 "$archive" | awk '{print $1}')"
if [[ "$actual" != "$expected" ]]; then
    echo "checksum mismatch: expected $expected, got $actual" >&2
    exit 1
fi
echo "downloaded verified archive to $archive"
