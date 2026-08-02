#!/usr/bin/env bash
set -euo pipefail

download_page="https://www.isis.stfc.ac.uk/Pages/Empirical-Potential-Structure-Refinement.aspx"
historical_workshop_url="https://www.isis.stfc.ac.uk/OtherFiles/Disordered%20Materials/EPSR26-Workshop-2024.zip"
historical_distribution_url="https://www.isis.stfc.ac.uk/OtherFiles/Disordered%20Materials/EPSR26distribution-2022-06-06.zip"
workshop_url="${EPSR_WORKSHOP_URL:-${historical_workshop_url}}"
distribution_url="${EPSR_DISTRIBUTION_URL:-${historical_distribution_url}}"

if [[ "${1:-}" != "--accept-local-testing-terms" ]]; then
  echo "Nothing was downloaded."
  echo
  echo "Read the current official conditions first:"
  echo "  ${download_page}"
  echo
  echo "The conditions include non-commercial use and restrict bundled raw/demo"
  echo "files to execution testing unless express publication permission is granted."
  echo "Only you can accept those conditions. If you accept them for this local"
  echo "testing download, rerun:"
  echo "  $0 --accept-local-testing-terms"
  exit 2
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
download_dir="${script_dir}/reference/downloads"
mkdir -p "${download_dir}"

download_archive() {
  local output_name="$1"
  local source_url="$2"
  local output_path="${download_dir}/${output_name}"
  if ! curl --fail --location --output "${output_path}" "${source_url}"; then
    rm -f "${output_path}"
    echo >&2
    echo "Could not download ${output_name} from:" >&2
    echo "  ${source_url}" >&2
    echo "The legacy ISIS links returned HTTP 404 after the website migration." >&2
    echo "Ask the ISIS Disordered Materials Group for the current official URLs," >&2
    echo "then set EPSR_WORKSHOP_URL and EPSR_DISTRIBUTION_URL when rerunning." >&2
    return 1
  fi
}

download_archive "EPSR26-Workshop-2024.zip" "${workshop_url}"
download_archive "EPSR26distribution-2022-06-06.zip" "${distribution_url}"

if command -v sha256sum >/dev/null 2>&1; then
  (
    cd "${download_dir}"
    sha256sum EPSR26-Workshop-2024.zip EPSR26distribution-2022-06-06.zip > SHA256SUMS
  )
else
  (
    cd "${download_dir}"
    shasum -a 256 EPSR26-Workshop-2024.zip EPSR26distribution-2022-06-06.zip > SHA256SUMS
  )
fi

echo "Downloaded for local execution testing and recorded SHA-256 hashes in:"
echo "  ${download_dir}/SHA256SUMS"
echo "Do not redistribute or use the raw/demo files in a publication without"
echo "express permission from STFC/ISIS."
