#!/usr/bin/env sh
set -eu

mkdir -p reference
base="https://raw.githubusercontent.com/diffpy/cmi_exchange/b76dda7dd899b0a391009225d11cf5651295e887/cmi_scripts/fitNiPDF"
curl --fail --location --output reference/ni.cif "$base/ni.cif"
curl --fail --location --output reference/ni-q27r60-xray.gr "$base/ni-q27r60-xray.gr"
curl --fail --location --output reference/Ni.stru \
  "https://raw.githubusercontent.com/diffpy/diffpy.pdfgui/9496a699c130efaa9ce620c29e7fe95b88fcad6f/tests/testdata/Ni.stru"

printf '%s  %s\n' \
  ee4d128747da0e0ab8cb0deca8170d9f39d612f5756a08007348d1d96f31b1b3 reference/ni.cif \
  87c484a36e2a6f987381b586fb7e124139830a897c2b7d4c7c7857c7128e93f2 reference/ni-q27r60-xray.gr \
  fd5a6ef34f723e6d06c6edccd4fd48344972de87f08adcccb97ad6a0d9ca29ab reference/Ni.stru \
  | shasum -a 256 --check
