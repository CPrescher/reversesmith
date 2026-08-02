# Official EPSR amorphous-silica benchmark

Amorphous silica is one of the examples selected by the EPSR authors in the
EPSR manual. The comparison must use the official EPSR26 implementation and
its reference-potential plus empirical-potential workflow. rsmith's internal
experimental EPSR mode is not an acceptable stand-in for the comparator.

The official download page states that bundled raw data and demonstration files
may not be used in publicity or publication without express permission. They
may be used locally for execution testing, but they are not vendored here and
must not enter a paper archive until STFC/ISIS grants permission. The official
page, rather than this repository, is authoritative for the current terms.

On 2026-08-02 the user accepted those conditions for local execution testing.
The historical official ZIP targets were then resolved from the archived ISIS
page, but both returned HTTP 404 from the migrated live site. Neither Internet
Archive nor Common Crawl held the ZIP payload. This is an upstream acquisition
blocker, not permission to replace the files with an unverified mirror. See
`ACQUISITION_LOG.md`.

The EPSR26 manual's silica exercise fixes several important starting points:

- number density: 0.06834 atoms/A^3;
- equilibrate the initially random network at 10,000 K before cooling to 300 K;
- neutron normalization: sum of neutron weights (`nrtype = 3`);
- X-ray normalization: single-atom scattering;
- use both neutron and X-ray total-scattering data.

Required sequence:

1. obtain written publication/reuse permission using `PERMISSION_REQUEST.md`;
2. obtain working official package URLs from ISIS, then acquire the files with
   `fetch_reference.sh --accept-local-testing-terms` (URL overrides are
   supported);
3. record the archive hashes, EPSR26 build, reference potentials, density, atom
   count, equilibration length, refinement length, and analysis commands;
4. reproduce native EPSR with ten independent seeds;
5. reproduce that ensemble with rsmith's EPSR mode before making any
   superiority comparison;
6. repeat with rsmith pure RMC, pair-potential HRMC, and MLIP-HRMC at matched
   fit quality;
7. compare fitted totals, partial RDFs, coordination, O-Si-O and Si-O-Si
   angles, ring topology, wall time, and seed-to-seed variability.

The comparison design is preregistered in `protocol.toml`. Numeric stochastic
equivalence limits will be frozen only after the native EPSR replicate spread
is known and before inspecting any rsmith result. This avoids defining
"equivalent" from the answer we want.

## Current implementation blockers

The canonical example is deliberately useful because it exposes two real gaps
rather than allowing an easy but misleading comparison:

1. rsmith can fit multiple datasets during RMC, but its empirical-potential
   update currently uses only one total-scattering dataset. Native EPSR silica
   uses neutron and X-ray data together.
2. the exact EPSR reference potential and mixing conventions must be imported
   faithfully. We will choose between a native implementation and tabulated
   pair potentials only after inspecting the permitted upstream input files.

These gaps must be closed before `rsmith_epsr` can be called a reproduction of
the native EPSR calculation.

The scientific claim is not "smaller chi-squared than EPSR". The meaningful
test is whether MLIP-HRMC gives better held-out or chemically diagnostic
structure at the same data agreement.

## Files

- `protocol.toml`: frozen methods, stages, contrasts, seeds, and observables;
- `PERMISSION_REQUEST.md`: ready-to-send publication/reuse request;
- `fetch_reference.sh`: explicit, locally gated acquisition with SHA-256 log;
- `reference/README.md`: provenance and redistribution rules for upstream data.
