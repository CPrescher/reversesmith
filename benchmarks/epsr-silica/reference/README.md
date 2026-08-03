# EPSR26 upstream files

No upstream archive is committed here.

The authoritative source and current conditions are on the ISIS
[EPSR download page](https://www.isis.stfc.ac.uk/Pages/Empirical-Potential-Structure-Refinement.aspx).
At the time this benchmark was designed, the relevant tutorial archive was
`EPSR26-Workshop-2024.zip` and the implementation distribution was
`EPSR26distribution-2022-06-06.zip`.

Those exact links returned HTTP 404 on 2026-08-02 after redirecting to the
migrated ISIS site. The former page and link labels remain visible in search
indexes and archived HTML, but the ZIP payloads were not preserved by Internet
Archive or Common Crawl. See `../ACQUISITION_LOG.md` for the audit trail.

The workshop raw and demonstration files may be used locally for execution
testing under the official download conditions, but publication use requires
express permission. Do not copy them into the repository, a paper supplement,
or a public data archive until that permission explicitly permits it.

After an authorized local download, `reference/downloads/SHA256SUMS` records
the exact archives used. An existing authorized EPSR26 installation can instead
be imported with `../import_local_reference.sh --accept-local-testing-terms`;
its files and complete hash manifest remain below ignored `reference/local/`.
If redistribution remains prohibited, the committed key checksums, official
link, and conversion scripts form the public provenance trail.

## Ambient GAP and Pedone model endpoints

The user's private `CPrescher/SiO2_glass` repository contains two additional
3,000-atom, 300 K endpoints: a Pedone-pre-equilibrated GAP quench and a
Pedone-only fast quench. The benchmark pins commit
`b05590846c87eb58cf1ed2e09a6787c1d67f9e53` and the two structure hashes in
`../expected/ambient-model-endpoints.toml`.

That private repository currently has no repository-level licence. Therefore
the structures are not vendored here. Use `../import_ambient_models.sh` with an
authenticated local checkout; the coordinates remain under ignored
`reference/local/ambient-models/`. Only provenance, derived diagnostics, and
reproduction tooling are committed.

The endpoints are secondary model references, not experimental truth. Each is
one configuration from the same deterministic seed and a fast `10^14 K/s`
quench. Their equilibrium densities also differ: 2.282 g/cm3 for GAP and
2.560 g/cm3 for Pedone, compared with 2.2 g/cm3 nominal experimental glass.
Coordination, angle, and network-topology comparisons are useful; raw RDF
distance must be interpreted with the density difference stated explicitly.

## Public silica GAP parameterization

The private structure repository contains a GAP-generated endpoint but not the
potential files. The benchmark therefore obtains the Erhard et al. silica GAP
from the primary Zenodo record, DOI `10.5281/zenodo.6353684`. The deposited
archive `sio2_potential_data.zip` is 157,430,123 bytes with MD5
`6a16de69b5e17fd18160d9b55972a3e1`; the extracted `silica_gap.xml` has SHA-256
`5d470c3cde09a26c7919caa556c84902a04a35e54fd87f27bd4abb2921150799`.

Run `../scripts/fetch_public_gap.py` to download, verify, and extract the XML
and its five sparse companions below ignored `reference/local/public-gap/`.
The files are not vendored because the two SOAP sparse tables alone occupy
about 334 MB. The model label used by QUIP is
`GAP_2021_4_19_120_7_32_55_336`.

## Public 2024 Si-O ACE parameterization

The native PACE comparison uses the Si-O ACE model from Erhard, Rohrer, Albe,
and Deringer, *Nature Communications* **15**, 1927 (2024), article DOI
`10.1038/s41467-024-45840-9` and dataset DOI
`10.5281/zenodo.10419194`. Only the public `potential.zip` archive is needed;
the much larger deposited trajectory archive is not downloaded.

Run `../scripts/fetch_public_ace2024.py` to download and checksum-gate the
archive. The ZIP has MD5 `c8eb4cd111af60d10c15bc3f1de9adbb` and SHA-256
`28b29becd2c3185c6a44e872f304af7689b30b22842f21ec91e52f3641dd72cb`.
The extracted `SiOx_potential.yace` has SHA-256
`c8f00d8f0cbc131b0298b79260ba8098975624363c9d178223e51e48f025e97a`.
It remains under ignored `reference/local/public-ace2024/`; provenance,
independent oracle values, and reproduction tooling are committed.

The model uses the PACE `SBessel` radial basis. Native rsmith energies for the
3,000-atom cross-start and for a nonzero single-atom displacement agree with
the official `python-ace` implementation to the six decimals printed by the
CLI. This checks model parsing and energy evaluation, not force correctness or
whether a configuration lies inside the model's training domain.
