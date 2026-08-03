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
