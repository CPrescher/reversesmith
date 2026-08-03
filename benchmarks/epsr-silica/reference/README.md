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
