# EPSR26 acquisition log

## 2026-08-02

The user explicitly accepted the official EPSR conditions for local execution
testing. No acceptance was inferred on their behalf.

The legacy official page named these packages:

- `EPSR26-Workshop-2024.zip`
- `EPSR26distribution-2022-06-06.zip`

The original link targets were recovered from an archived copy of the official
ISIS page:

- `https://www.isis.stfc.ac.uk/OtherFiles/Disordered%20Materials/EPSR26-Workshop-2024.zip`
- `https://www.isis.stfc.ac.uk/OtherFiles/Disordered%20Materials/EPSR26distribution-2022-06-06.zip`

Both live targets returned HTTP 404 after the ISIS website migration. The
legacy EPSR page itself redirected to an unrelated news item on the new site.
The current ISIS site search did not locate replacement EPSR pages.

Archive checks:

- Internet Archive contains snapshots of the EPSR HTML page through 2026, but
  its CDX index contains no successful capture of either ZIP payload.
- Common Crawl contains no capture of the workshop ZIP at the exact official
  URL checked against its 2025 crawl index.
- public web and GitHub code searches found descriptions and manuals, but no
  verifiable mirror of either official package.

No substitute archive was downloaded. The next acquisition step is to ask the
ISIS Disordered Materials Group for working official links or copies using
`PERMISSION_REQUEST.md`. If URLs are supplied, run:

```bash
EPSR_WORKSHOP_URL="<official URL>" \
EPSR_DISTRIBUTION_URL="<official URL>" \
./fetch_reference.sh --accept-local-testing-terms
```

The fetcher will write SHA-256 hashes for both received archives.

## 2026-08-03

The complete `DTBsilicaNX` worked example was located inside the user's
existing local EPSR26 installation under
`workshop2022/EPSR_WorkedExamples_2022/DTBsilicaNX`. This is an authorized
local-testing source, not a redistributable replacement archive.

`import_local_reference.sh` now copies that directory into the ignored
`reference/local/upstream` sandbox and records a complete SHA-256 manifest.
Key hashes are committed in `expected/native-forward.toml`, while all tutorial
inputs and native outputs remain ignored. The exact `.ato` hash is
`f6a5528c0f8adea3502a718e470b3b3e44b589d76556e95618686da9925882f0`.
