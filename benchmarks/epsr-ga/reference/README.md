# Local EPSR26 reference files

The EPSR26 workshop files are not committed to this repository. Their bundled
disclaimer permits use only as tutorial/example data, and publication reuse
requires separate permission from the data owner.

To create a private test copy from an existing EPSR26 installation, run:

```bash
EPSR26_ROOT=/path/to/EPSR ./import_local_reference.sh --accept-local-testing-terms
```

The command copies the `LiquidGa50C` worked example into the ignored
`reference/local/upstream/` directory and records SHA-256 checksums. Native EPSR
runs are made from another ignored copy so the installed tutorial is never
modified.
