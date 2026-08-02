#!/usr/bin/env python3
"""Extract r, G(r), sigma from the PDFgetX2 four-column tutorial file."""

from pathlib import Path

source = Path("reference/ni-q27r60-xray.gr")
target = Path("reference/ni-xray-r-G-sigma.dat")
rows = []
for raw in source.read_text().splitlines():
    fields = raw.split()
    if len(fields) < 4:
        continue
    try:
        values = [float(value) for value in fields[:4]]
    except ValueError:
        continue
    rows.append((values[0], values[1], values[3]))
if not rows:
    raise SystemExit(f"no four-column numeric rows found in {source}")
target.write_text(
    "# r(A) G(r)(1/A^2) sigma\n"
    + "".join(f"{r:.8g} {g:.12g} {sigma:.12g}\n" for r, g, sigma in rows)
)
print(f"wrote {len(rows)} points to {target}")
