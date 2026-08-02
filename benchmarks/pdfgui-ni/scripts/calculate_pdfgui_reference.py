#!/usr/bin/env python3
"""Calculate the pinned Ni curve with PDFgui's PDFfit2 engine."""

from pathlib import Path

try:
    from diffpy.pdffit2 import PdfFit
except ImportError as error:
    raise SystemExit(
        "diffpy.pdffit2 is required to generate the independent reference curve"
    ) from error

structure = Path("reference/Ni.stru")
output = Path("reference/pdfgui-calculated.dat")
calculator = PdfFit(create_intro=False)
calculator.read_struct(str(structure))
calculator.alloc("X", 27.0, 0.08, 0.01, 20.0, 2000)
calculator.calc()
output.write_text(
    "# PDFfit2/PDFgui calculated X-ray G(r): Qmax=27 1/A, Qdamp=0.08 1/A\n"
    + "".join(
        f"{r:.8f} {value:.12g}\n"
        for r, value in zip(calculator.getR(), calculator.getpdf_fit())
    )
)
print(f"wrote {len(calculator.getR())} points to {output}")
