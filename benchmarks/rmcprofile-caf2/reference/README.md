# Native RMCProfile Figure 11 artifacts

The original Figure 11 files are not attached to the open-access article and
are not included in the official RMCProfile7 distribution.

Sources checked on 2026-08-02:

- Article DOI: <https://doi.org/10.1107/S1600576724004175>
- Official RMCProfile7 downloads:
  <https://sourceforge.net/projects/rmcprofile/files/Version_7/>
- `RMCProfile_V7b35_Linux_ARM64.tgz`, SHA-256
  `828b9292167a6691e02d1a9073006a42be7300ce770309702a827e5886b105c3`

The archive contains SrTiO3, SF6, and GaPO4 tutorial cases, but no path or file
name containing `CaF2`, `calcium`, or `Figure 11`. The publication names
Wojciech A. Slawinski as the main developer and directs enquiries to
`wslawinski@chem.uw.edu.pl`.

Request the following publication artifacts from the corresponding author:

1. the exact CaF2 `.rmc6f`/`.rmc7` configuration used in Figure 11;
2. the RMCProfile7 control file and exact executable revision;
3. histogram PDF output;
4. back-Fourier PDF outputs for every plotted Qmax;
5. lattice, ADP, scattering-length, RDF-grid, and resolution settings;
6. permission to redistribute the files in a reproducibility archive.

Until those files are supplied, `../scripts/run_qmax_series.py` remains a
transparent reconstruction of the published protocol, not a native-output
cross-code comparison.
