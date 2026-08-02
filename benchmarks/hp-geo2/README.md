# High-pressure GeO2 discriminator benchmark

This is the proposed superiority case, but it begins only after the canonical
PDFgui, RMCProfile, and EPSR cases pass. It must be framed as a falsifiable
comparison of structural inference, not as an assertion that EPSR "fails".

For every pressure point, compare official EPSR/Dissolve, pure RMC, pair-HRMC,
and MLIP-HRMC using identical fitted contrasts, fit ranges, density, atom count,
and starting-state policy. Use at least ten independent replicas. The primary
tests are:

- held-out isotope contrast or held-out Q/r interval;
- GeO4/GeO5/GeO6 fractions and transition pressure;
- Ge-O coordination and bond-length distributions;
- O-Ge-O and Ge-O-Ge angles, edge-sharing fraction, and ring statistics;
- MLIP energy/force plausibility and out-of-domain uncertainty;
- fit-quality-versus-physical-plausibility Pareto curves;
- end-to-end wall time including potential initialization.

An EPSR weakness can be claimed only if it is reproducible across seeds and
reasonable reference-potential/empirical-potential settings, while MLIP-HRMC
passes a held-out experimental or independent simulation test. If only one EPSR
setup is poor, the result is parameter sensitivity, not method failure.

rsmith now supports explicit per-point uncertainties, Qdamp, and per-dataset
neutron scattering-length overrides needed for this case. Experimental data,
pressure metadata, a pressure-valid MLIP, its training-domain audit, and the
official comparator inputs remain to be assembled.
