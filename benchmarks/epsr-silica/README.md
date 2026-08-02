# Official EPSR amorphous-silica benchmark

Amorphous silica is one of the examples selected by the EPSR authors in the
EPSR manual. The comparison must use the official EPSR26 implementation and
its reference-potential plus empirical-potential workflow. rsmith's internal
experimental EPSR mode is not an acceptable stand-in for the comparator.

The official download page states that bundled raw data and demonstration files
may not be used in publicity or publication without express permission. They
may be used locally for execution testing, but they are not vendored here and
must not enter a paper archive until STFC/ISIS grants permission.

Required sequence:

1. obtain written publication/reuse permission for the EPSR26 workshop silica
   files;
2. record EPSR26 build, reference potentials, density, atom count, equilibration
   length, refinement length, and analysis commands;
3. convert the final EPSR ensemble without changing coordinates or labels;
4. compare fitted totals, partial RDFs, Si coordination, O-Si-O and Si-O-Si
   angles, energy, wall time, and seed-to-seed variability;
5. repeat with rsmith pure RMC, pair-potential HRMC, and MLIP-HRMC at matched
   fit quality.

The scientific claim is not "smaller chi-squared than EPSR". The meaningful
test is whether MLIP-HRMC gives better held-out or chemically diagnostic
structure at the same data agreement.
