# RMCProfile MLIP CeO2-x benchmark

This reproduces the benchmark RMCProfile uses to demonstrate its own MLIP
integration, so it is the fairest MLIP-to-MLIP comparison. The published model
is a 6144-atom CeO1.71 configuration with synthetic X-ray and neutron PDFs,
ReaxFF-derived vacancy ordering, a FLARE Gaussian-process potential, and oxygen
vacancy/site-swap moves. The key observable is the evolution of <100> vacancy
pairs with and without the potential constraint.

Primary reference: P. Cuillier, M. G. Tucker and Y. Zhang, *J. Appl. Cryst.*
**57** (2024), DOI `10.1107/S1600576724009282`.

rsmith cannot yet execute this case because it has displacement moves but not
vacancy/site swaps. Implementation acceptance requires exact composition
preservation, incremental histogram correctness after swaps, energy-backend
accept/reject transactions, and detailed-balance tests. Upstream configuration,
synthetic PDFs, FLARE model, and vacancy analysis code must also be acquired
from the authors or supplementary archive before comparison.
