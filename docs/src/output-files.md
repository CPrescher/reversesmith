# Output Files

By default, all output is written to the config file's directory. Use `--output-dir DIR` to redirect output to a different directory.

## RMC refinement

Starting structure (computed before refinement begins):

| File | Description |
|------|-------------|
| `start_xray_sq.dat` | Total X-ray S(Q) of the starting structure |
| `start_neutron_sq.dat` | Total neutron S(Q) (only when `[data.neutron_sq]` is configured) |
| `start_gr.dat` | All partial g(r) functions of the starting structure |
| `start_total_gr.dat` | Total X-ray g(r) via inverse FT (when g(r) data is configured) |

Refined structure (after refinement completes):

| File | Description |
|------|-------------|
| `refined.xyz` | Refined atomic structure (extended XYZ with lattice) |
| `refined_xray_sq.dat` | Final computed total X-ray S(Q) |
| `refined_neutron_sq.dat` | Final computed total neutron S(Q) (only when `[data.neutron_sq]` is configured) |
| `refined_gr.dat` | All partial g(r) functions of the refined structure |
| `refined_total_gr.dat` | Final computed total X-ray g(r) via inverse FT |
| `checkpoint.dat` | Periodic checkpoint; use `--resume` to continue from it |
| `rsmith.log` | Full log of the run |

## `--analyze`

Analysis uses a separate log file (`rsmith_analyze.log`) to avoid overwriting the RMC refinement log.

Without explicit path (dual mode):

| File | Description |
|------|-------------|
| `analysis_starting_cn.dat` | CN histograms for starting structure |
| `analysis_starting_angles.dat` | Angle distributions for starting structure |
| `analysis_refined_cn.dat` | CN histograms for refined structure |
| `analysis_refined_angles.dat` | Angle distributions for refined structure |

With explicit path:

| File | Description |
|------|-------------|
| `analysis_analysis_cn.dat` | CN histograms |
| `analysis_analysis_angles.dat` | Angle distributions |
| `rsmith_analyze.log` | Log of the analysis run |
