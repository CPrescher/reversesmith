# Output Files

By default, all output is written to the config file's directory. Use `--output-dir DIR` to redirect output to a different directory.

## RMC refinement

Starting structure (computed before refinement begins):

| File | Description |
|------|-------------|
| `start_sq.dat` | Total X-ray S(Q) of the starting structure |
| `start_gr.dat` | All partial g(r) functions of the starting structure |
| `start_total_gr.dat` | Total X-ray g(r) via inverse FT (when g(r) data is configured) |

Refined structure (after refinement completes):

| File | Description |
|------|-------------|
| `refined.xyz` | Refined atomic structure (extended XYZ with lattice) |
| `refined_sq.dat` | Final computed total X-ray S(Q) |
| `refined_gr.dat` | All partial g(r) functions of the refined structure |
| `refined_total_gr.dat` | Final computed total X-ray g(r) via inverse FT |
| `checkpoint.dat` | Periodic checkpoint; use `--resume` to continue from it |
| `reversesmith.log` | Full log of the run |

## `--compute-sq-only`

| File | Description |
|------|-------------|
| `computed_sq.dat` | Total X-ray S(Q) of the input structure |
| `computed_gr.dat` | All partial g(r) functions |
| `computed_total_gr.dat` | Total X-ray g(r) via inverse FT |

## `--analyze`

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
