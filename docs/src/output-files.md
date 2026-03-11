# Output Files

## RMC refinement

| File | Description |
|------|-------------|
| `refined.xyz` | Refined atomic structure (extended XYZ with lattice) |
| `refined_sq.dat` | Final computed total X-ray S(Q) |
| `refined_total_gr.dat` | Final computed total X-ray g(r) via inverse FT |
| `checkpoint.dat` | Checkpoint for resuming (text format) |

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

All output files are written to the same directory as the config file.
