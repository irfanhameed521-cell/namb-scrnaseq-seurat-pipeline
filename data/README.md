# Data

Raw sequencing data is **not stored in this repository** (files are large and are excluded by `.gitignore`).

Place the following 10X Genomics `filtered_feature_bc_matrix` HDF5 files in `data/raw/`:

| File name | Sample |
|-----------|--------|
| `female_filtered_feature_bc_matrix.h5` | Female (batch 1) |
| `male_filtered_feature_bc_matrix.h5`   | Male (batch 2) |
| `mixed_filtered_feature_bc_matrix.h5`  | Mixed (batch 3) |

If your files have different names, edit the `files` entry in `R/config.R`.
(The original R Markdown notebook used `filtered_feature_bc_matrix.h5` for the mixed sample; the cleaned R script used `mixed_filtered_feature_bc_matrix.h5`. This repo follows the cleaned script.)
