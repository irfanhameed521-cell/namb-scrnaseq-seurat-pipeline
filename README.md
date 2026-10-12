# NAMB scRNA-seq Analysis Pipeline (Seurat)

A reproducible single-cell RNA-seq workflow in R (Seurat) that loads three 10X Genomics datasets (**female**, **male**, **mixed**), merges them, performs quality control, and produces clusters and UMAP visualisations.


## What the pipeline does

1. Loads three 10X `.h5` files and creates a Seurat object for each
2. Merges them into one object (`chatall`)
3. Calculates mitochondrial percentage and runs QC plots
4. Filters low-quality cells (`nFeature_RNA > 100` and `percent.mt < 2`)
5. Normalises, finds 2000 highly variable genes, and scales the data
6. Runs PCA and shows the elbow plot to choose the number of PCs (17 used)
7. Builds the neighbour graph, finds clusters (resolution 1), runs UMAP
8. Reorders clusters by relatedness (cluster tree) and plots UMAPs by cluster and by sample

See [`docs/pipeline_overview.md`](docs/pipeline_overview.md) for a step-by-step table.

## Repository structure

```
namb-scrnaseq-seurat-pipeline/
├── README.md
├── LICENSE
├── .gitignore
├── namb-scrnaseq-seurat-pipeline.Rproj   # open this in RStudio
├── R/
│   ├── config.R                # all editable settings (paths, thresholds, PCs, resolution)
│   └── utils.R                 # helper functions (load samples, save plots, file checks)
├── scripts/
│   ├── 00_install_packages.R   # one-time package installation
│   └── 01_run_pipeline.R       # full analysis as a script
├── notebooks/
│   └── analysis.Rmd      # same analysis as an R Markdown notebook
├── data/
│   ├── README.md               # which input files to provide
│   └── raw/                    # put your .h5 files here (not tracked by git)
├── results/
│   ├── figures/                # PNG plots (generated)
│   ├── tables/                 # CSV summaries + session info (generated)
│   └── objects/                # saved Seurat object, .rds (generated, not tracked)
└── docs/
    └── pipeline_overview.md
```

## Requirements

- R >= 4.1
- R packages: Seurat, dplyr, tidyverse, ggplot2, cowplot, patchwork, hdf5r, here, rmarkdown, knitr
- `hdf5r` needs the HDF5 system library (Windows and macOS binaries usually include it; on Linux install `libhdf5-dev`)

Install everything with:

```r
source("scripts/00_install_packages.R")
```

## Input data

Raw data is not included in the repository. Put these files in `data/raw/`:

| File | Sample |
|------|--------|
| `female_filtered_feature_bc_matrix.h5` | Female |
| `male_filtered_feature_bc_matrix.h5` | Male |
| `mixed_filtered_feature_bc_matrix.h5` | Mixed |

If your file names differ, change them in `R/config.R`.

## How to run

**Option A: script (recommended)**

```bash
git clone https://github.com/<your-username>/namb-scrnaseq-seurat-pipeline.git
cd namb-scrnaseq-seurat-pipeline
Rscript scripts/00_install_packages.R   # first time only
Rscript scripts/01_run_pipeline.R
```

**Option B: RStudio**

1. Open `namb-scrnaseq-seurat-pipeline.Rproj`
2. Run `scripts/00_install_packages.R` once
3. Open `scripts/01_run_pipeline.R` and click **Source**, or open `notebooks/week3_analysis.Rmd` and click **Knit**

No `setwd()` is needed: paths are resolved relative to the project root using the `here` package.

## Configuration

Edit `R/config.R` to change parameters:

| Setting | Default | Meaning |
|---------|---------|---------|
| `min_cells` / `min_features` | 3 / 200 | Filters used when creating each Seurat object |
| `mt_pattern` | `"^mt-"` | Mitochondrial gene pattern (use `"^MT-"` for human genes) |
| `qc_min_features` | 100 | Keep cells with more than this many genes |
| `qc_max_percent_mt` | 2 | Keep cells with less than this % mitochondrial reads |
| `n_variable_features` | 2000 | Number of highly variable genes |
| `n_pcs` | 17 | Number of principal components used |
| `cluster_resolution` | 1 | Clustering resolution |

## Outputs

After a run you will find:

- `results/figures/`: QC plots, variable features, PCA plots, elbow plot, cluster tree, UMAPs
- `results/tables/`: `cell_counts_qc.csv`, `top10_variable_genes.csv`, `cluster_sizes.csv`, `session_info.txt`
- `results/objects/chatall_clustered.rds`: the final Seurat object (not tracked by git)

## Notes

- The mixed-sample file name differs between the original notebook (`filtered_feature_bc_matrix.h5`) and the cleaned script (`mixed_filtered_feature_bc_matrix.h5`). This repo uses the cleaned script's name.
- `BuildClusterTree(reorder = TRUE)` renumbers clusters, so final cluster numbers differ from the raw `FindClusters` output.
- Check the elbow plot before reusing `n_pcs = 17` on different data.

## License

MIT, see [LICENSE](LICENSE). Replace `<YOUR NAME>` in the license file with your name.
