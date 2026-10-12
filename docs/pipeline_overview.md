# Pipeline overview

| # | Step | Seurat function(s) | Key parameter (in `R/config.R`) |
|---|------|--------------------|-------------------------------|
| 1 | Load each 10X `.h5` file | `Read10X_h5`, `CreateSeuratObject` | `min_cells = 3`, `min_features = 200` |
| 2 | Merge female, male, mixed | `merge` | `cell_ids` |
| 3 | Mitochondrial % | `PercentageFeatureSet` | `mt_pattern = "^mt-"` |
| 4 | QC plots | `VlnPlot`, `FeatureScatter` | |
| 5 | Filter cells | `subset` | `qc_min_features = 100`, `qc_max_percent_mt = 2` |
| 6 | Normalise | `NormalizeData` | |
| 7 | Highly variable genes | `FindVariableFeatures` (vst) | `n_variable_features = 2000` |
| 8 | Scale | `ScaleData` (all genes) | |
| 9 | PCA | `RunPCA`, `ElbowPlot`, `DimHeatmap` | |
| 10 | Neighbours + clusters | `FindNeighbors`, `FindClusters` | `n_pcs = 17`, `cluster_resolution = 1` |
| 11 | UMAP | `RunUMAP` | `n_pcs` |
| 12 | Cluster tree (reorders cluster numbers) | `BuildClusterTree` | |
| 13 | UMAP plots by cluster / sample | `DimPlot` | |

## Notes

- **Mouse vs human genes:** the default mitochondrial pattern `^mt-` matches mouse-style gene names. For human data change it to `^MT-`.
- **Choosing PCs:** 17 PCs were used in the original analysis. Check `results/figures/09_elbow_plot.png` before reusing the value on new data.
- **Cluster renumbering:** `BuildClusterTree(reorder = TRUE, reorder.numeric = TRUE)` renumbers clusters by relatedness, so cluster numbers differ from the raw `FindClusters` output.
- **Reproducibility:** `results/tables/session_info.txt` records the exact package versions used in each run.
