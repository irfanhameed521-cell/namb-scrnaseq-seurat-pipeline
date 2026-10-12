# ============================================================
# 01_run_pipeline.R
# NAMB scRNA-seq analysis (Week 3): load -> merge -> QC ->
# normalise -> variable genes -> scale -> PCA -> cluster -> UMAP
#
# Run from the repository root:
#   Rscript scripts/01_run_pipeline.R
# or open the .Rproj in RStudio and click "Source".
#
# All settings live in R/config.R.
# ============================================================

# ------------------------------------------------------------
# STEP 1: Packages, settings and helpers
# ------------------------------------------------------------
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(hdf5r)   # required by Read10X_h5()
library(here)

source(here("R", "config.R"))
source(here("R", "utils.R"))

check_input_files(cfg)
make_output_dirs(cfg)


# ------------------------------------------------------------
# STEP 2: Load the three datasets
# ------------------------------------------------------------
F_namb <- load_sample("female", cfg); F_namb
M_namb <- load_sample("male",   cfg); M_namb
Namb   <- load_sample("mixed",  cfg); Namb


# ------------------------------------------------------------
# STEP 3: Merge all three datasets
# ------------------------------------------------------------
chatall <- merge(
  F_namb,
  y = list(M_namb, Namb),
  add.cell.ids = cfg$cell_ids,
  project = "Chat_namb"
)

# Mitochondrial percentage (calculated on the merged object)
chatall[["percent.mt"]] <- PercentageFeatureSet(chatall, pattern = cfg$mt_pattern)

chatall
n_cells_before_qc <- length(WhichCells(chatall))
n_cells_before_qc


# ------------------------------------------------------------
# STEP 4: Quality control (before filtering)
# ------------------------------------------------------------
# nFeature_RNA = genes detected per cell
# nCount_RNA   = total UMIs per cell
# percent.mt   = % of reads from mitochondrial genes
p_qc_before <- VlnPlot(
  chatall,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3, pt.size = 0.001
)
save_gg(p_qc_before, "01_qc_violin_before_filtering", cfg, width = 12, height = 5)

plot1 <- FeatureScatter(chatall, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(chatall, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
save_gg(plot1 + plot2, "02_qc_scatter", cfg, width = 12, height = 5)


# ------------------------------------------------------------
# STEP 5: Filter low-quality cells
# ------------------------------------------------------------
chatall <- subset(
  chatall,
  subset = nFeature_RNA > cfg$qc_min_features &
           percent.mt   < cfg$qc_max_percent_mt
)

n_cells_after_qc <- length(WhichCells(chatall))
n_cells_after_qc

p_qc_after <- VlnPlot(
  chatall,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3, pt.size = 0.001
)
save_gg(p_qc_after, "03_qc_violin_after_filtering", cfg, width = 12, height = 5)

write.csv(
  data.frame(
    stage = c("before_qc", "after_qc"),
    n_cells = c(n_cells_before_qc, n_cells_after_qc)
  ),
  file.path(cfg$table_dir, "cell_counts_qc.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# STEP 6: Normalise
# ------------------------------------------------------------
chatall <- NormalizeData(chatall)


# ------------------------------------------------------------
# STEP 7: Highly variable genes
# ------------------------------------------------------------
chatall <- FindVariableFeatures(
  chatall,
  selection.method = "vst",
  nfeatures = cfg$n_variable_features
)

top10 <- head(VariableFeatures(chatall), 10)
top10

plot1 <- VariableFeaturePlot(chatall)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_gg(plot2, "04_variable_features", cfg, width = 8, height = 6)

write.csv(
  data.frame(top10_variable_genes = top10),
  file.path(cfg$table_dir, "top10_variable_genes.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# STEP 8: Scale
# ------------------------------------------------------------
all.genes <- rownames(chatall)
chatall <- ScaleData(chatall, features = all.genes)
length(WhichCells(chatall))


# ------------------------------------------------------------
# STEP 9: PCA
# ------------------------------------------------------------
chatall <- RunPCA(chatall, features = VariableFeatures(object = chatall))

print(chatall[["pca"]], dims = 1:5, nfeatures = 5)

save_gg(VizDimLoadings(chatall, dims = 1:2, reduction = "pca"),
        "05_pca_loadings", cfg, width = 10, height = 7)

save_gg(DimPlot(chatall, reduction = "pca"),
        "06_pca_dimplot", cfg)

save_gg(DimHeatmap(chatall, dims = 1, cells = 500, balanced = TRUE),
        "07_pca_heatmap_PC1", cfg)

save_gg(DimHeatmap(chatall, dims = 1:12, cells = 500, balanced = TRUE),
        "08_pca_heatmap_PC1-12", cfg, width = 14, height = 14)

# Elbow plot: use it to decide how many PCs to keep (cfg$n_pcs)
save_gg(ElbowPlot(chatall), "09_elbow_plot", cfg)


# ------------------------------------------------------------
# STEP 10: Choose number of PCs
# ------------------------------------------------------------
x <- cfg$n_pcs   # original analysis used 17

if (ncol(Embeddings(chatall, reduction = "pca")) < x) {
  stop("Fewer than ", x, " PCs are available. Lower 'n_pcs' in R/config.R after checking the elbow plot.")
}


# ------------------------------------------------------------
# STEP 11: Neighbours and clusters
# ------------------------------------------------------------
chatall <- FindNeighbors(chatall, dims = 1:x)
chatall <- FindClusters(chatall, resolution = cfg$cluster_resolution)

head(Idents(chatall), 5)
length(WhichCells(chatall))
length(levels(Idents(chatall)))   # number of clusters
table(Idents(chatall))            # cells per cluster


# ------------------------------------------------------------
# STEP 12: UMAP
# ------------------------------------------------------------
chatall <- RunUMAP(chatall, dims = 1:x)


# ------------------------------------------------------------
# STEP 13: Reorder clusters by relatedness + cluster tree
# ------------------------------------------------------------
# NOTE: reorder = TRUE renumbers the clusters. This happens in the
# original notebook BEFORE the UMAP plots, so we do the same here.
chatall <- BuildClusterTree(chatall, reorder = TRUE, reorder.numeric = TRUE)

save_base(function() PlotClusterTree(chatall), "10_cluster_tree", cfg)


# ------------------------------------------------------------
# STEP 14: Visualise clusters on the UMAP
# ------------------------------------------------------------
save_gg(DimPlot(chatall, reduction = "umap", label = TRUE),
        "11_umap_clusters", cfg)

save_gg(DimPlot(chatall, reduction = "umap", label = TRUE, split.by = "orig.ident"),
        "12_umap_split_by_sample", cfg, width = 16, height = 6)

save_gg(DimPlot(chatall, reduction = "umap", label = TRUE, group.by = "orig.ident"),
        "13_umap_grouped_by_sample", cfg)


# ------------------------------------------------------------
# STEP 15: Save results
# ------------------------------------------------------------
cluster_sizes <- as.data.frame(table(cluster = Idents(chatall)))
write.csv(cluster_sizes, file.path(cfg$table_dir, "cluster_sizes.csv"), row.names = FALSE)

saveRDS(chatall, file.path(cfg$obj_dir, "chatall_clustered.rds"))

writeLines(capture.output(sessionInfo()), file.path(cfg$table_dir, "session_info.txt"))

length(WhichCells(chatall))   # final cell count

message("Pipeline finished. Figures: ", cfg$fig_dir)
