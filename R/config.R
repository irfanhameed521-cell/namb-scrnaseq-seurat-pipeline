# ============================================================
# config.R - all editable settings for the pipeline
# ============================================================
# Every path is relative to the repository root (via the `here`
# package), so nothing here needs to be changed when you move the
# project to another computer. Only edit the values you want to tune.

library(here)

cfg <- list(

  # ---- Input data ----------------------------------------
  # Put the three 10X .h5 files inside data/raw/
  data_dir = here("data", "raw"),
  files = c(
    female = "female_filtered_feature_bc_matrix.h5",
    male   = "male_filtered_feature_bc_matrix.h5",
    mixed  = "mixed_filtered_feature_bc_matrix.h5"
  ),

  # Seurat project names given to each sample (becomes `orig.ident`)
  project_names = c(
    female = "Namb_female",
    male   = "Namb_male",
    mixed  = "Namb_mixed"
  ),

  # Prefixes added to cell barcodes when merging
  cell_ids = c("Female", "Male", "Mixed"),

  # ---- Object creation (per sample) ------------------------
  min_cells    = 3,
  min_features = 200,

  # ---- Quality control -----------------------------------
  # "^mt-" = mouse-style gene names (mt-Nd1, mt-Co1 ...)
  # Use "^MT-" for human-style names (MT-ND1, MT-CO1 ...)
  mt_pattern        = "^mt-",
  qc_min_features   = 100,   # keep cells with nFeature_RNA > this
  qc_max_percent_mt = 2,     # keep cells with percent.mt < this

  # ---- Feature selection / dimensionality reduction ---------
  n_variable_features = 2000,
  n_pcs               = 17,   # chosen from the ElbowPlot

  # ---- Clustering ----------------------------------------
  cluster_resolution = 1,

  # ---- Output folders ------------------------------------
  fig_dir   = here("results", "figures"),
  table_dir = here("results", "tables"),
  obj_dir   = here("results", "objects")
)
