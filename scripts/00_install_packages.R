# ============================================================
# 00_install_packages.R
# Run this ONCE to install everything the pipeline needs.
#   Rscript scripts/00_install_packages.R
# ============================================================

cran_packages <- c(
  "Seurat",     # single-cell analysis
  "dplyr",
  "tidyverse",
  "ggplot2",
  "cowplot",
  "patchwork",
  "hdf5r",      # needed by Read10X_h5(); requires the HDF5 system library
  "here",       # repository-relative file paths
  "rmarkdown",  # to knit the notebook
  "knitr"
)

to_install <- setdiff(cran_packages, rownames(installed.packages()))

if (length(to_install) > 0) {
  message("Installing: ", paste(to_install, collapse = ", "))
  install.packages(to_install)
} else {
  message("All required packages are already installed.")
}
