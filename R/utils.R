# ============================================================
# utils.R - small helper functions used by the pipeline
# ============================================================

# Stop early with a clear message if any input file is missing
check_input_files <- function(cfg) {
  paths <- file.path(cfg$data_dir, cfg$files)
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(
      "The following input file(s) were not found:\n  ",
      paste(missing, collapse = "\n  "),
      "\n\nPlace your 10X .h5 files in 'data/raw/' (see data/README.md), ",
      "or edit the file names in R/config.R.",
      call. = FALSE
    )
  }
  invisible(paths)
}

# Create output folders if they do not exist
make_output_dirs <- function(cfg) {
  for (d in c(cfg$fig_dir, cfg$table_dir, cfg$obj_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

# Load one 10X .h5 file and create a Seurat object
load_sample <- function(sample, cfg) {
  counts <- Seurat::Read10X_h5(
    filename        = file.path(cfg$data_dir, cfg$files[[sample]]),
    use.names       = TRUE,
    unique.features = TRUE
  )
  Seurat::CreateSeuratObject(
    counts       = counts,
    project      = cfg$project_names[[sample]],
    min.cells    = cfg$min_cells,
    min.features = cfg$min_features
  )
}

# Save a ggplot / patchwork figure to results/figures/ and also show it
save_gg <- function(plot, name, cfg, width = 8, height = 6) {
  ggplot2::ggsave(
    filename = file.path(cfg$fig_dir, paste0(name, ".png")),
    plot = plot, width = width, height = height, dpi = 300
  )
  print(plot)
  invisible(plot)
}

# Save a base-R style plot (e.g. PlotClusterTree) to results/figures/
# `plot_fun` is a function with no arguments that draws the plot.
save_base <- function(plot_fun, name, cfg, width = 8, height = 6) {
  grDevices::png(
    filename = file.path(cfg$fig_dir, paste0(name, ".png")),
    width = width, height = height, units = "in", res = 300
  )
  plot_fun()
  grDevices::dev.off()
  plot_fun()   # also draw on screen
  invisible(NULL)
}
