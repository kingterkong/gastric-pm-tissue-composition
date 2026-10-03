# Staged dependency installation into the project-local renv library.
# Usage: Rscript --vanilla 03_scripts/R/install_deps.R <stage>
# Stages: base, bioc_bulk, singlecell, docs. The old project library is never modified.
args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) args[[1]] else "base"
proj <- normalizePath(".")
stopifnot(basename(proj) == "project_pm")
source(file.path(proj, "renv", "activate.R"))
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
pk <- list(
  base = c("data.table", "ggplot2", "patchwork", "cowplot", "ggrepel", "viridisLite", "scales",
           "dplyr", "tidyr", "readr", "stringr", "purrr", "tibble", "jsonlite", "yaml", "digest",
           "openxlsx", "writexl", "readxl", "Matrix", "boot", "httr", "xml2", "rentrez", "curl",
           "metafor", "R.utils", "BiocManager", "ragg", "svglite", "systemfonts", "knitr", "rmarkdown"),
  bioc_bulk = c("bioc::limma", "bioc::edgeR", "bioc::GEOquery", "bioc::Biobase",
                "bioc::AnnotationDbi", "bioc::org.Hs.eg.db", "bioc::singscore", "bioc::GSEABase",
                "bioc::ComplexHeatmap"),
  singlecell = c("bioc::SingleCellExperiment", "bioc::scuttle", "bioc::UCell", "bioc::scDblFinder",
                 "SeuratObject", "Seurat", "hdf5r", "irlba", "uwot"),
  docs = c("officer", "flextable", "pdftools", "magick", "qpdf")
)
todo <- pk[[stage]]
stopifnot(!is.null(todo))
message("Installing stage ", stage, ": ", paste(todo, collapse = ", "))
renv::install(todo, prompt = FALSE)
ip <- rownames(installed.packages(lib.loc = .libPaths()[1]))
want <- sub("^bioc::", "", todo)
message("Missing after install: ", paste(setdiff(want, ip), collapse = ", "))
