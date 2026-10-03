# Stage 15: compare an isolated rebuild against the accepted results.
# Usage: Rscript --vanilla 03_scripts/15_rebuild_compare.R <rebuild_dir>
# Scientific agreement (numbers) and document rendering are reported separately.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
args <- commandArgs(trailingOnly = TRUE)
rb <- normalizePath(if (length(args)) args[[1]] else file.path(pm_project_root(), "_rebuilds", "clean_20260926"), mustWork = TRUE)
step_begin("15_rebuild_compare", params = list(rebuild_dir = basename(rb)))

TOL <- 1e-8
targets <- c("04_results/primary/primary_endpoint.csv", "04_results/primary/secondary_families.csv",
             "04_results/primary/GSE314812_patient_deltas.csv", "04_results/primary/adjusted_models.csv",
             "04_results/sensitivity/primary_sensitivity.csv", "04_results/sensitivity/leave_one_patient_out.csv",
             "04_results/sensitivity/single_composition_consistency.csv", "04_results/sensitivity/purity_equivalence.csv",
             "04_results/validation/GSE237876_frozen_endpoint.csv", "04_results/cell_source/host_compartment_deltas.csv",
             "04_results/audit/source_counts.csv", "04_results/audit/matrix_type_audit.csv",
             "00_admin/sample_crosswalk.csv", "00_admin/eligibility_flow.csv",
             "02_data_processed/frozen_programs/program_definitions.csv")

cmp <- rbindlist(lapply(targets, function(rel) {
  a <- file.path(pm_project_root(), rel); b <- file.path(rb, rel)
  if (!file.exists(a) || !file.exists(b))
    return(data.table(file = rel, status = "missing", detail = paste0("accepted:", file.exists(a), " rebuild:", file.exists(b))))
  x <- fread(a); y <- fread(b)
  if (!identical(dim(x), dim(y)))
    return(data.table(file = rel, status = "shape_differs", detail = sprintf("%s vs %s", paste(dim(x), collapse = "x"), paste(dim(y), collapse = "x"))))
  if (!identical(sort(names(x)), sort(names(y))))
    return(data.table(file = rel, status = "columns_differ", detail = paste(setdiff(names(x), names(y)), collapse = ",")))
  setcolorder(y, names(x))
  num <- names(x)[vapply(x, is.numeric, TRUE)]
  chr <- setdiff(names(x), num)
  maxdiff <- if (length(num)) max(abs(as.matrix(x[, ..num]) - as.matrix(y[, ..num])), na.rm = TRUE) else 0
  na_same <- all(is.na(as.matrix(x[, ..num])) == is.na(as.matrix(y[, ..num])))
  chr_same <- all(vapply(chr, function(c) identical(as.character(x[[c]]), as.character(y[[c]])), TRUE))
  ok <- (is.finite(maxdiff) && maxdiff <= TOL) && na_same && chr_same
  data.table(file = rel, status = ifelse(ok, "identical", "differs"),
             detail = sprintf("max numeric difference %.3g; NA pattern %s; text columns %s",
                              maxdiff, ifelse(na_same, "same", "differ"), ifelse(chr_same, "same", "differ")))
}))
o1 <- pm_out("_rebuilds", paste0(basename(rb), "_scientific_comparison.csv")); write_csv_atomic(cmp, o1)

figs_a <- list.files(pm_in("05_figures", "main"), pattern = "\\.pdf$")
figs_b <- list.files(file.path(rb, "05_figures", "main"), pattern = "\\.pdf$")
doc <- data.table(item = "main figures regenerated", accepted = length(figs_a), rebuild = length(figs_b),
                  note = "PDF bytes are not compared: timestamps make them differ without any scientific difference")
o2 <- pm_out("_rebuilds", paste0(basename(rb), "_document_check.csv")); write_csv_atomic(doc, o2)

rep <- c(
  sprintf("# Clean rebuild comparison (%s)", basename(rb)), "",
  sprintf("Rebuilt into an isolated output root from the same read-only raw inputs, frozen configuration and locked environment on %s.", format(Sys.Date())),
  "", sprintf("Files compared: %d. Identical: %d. Differing: %d. Missing: %d.",
              nrow(cmp), cmp[status == "identical", .N], cmp[status == "differs", .N], cmp[status == "missing", .N]),
  sprintf("Numeric tolerance: %g.", TOL), "",
  "| file | status | detail |", "|---|---|---|",
  cmp[, sprintf("| %s | %s | %s |", file, status, detail)], "",
  "## What was and was not rebuilt", "",
  "Stages rebuilt: preflight, metadata parsing, sample audit, bulk preparation, single-cell preparation, programme freeze,",
  "primary analysis, independent-source analysis, composition, measurement calibration, optional-module status, synthesis and figures.",
  "",
  "Not rebuilt in this run: manuscript assembly and reference verification, which read the result files and require network access",
  "for PubMed verification. Their inputs are the files compared above.",
  "",
  "Raw data were not re-downloaded; the rebuild reads the same read-only files whose SHA-256 checksums are recorded in",
  "00_admin/data_manifest.csv. This is a rebuild from archived processed inputs, not from raw sequencing reads.",
  "",
  sprintf("Figures regenerated in the rebuild: %d of %d.", length(figs_b), length(figs_a)),
  "PDF files are not byte-compared because embedded timestamps differ without any scientific difference."
)
o3 <- pm_out("_rebuilds", paste0(basename(rb), "_report.md")); write_text_atomic(rep, o3)
print(cmp[, .(file = substr(file, 1, 58), status)])
if (cmp[status != "identical", .N]) warning("rebuild differences present; see the comparison table")
step_end(outputs = c(o1, o2, o3))
