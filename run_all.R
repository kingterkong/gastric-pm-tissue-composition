#!/usr/bin/env Rscript
# project_pm main entry point.
# Usage (from project_pm/):
#   Rscript --vanilla run_all.R --stage preflight|audit|prepare|freeze|analyze|figures|manuscript|verify|posthoc|all
#                               [--resume] [--output-root <dir under project_pm or project_pm/_rebuilds>]
# Each stage script runs in a fresh R process; --resume skips a script only when its
# success marker is present and all recorded input/output/script/config/protocol hashes match.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value for ", flag)
  args[[i + 1]]
}
stage <- get_arg("--stage", "all")
resume <- "--resume" %in% args
script_path <- tryCatch(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]), error = function(e) NA)
proj <- if (!is.na(script_path)) dirname(script_path) else normalizePath(".")
out_root <- get_arg("--output-root", proj)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
out_root <- normalizePath(out_root)
Sys.setenv(PM_PROJECT_ROOT = proj, PM_OUTPUT_ROOT = out_root, RENV_PROJECT = proj)

source(file.path(proj, "renv", "activate.R"))
source(file.path(proj, "03_scripts", "R", "utils.R"))
check_output_root()

stages <- list(
  preflight  = c("00_preflight"),
  audit      = c("01_literature_and_metadata", "02_audit_samples"),
  prepare    = c("03_download_inputs", "04_prepare_bulk", "05_prepare_singlecell"),
  freeze     = c("06_define_and_freeze_programs"),
  analyze    = c("07_primary_paired_analysis", "08_independent_validation", "09_composition_and_coupling",
                 "09b_measurement_model", "10_optional_spatial_clinical", "11_sensitivity_and_synthesis"),
  figures    = c("12_generate_figures_tables"),
  manuscript = c("13a_references", "13_build_manuscript"),
  verify     = c("14_verify_delivery"),
  posthoc    = c("17_posthoc_robustness", "18_generate_fig_submission_assets")
)
if (!stage %in% c(names(stages), "all")) stop("Unknown --stage: ", stage)
todo <- if (stage == "all") unlist(stages, use.names = FALSE) else stages[[stage]]

rscript <- file.path(R.home("bin"), "Rscript")
for (s in todo) {
  f <- file.path(proj, "03_scripts", paste0(s, ".R"))
  if (!file.exists(f)) stop("Stage script not found: ", f)
  if (resume && step_is_current(s)) { cat(sprintf("[resume] %s is current; skipped\n", s)); next }
  invalidate_step(s)
  cat(sprintf("==> %s\n", s))
  st <- system2(rscript, c("--vanilla", shQuote(f)), env = c(sprintf("PM_PROJECT_ROOT=%s", shQuote(proj)),
                                                              sprintf("PM_OUTPUT_ROOT=%s", shQuote(out_root)),
                                                              sprintf("RENV_PROJECT=%s", shQuote(proj))))
  if (!identical(st, 0L)) {
    # invalidate everything downstream of the failed step
    down <- todo[seq(match(s, todo), length(todo))]
    for (d in down) invalidate_step(d)
    stop(sprintf("Stage script %s failed with exit status %s; downstream markers invalidated", s, st))
  }
}
cat("run_all.R finished: ", paste(todo, collapse = ", "), "\n")
