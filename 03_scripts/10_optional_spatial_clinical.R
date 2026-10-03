# Stage 10: optional modules (spatial single case; clinical peritoneal recurrence).
# Both are OFF by protocol default. This script records why, so a closed module leaves an
# explicit reason rather than silently missing output (plan section 9).
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
cfg <- pm_config()
step_begin("10_optional_spatial_clinical", inputs = pm_in("00_admin", "protocol_lock.json"), uses_protocol = TRUE)

modules <- data.table(
  module = c("spatial_single_case_GSE251950", "clinical_peritoneal_recurrence_ACRG"),
  enabled = c(isTRUE(cfg$modules$spatial %||% FALSE), isTRUE(cfg$modules$clinical %||% FALSE)),
  default_state = "off",
  reason_off = c(paste("GSE251950 contains a single primary-peritoneal pair (patient GC6). One patient cannot add to",
                       "patient-level evidence, and that section was already reanalysed by an unrefereed preprint.",
                       "The archive is reused read-only from the previous project if the module is ever enabled."),
                 paste("The primary question concerns transcriptional change in established peritoneal lesions, not future",
                       "risk. ACRG has a usable M0 risk set (273 patients, 42 first peritoneal recurrences) but with about",
                       "40-50 events and no second time-to-event cohort it cannot support a prediction claim; running it",
                       "would add a different question rather than evidence for this one.")),
  data_available_if_enabled = c("../project/01_data_raw/GSE251950/GSE251950_RAW.tar (read-only reuse, hash recorded)",
                                "GSE62254 series matrix + Cristescu 2015 Supplementary Data 1 (patient-level, 300/300 matched)"))
o <- pm_out("04_results", "optional", "module_status.csv"); write_csv_atomic(modules, o)
cat("Optional modules (both off by protocol default):\n"); print(modules[, .(module, enabled, default_state)])
step_end(outputs = o)
