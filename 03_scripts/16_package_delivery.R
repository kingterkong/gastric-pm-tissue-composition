# Stage 16: build the delivery packages and their manifest.
# Raw expression matrices are never bundled; the manifest carries their URLs and checksums.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
cfg <- pm_config()
step_begin("16_package_delivery", inputs = pm_in("00_admin", "final_checklist.csv"))

DATE <- "20260926"
zip_dir <- pm_out("08_delivery", "x"); zip_dir <- dirname(zip_dir)
stage <- file.path(tempdir(), paste0("pm_delivery_", DATE))
unlink(stage, recursive = TRUE); dir.create(stage, recursive = TRUE)

copy_into <- function(files, subdir) {
  d <- file.path(stage, subdir); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  files <- files[file.exists(files)]
  file.copy(files, d, overwrite = TRUE, recursive = TRUE)
  files
}
zip_it <- function(subdir, out) {
  wd <- getwd(); on.exit(setwd(wd))
  setwd(stage)
  unlink(out)
  system2("zip", c("-rq", shQuote(out), shQuote(subdir)), stdout = TRUE, stderr = TRUE)
  out
}

JS <- pm_in("06_manuscript", "journal_submission")
IR <- pm_in("06_manuscript", "internal_review")

# 1. journal submission materials
copy_into(list.files(JS, full.names = TRUE), "journal_submission")
copy_into(list.files(pm_in("05_figures", "main"), full.names = TRUE), "journal_submission/figures")
copy_into(list.files(pm_in("04_results", "manuscript"), pattern = "^Table.*\\.csv$|main_tables\\.xlsx$", full.names = TRUE), "journal_submission/tables")
copy_into(c(pm_in("07_references", "references.bib"), pm_in("07_references", "reference_verification.csv"),
            pm_in("07_references", "claim_citation_map.csv")), "journal_submission/references")
z1 <- zip_it("journal_submission", file.path(zip_dir, sprintf("journal_submission_materials_%s.zip", DATE)))

# 2. internal review materials
copy_into(list.files(IR, full.names = TRUE), "internal_review")
copy_into(c(pm_in("00_admin", "evidence_status.csv"), pm_in("00_admin", "claims_ledger.csv"),
            pm_in("00_admin", "hypothesis_registry.csv"), pm_in("00_admin", "deviations.md"),
            pm_in("00_admin", "final_checklist.md"), pm_in("00_admin", "final_checklist.csv")), "internal_review/records")
z2 <- zip_it("internal_review", file.path(zip_dir, sprintf("internal_review_%s.zip", DATE)))

# 3. reproducibility bundle: code, config, frozen definitions, manifests, logs, reports
copy_into(list.files(pm_src("03_scripts"), full.names = TRUE, recursive = TRUE, include.dirs = FALSE), "reproducibility/03_scripts")
copy_into(c(pm_src("run_all.R"), pm_src("config.yml"), pm_src("renv.lock"), pm_src("README.md"), pm_src("RESUME.md")), "reproducibility")
copy_into(list.files(pm_in("02_data_processed", "frozen_programs"), full.names = TRUE), "reproducibility/frozen_programs")
copy_into(c(pm_in("00_admin", "protocol_v1.md"), pm_in("00_admin", "protocol_lock.json"),
            pm_in("00_admin", "data_manifest.csv"), pm_in("00_admin", "sample_crosswalk.csv"),
            pm_in("00_admin", "overlap_register.csv"), pm_in("00_admin", "eligibility_flow.csv"),
            pm_in("00_admin", "decision_log.md"), pm_in("00_admin", "topic_decision.md"),
            pm_in("00_admin", "novelty_matrix.csv"), pm_in("00_admin", "environment_report.md"),
            pm_in("00_admin", "environment_report_auto.md"), pm_in("00_admin", "access_log.csv"),
            pm_in("00_admin", "raw_file_inventory.csv"), pm_in("00_admin", "analysis_log.md"),
            pm_in("00_admin", "STATE.json")), "reproducibility/00_admin")
copy_into(list.files(pm_in("04_results"), full.names = TRUE, recursive = TRUE, pattern = "\\.csv$"), "reproducibility/04_results")
copy_into(list.files(pm_in("_rebuilds"), full.names = TRUE, pattern = "_report\\.md$|_comparison\\.csv$|_check\\.csv$"), "reproducibility/rebuild")
copy_into(list.files(pm_in("00_admin", "execution_logs"), full.names = TRUE, pattern = "\\.log$"), "reproducibility/execution_logs")
z3 <- zip_it("reproducibility", file.path(zip_dir, sprintf("reproducibility_bundle_%s.zip", DATE)))

zips <- c(z1, z2, z3)
man <- rbindlist(lapply(zips, function(z) data.table(
  file = basename(z), purpose = fcase(grepl("journal", z), "materials prepared for journal submission",
                                      grepl("internal", z), "internal assessment; not for submission",
                                      default = "code, frozen definitions, manifests and logs for reproduction"),
  format = "zip", bytes = file.info(z)$size, sha256 = sha256_file(z),
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), version = "current", check_status = "built and hashed")))
extra <- data.table(file = "START_HERE_中文交付说明.md", purpose = "two-minute delivery note in Chinese",
                    format = "markdown", bytes = file.info(pm_in("08_delivery", "START_HERE_中文交付说明.md"))$size,
                    sha256 = sha256_file(pm_in("08_delivery", "START_HERE_中文交付说明.md")),
                    generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), version = "current", check_status = "written")
man <- rbind(man, extra, fill = TRUE)
man[, note := "raw expression matrices are not bundled; their URLs and SHA-256 are in 00_admin/data_manifest.csv"]
o <- pm_out("08_delivery", "delivery_manifest.csv"); write_csv_atomic(man, o)
print(man[, .(file, mb = round(bytes / 1024^2, 2), sha = substr(sha256, 1, 12))])
unlink(stage, recursive = TRUE)
step_end(outputs = c(zips, o))
