# Stage 1: parse public GEO metadata (SOFT family files) into sample tables, merge
# per-source access logs, and inventory every raw file with size and SHA-256.
# Reads only 01_data_raw; writes 02_data_processed/metadata and 00_admin tables.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
source(pm_src("03_scripts", "R", "geo_utils.R"))

softs <- list.files(pm_raw(), pattern = "_family\\.soft\\.gz$", recursive = TRUE, full.names = TRUE)
softs <- softs[!grepl("/literature/", softs)]
step_begin("01_literature_and_metadata", inputs = softs)

outs <- character()
for (f in softs) {
  gse <- sub("_family\\.soft\\.gz$", "", basename(f))
  smp <- parse_soft_samples(f)
  smp[, dataset_id := gse]
  setcolorder(smp, c("dataset_id", "gsm"))
  # drop submitter contact fields (not needed for analysis)
  drop <- grep("^contact_", names(smp), value = TRUE)
  if (length(drop)) smp[, (drop) := NULL]
  o <- pm_out("02_data_processed", "metadata", paste0(gse, "_geo_samples.csv"))
  write_csv_atomic(smp, o); outs <- c(outs, o)
  ser <- parse_soft_series(f)
  ser <- ser[!grepl("^contact_", field)]
  o2 <- pm_out("02_data_processed", "metadata", paste0(gse, "_geo_series.csv"))
  write_csv_atomic(ser, o2); outs <- c(outs, o2)
}

# Merge access logs written during the audit (one file per audit source).
acc <- list.files(pm_src("00_admin", "audit_notes"), pattern = "^access_.*\\.csv$", full.names = TRUE)
acc_dt <- rbindlist(lapply(acc, function(a) {
  d <- tryCatch(fread(a, colClasses = "character", fill = TRUE), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  d[, audit_source := sub("^access_(.*)\\.csv$", "\\1", basename(a))]
  d
}), fill = TRUE)
if (nrow(acc_dt)) {
  o <- pm_out("00_admin", "access_log.csv"); write_csv_atomic(acc_dt, o); outs <- c(outs, o)
}

# Raw file inventory (all files under 01_data_raw, excluding VCS internals).
raw_files <- list.files(pm_raw(), recursive = TRUE, full.names = TRUE, all.files = FALSE)
raw_files <- raw_files[!grepl("/\\.git/", raw_files) & !grepl("\\.part$", raw_files)]
inv <- data.table(
  local_path = sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", pm_project_root()), "/"), "", raw_files),
  bytes = file.info(raw_files)$size,
  sha256 = vapply(raw_files, sha256_file, "")
)
o <- pm_out("00_admin", "raw_file_inventory.csv"); write_csv_atomic(inv, o); outs <- c(outs, o)
step_end(outputs = outs)
