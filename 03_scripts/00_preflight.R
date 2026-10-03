# Stage 0: environment probe and protection snapshot of the old project.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
step_begin("00_preflight")

cfg <- pm_config()
need <- unique(unlist(cfg$packages))
have <- rownames(installed.packages())
miss <- setdiff(need, have)

disk <- tryCatch(system2("df", c("-k", shQuote(pm_project_root())), stdout = TRUE), error = function(e) NA)
free_gb <- tryCatch(as.numeric(strsplit(trimws(disk[2]), "\\s+")[[1]][4]) / 1024^2, error = function(e) NA)
mem_gb <- tryCatch(as.numeric(system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE)) / 1024^3, error = function(e) NA)
ncpu <- parallel::detectCores(logical = TRUE)
tool <- function(p) if (nzchar(p) && file.exists(p)) p else NA_character_
pandoc <- tool(cfg$tools$pandoc); soffice <- tool(cfg$tools$soffice)

# Old project snapshot: relative path, size and mtime (cheap, detects any modification).
old <- pm_old_project()
snap <- if (dir.exists(old)) {
  f <- list.files(old, recursive = TRUE, all.files = TRUE, full.names = TRUE, no.. = TRUE)
  f <- f[!grepl("/renv/(library|staging|sandbox)/", f)]
  f <- f[!grepl("(^|/)\\.DS_Store$", f)]   # macOS Finder metadata; touched by the OS, carries no research content
  fi <- file.info(f)
  data.table(path = sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", old), "/"), "", f),
             size = as.numeric(fi$size), mtime = format(fi$mtime, "%Y-%m-%d %H:%M:%OS3"))
} else data.table(path = character(), size = numeric(), mtime = character())
snap_file <- pm_out("00_admin", "old_project_snapshot.csv")
if (!file.exists(snap_file)) write_csv_atomic(snap, snap_file)   # baseline written once; later runs compare
base <- fread(snap_file, colClasses = c(path = "character", size = "numeric", mtime = "character"))
changed <- if (nrow(snap)) {
  m <- merge(base, snap, by = "path", all = TRUE, suffixes = c(".base", ".now"))
  m[is.na(size.base) | is.na(size.now) | size.base != size.now | mtime.base != mtime.now]
} else base[0]

rep <- c(
  "# Automated environment report", "",
  sprintf("- generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("- R: %s; platform %s", R.version.string, R.version$platform),
  sprintf("- logical CPUs: %d; configured workers: %d", ncpu, cfg$resources$workers),
  sprintf("- memory GB: %.1f", mem_gb),
  sprintf("- free disk GB at project volume: %.1f (reserve %d GB)", free_gb, cfg$resources$reserve_disk_gb),
  sprintf("- pandoc: %s", pandoc %||% "not found"),
  sprintf("- soffice: %s", soffice %||% "not found"),
  sprintf("- required packages missing: %s", if (length(miss)) paste(miss, collapse = ", ") else "none"),
  sprintf("- old project files in snapshot: %d; changed since baseline: %d", nrow(base), nrow(changed))
)
out <- pm_out("00_admin", "environment_report_auto.md")
write_text_atomic(rep, out)
if (nrow(changed)) {
  write_csv_atomic(changed, pm_out("00_admin", "old_project_changes_detected.csv"))
  warning("Old project files differ from baseline snapshot; see old_project_changes_detected.csv")
}
if (length(miss)) stop("Missing required packages: ", paste(miss, collapse = ", "))
if (!is.na(free_gb) && free_gb < cfg$resources$reserve_disk_gb) stop("Free disk below reserve")
step_end(outputs = c(out, snap_file))
