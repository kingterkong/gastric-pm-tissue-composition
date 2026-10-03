# Stage 3: download the analysis inputs declared in config.yml (downloads:) with bounded
# retries, verify content and SHA-256, and write 00_admin/data_manifest.csv. Files already
# present with the recorded hash are not re-downloaded. Read-only reuse of an old-project
# file is recorded (path + hash) without copying or writing back.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
source(pm_src("03_scripts", "R", "geo_utils.R"))
cfg <- pm_config()
dl <- cfg$downloads
step_begin("03_download_inputs", params = list(n_files = length(dl)))

prev_manifest <- pm_src("00_admin", "data_manifest.csv")
prev <- if (file.exists(prev_manifest)) fread(prev_manifest, colClasses = "character") else data.table()

rows <- list()
for (d in dl) {
  if (!isTRUE(d$enabled %||% TRUE)) {
    rows[[length(rows) + 1]] <- data.table(accession = d$accession, file_name = basename(d$dest), official_url = d$url,
                                           local_path = d$dest, status = paste("disabled:", d$disabled_reason %||% "not needed"),
                                           purpose = d$purpose, role = d$role)
    next
  }
  if (!is.null(d$reuse_readonly)) {
    src <- normalizePath(file.path(pm_project_root(), d$reuse_readonly), mustWork = TRUE)
    h <- sha256_file(src)
    if (!is.null(d$sha256) && !identical(h, d$sha256)) stop("Read-only reuse hash mismatch: ", src)
    rows[[length(rows) + 1]] <- data.table(accession = d$accession, file_name = basename(src), official_url = d$url,
                                           local_path = d$reuse_readonly, bytes = file.info(src)$size, sha256 = h,
                                           status = "reused_readonly", purpose = d$purpose, role = d$role,
                                           download_time = NA_character_)
    next
  }
  dest <- pm_raw(d$dest)
  exp_h <- d$sha256 %||% NA_character_
  if (is.na(exp_h) && nrow(prev) && "local_path" %in% names(prev)) {
    hit <- prev[local_path == file.path("01_data_raw", d$dest) & nzchar(sha256)]
    if (nrow(hit)) exp_h <- hit$sha256[1]
  }
  if (file.exists(dest) && is.na(exp_h)) {
    r <- data.table(url = d$url, local_path = dest, bytes = file.info(dest)$size, sha256 = sha256_file(dest),
                    status = "present_unrecorded_hash", attempts = 0L, time = NA_character_)
  } else {
    r <- download_verified(d$url, dest, expected_sha256 = exp_h, max_time = d$max_time %||% cfg$resources$large_file_attempt_timeout_s,
                           attempts = cfg$resources$download_attempts, min_bytes = d$min_bytes %||% 1)
  }
  if (grepl("^failed", r$status)) stop("Download failed for ", d$url, ": ", r$status)
  rows[[length(rows) + 1]] <- data.table(accession = d$accession, file_name = basename(dest), official_url = d$url,
                                         local_path = file.path("01_data_raw", d$dest), bytes = r$bytes, sha256 = r$sha256,
                                         status = r$status, purpose = d$purpose, role = d$role, download_time = r$time)
}
man <- rbindlist(rows, fill = TRUE)
extra <- rbindlist(lapply(cfg$downloads, function(d) data.table(
  file_name = basename(d$dest %||% d$reuse_readonly), source_study = d$source_study %||% NA_character_,
  release = d$release %||% NA_character_, compression = if (grepl("\\.(gz|tar|zip)$", d$dest %||% d$reuse_readonly)) sub(".*\\.", "", d$dest %||% d$reuse_readonly) else "none",
  data_type = d$data_type %||% NA_character_, patient_key_source = d$patient_key_source %||% NA_character_,
  access_status = d$access_status %||% "public", redistribution_note = d$redistribution_note %||% "public GEO/journal file; cite source; not bundled in delivery"
)), fill = TRUE)
man <- cbind(man, extra[, !"file_name"])
o <- pm_out("00_admin", "data_manifest.csv")
write_csv_atomic(man, o)
step_end(outputs = o)
