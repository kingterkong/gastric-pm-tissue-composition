# Shared infrastructure for project_pm scripts: paths, write whitelist, hashing,
# atomic writes, step logging and resume markers. Sourced by every stage script.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(digest)
  library(yaml)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# ---- roots -------------------------------------------------------------------
pm_project_root <- function() {
  r <- Sys.getenv("PM_PROJECT_ROOT", "")
  if (!nzchar(r)) {
    here <- normalizePath(getwd(), mustWork = TRUE)
    while (!file.exists(file.path(here, "run_all.R")) || !file.exists(file.path(here, "config.yml"))) {
      parent <- dirname(here)
      if (identical(parent, here)) stop("Cannot locate project_pm root (run_all.R + config.yml)")
      here <- parent
    }
    r <- here
  }
  normalizePath(r, mustWork = TRUE)
}

pm_output_root <- function() {
  o <- Sys.getenv("PM_OUTPUT_ROOT", "")
  if (!nzchar(o)) o <- pm_project_root()
  dir.create(o, recursive = TRUE, showWarnings = FALSE)
  normalizePath(o, mustWork = TRUE)
}

# Raw data are always read from the project raw directory (read-only for all
# stages except the download stage), even during an isolated rebuild.
pm_raw <- function(...) file.path(pm_project_root(), "01_data_raw", ...)
pm_old_project <- function() normalizePath(file.path(pm_project_root(), "..", "project"), mustWork = FALSE)

pm_out <- function(...) {
  p <- file.path(pm_output_root(), ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}
pm_in <- function(...) file.path(pm_output_root(), ...)          # read an output of an earlier stage
pm_src <- function(...) file.path(pm_project_root(), ...)          # read project source (scripts, config, curated inputs)

# Read a record that a rebuild does not regenerate (for example the download manifest):
# prefer the current output root, fall back to the project root, and fail loudly if neither has it.
pm_record <- function(...) {
  a <- file.path(pm_output_root(), ...)
  if (file.exists(a)) return(a)
  b <- file.path(pm_project_root(), ...)
  if (file.exists(b)) return(b)
  stop("record not found in output root or project root: ", file.path(...))
}

# ---- write whitelist ------------------------------------------------------------
.real <- function(p) {
  # resolve the deepest existing ancestor (symlinks included) and re-append the rest
  p <- path.expand(p)
  tail <- character()
  while (!file.exists(p) && !identical(dirname(p), p)) {
    tail <- c(basename(p), tail)
    p <- dirname(p)
  }
  r <- normalizePath(p, mustWork = TRUE)
  if (length(tail)) r <- do.call(file.path, as.list(c(r, tail)))
  r
}
.under <- function(p, root) {
  p <- .real(p); root <- .real(root)
  startsWith(paste0(p, "/"), paste0(root, "/"))
}

assert_writable <- function(path, allow_raw = FALSE) {
  rp <- .real(path)
  old <- pm_old_project()
  if (file.exists(old) && .under(rp, old)) stop("Refusing to write into old project: ", rp)
  raw <- pm_raw()
  if (.under(rp, raw) && !allow_raw) stop("Refusing to write into raw data directory outside download stage: ", rp)
  ok <- .under(rp, pm_output_root()) || (allow_raw && .under(rp, raw))
  if (!ok) stop("Path outside write whitelist: ", rp)
  # in an isolated rebuild the output root must not be the working project tree itself
  invisible(rp)
}

check_output_root <- function() {
  o <- pm_output_root(); proj <- pm_project_root()
  old <- pm_old_project()
  if (file.exists(old) && .under(o, old)) stop("--output-root may not point into old project")
  if (.under(o, pm_raw())) stop("--output-root may not point into 01_data_raw")
  if (!(identical(.real(o), .real(proj)) || .under(o, file.path(proj, "_rebuilds"))))
    stop("--output-root must be project_pm or a directory under project_pm/_rebuilds")
  invisible(o)
}

# ---- hashing and atomic writes ------------------------------------------------
sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
}
sha256_obj <- function(x) digest::digest(x, algo = "sha256")

atomic_write <- function(path, writer, allow_raw = FALSE) {
  assert_writable(path, allow_raw = allow_raw)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writer(tmp)
  if (!file.exists(tmp) || file.info(tmp)$size == 0) stop("Atomic write produced empty file: ", path)
  ok <- file.rename(tmp, path)
  if (!ok) stop("Atomic rename failed: ", path)
  invisible(path)
}
write_csv_atomic <- function(x, path, ...) atomic_write(path, function(tmp) data.table::fwrite(x, tmp, ...))
write_json_atomic <- function(x, path) atomic_write(path, function(tmp) jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", na = "null"))
write_text_atomic <- function(lines, path) atomic_write(path, function(tmp) writeLines(enc2utf8(lines), tmp, useBytes = TRUE))
save_rds_atomic <- function(x, path) atomic_write(path, function(tmp) saveRDS(x, tmp))

# ---- config ----------------------------------------------------------------------
pm_config <- function() yaml::read_yaml(pm_src("config.yml"))

# ---- step logging and resume markers ----------------------------------------
.step <- new.env()

step_begin <- function(step, inputs = character(), params = list(), seed = NULL, uses_protocol = FALSE) {
  check_output_root()
  .step$name <- step
  .step$uses_protocol <- uses_protocol
  .step$start <- Sys.time()
  .step$inputs <- inputs
  .step$params <- params
  .step$seed <- seed
  if (!is.null(seed)) set.seed(seed)
  logdir <- pm_out("00_admin", "execution_logs")
  .step$log <- file.path(logdir, sprintf("%s_%s.log", step, format(.step$start, "%Y%m%d_%H%M%S")))
  cat(sprintf("[%s] BEGIN %s\n", format(.step$start), step))
  invisible(TRUE)
}

.script_hash <- function(step) {
  f <- list.files(pm_src("03_scripts"), pattern = paste0("^", step, ".*\\.R$"), full.names = TRUE)
  u <- list.files(pm_src("03_scripts", "R"), pattern = "\\.R$", full.names = TRUE)
  sha256_obj(vapply(sort(c(f, u)), sha256_file, ""))
}

step_end <- function(outputs = character(), status = "success", notes = NULL) {
  end <- Sys.time()
  ins <- .step$inputs; outs <- outputs
  in_h <- setNames(vapply(ins, sha256_file, ""), ins)
  out_h <- setNames(vapply(outs, sha256_file, ""), outs)
  missing_out <- outs[is.na(out_h)]
  if (length(missing_out) && status == "success") stop("Declared outputs missing: ", paste(missing_out, collapse = ", "))
  cfg_h <- sha256_file(pm_src("config.yml"))
  lock <- pm_in("00_admin", "protocol_lock.json")
  prot_h <- if (isTRUE(.step$uses_protocol) && file.exists(lock)) sha256_file(lock) else NA_character_
  rel <- function(p) sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", pm_output_root()), "/?"), "", p)
  rec <- list(
    step = .step$name, status = status,
    command = paste(commandArgs(), collapse = " "),
    start = format(.step$start, "%Y-%m-%dT%H:%M:%S%z"), end = format(end, "%Y-%m-%dT%H:%M:%S%z"),
    seconds = round(as.numeric(difftime(end, .step$start, units = "secs")), 1),
    seed = .step$seed %||% NA, params = .step$params,
    inputs = as.list(setNames(unname(in_h), vapply(names(in_h), rel, ""))),
    outputs = as.list(setNames(unname(out_h), vapply(names(out_h), rel, ""))),
    script_hash = .script_hash(.step$name), config_hash = cfg_h, protocol_lock_hash = prot_h,
    r_version = R.version.string, notes = notes %||% NA
  )
  si <- capture.output(sessionInfo())
  write_text_atomic(c(toJSON(rec, auto_unbox = TRUE, pretty = TRUE, na = "null"), "", "## sessionInfo()", si), .step$log)
  if (status == "success") write_json_atomic(rec, pm_out("00_admin", "execution_logs", "markers", paste0(.step$name, ".json")))
  cat(sprintf("[%s] END %s (%s, %.1fs)\n", format(end), .step$name, status, rec$seconds))
  invisible(rec)
}

# A step is current when its marker exists, all recorded input/output hashes
# still match, and script/config/protocol hashes are unchanged.
step_is_current <- function(step) {
  m <- pm_in("00_admin", "execution_logs", "markers", paste0(step, ".json"))
  if (!file.exists(m)) return(FALSE)
  rec <- jsonlite::read_json(m)
  chk <- function(lst) all(vapply(names(lst), function(p) {
    fp <- if (startsWith(p, "/")) p else file.path(pm_output_root(), p)
    identical(sha256_file(fp), lst[[p]])
  }, TRUE))
  lock <- pm_in("00_admin", "protocol_lock.json")
  prot_ok <- is.null(rec$protocol_lock_hash) || identical(rec$protocol_lock_hash, sha256_file(lock))
  identical(rec$status, "success") && chk(rec$inputs) && chk(rec$outputs) &&
    identical(rec$script_hash, .script_hash(step)) &&
    identical(rec$config_hash, sha256_file(pm_src("config.yml"))) && prot_ok
}

invalidate_step <- function(step) {
  m <- pm_in("00_admin", "execution_logs", "markers", paste0(step, ".json"))
  if (file.exists(m)) { assert_writable(m); unlink(m) }
}

stopf <- function(...) stop(sprintf(...), call. = FALSE)
