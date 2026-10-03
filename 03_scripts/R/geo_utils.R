# GEO metadata helpers: parse SOFT family files into one-row-per-GSM tables while
# keeping every raw field, and a bounded, verified downloader used by the download stage.

parse_soft_samples <- function(soft_gz) {
  con <- gzfile(soft_gz, "r"); on.exit(close(con))
  x <- readLines(con, warn = FALSE, encoding = "UTF-8")
  starts <- grep("^\\^SAMPLE = ", x)
  if (!length(starts)) return(data.table())
  ends <- c(starts[-1] - 1L, length(x))
  rows <- lapply(seq_along(starts), function(i) {
    blk <- x[starts[i]:ends[i]]
    blk <- blk[startsWith(blk, "!Sample_")]
    key <- sub("^!Sample_([^ ]+) = .*$", "\\1", blk)
    val <- sub("^!Sample_[^ ]+ = ?", "", blk)
    out <- list(gsm = sub("^\\^SAMPLE = ", "", x[starts[i]]))
    # characteristics: "key: value" -> char_<key>
    ch <- key == "characteristics_ch1"
    if (any(ch)) {
      ck <- tolower(trimws(sub(":.*$", "", val[ch])))
      ck <- gsub("[^a-z0-9]+", "_", ck)
      cv <- trimws(sub("^[^:]*:", "", val[ch]))
      for (j in seq_along(ck)) {
        nm <- paste0("char_", ck[j])
        out[[nm]] <- if (is.null(out[[nm]])) cv[j] else paste(out[[nm]], cv[j], sep = " | ")
      }
      out$characteristics_raw <- paste(val[ch], collapse = " || ")
    }
    for (k in unique(key[!ch])) out[[k]] <- paste(val[key == k], collapse = " || ")
    out
  })
  rbindlist(rows, fill = TRUE)
}

parse_soft_series <- function(soft_gz) {
  con <- gzfile(soft_gz, "r"); on.exit(close(con))
  x <- readLines(con, warn = FALSE, encoding = "UTF-8")
  s <- grep("^\\^SERIES = ", x); e <- grep("^\\^(PLATFORM|SAMPLE) = ", x)
  e <- if (length(e)) min(e[e > s]) - 1L else length(x)
  blk <- x[s:e]; blk <- blk[startsWith(blk, "!Series_")]
  key <- sub("^!Series_([^ ]+) = .*$", "\\1", blk)
  val <- sub("^!Series_[^ ]+ = ?", "", blk)
  data.table(field = key, value = val)
}

# Download one file with bounded retries, resumable transfer, content-type sanity
# checks and SHA-256. Returns a manifest row. Never overwrites a verified file whose
# hash matches the expected value.
download_verified <- function(url, dest, expected_sha256 = NA_character_, max_time = 3600,
                              attempts = 3, min_bytes = 1, allow_html = FALSE) {
  assert_writable(dest, allow_raw = TRUE)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && !is.na(expected_sha256) && identical(sha256_file(dest), expected_sha256)) {
    return(data.table(url = url, local_path = dest, bytes = file.info(dest)$size, sha256 = expected_sha256,
                      status = "present_verified", attempts = 0L, time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
  }
  part <- paste0(dest, ".part")
  ok <- FALSE; a <- 0L; msg <- NA_character_
  while (!ok && a < attempts) {
    a <- a + 1L
    st <- system2("curl", c("-L", "-f", "-sS", "--connect-timeout", "30", "--max-time", as.character(max_time),
                            "--speed-time", "120", "--speed-limit", "1024", "-C", "-", "-o", shQuote(part), shQuote(url)),
                  stdout = TRUE, stderr = TRUE)
    code <- attr(st, "status") %||% 0L
    ok <- identical(as.integer(code), 0L) && file.exists(part) && file.info(part)$size >= min_bytes
    if (!ok) { msg <- paste(st, collapse = " "); Sys.sleep(5 * a) }
  }
  if (!ok) return(data.table(url = url, local_path = dest, bytes = NA_real_, sha256 = NA_character_,
                             status = paste("failed:", msg), attempts = a, time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
  head <- readBin(part, "raw", n = 512)
  is_html <- any(vapply(c("<!doctype html", "<!DOCTYPE html", "<!DOCTYPE HTML", "<html", "<HTML"),
                        function(p) length(grepRaw(p, head, fixed = TRUE)) > 0, TRUE))
  if (!allow_html && is_html) {
    unlink(part)
    return(data.table(url = url, local_path = dest, bytes = NA_real_, sha256 = NA_character_,
                      status = "failed: HTML page received", attempts = a, time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
  }
  if (grepl("\\.gz$", dest)) {
    gz_ok <- identical(system2("gzip", c("-t", shQuote(part)), stdout = FALSE, stderr = FALSE), 0L)
    if (!gz_ok) stop("Downloaded gzip failed integrity test: ", url)
  }
  if (grepl("\\.tar$", dest)) {
    tar_ok <- identical(system2("tar", c("-tf", shQuote(part)), stdout = FALSE, stderr = FALSE), 0L)
    if (!tar_ok) stop("Downloaded tar failed listing test: ", url)
  }
  if (grepl("\\.zip$", dest)) {
    zip_ok <- identical(system2("unzip", c("-tq", shQuote(part)), stdout = FALSE, stderr = FALSE), 0L)
    if (!zip_ok) stop("Downloaded zip failed integrity test: ", url)
  }
  h <- sha256_file(part)
  if (!is.na(expected_sha256) && !identical(h, expected_sha256)) stop("SHA-256 mismatch for ", url)
  file.rename(part, dest)
  Sys.chmod(dest, mode = "0444")
  data.table(url = url, local_path = dest, bytes = file.info(dest)$size, sha256 = h,
             status = "downloaded", attempts = a, time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
}
