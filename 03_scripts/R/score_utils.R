# Program scoring: fixed-direction, equal-weight, within-sample gene-rank scores on a
# frozen common gene universe (plan section 14.3 default). The universe is defined by
# expression-level rules that do not use any site/group information.

suppressPackageStartupMessages(library(Matrix))

# Low-expression filter that never sees the site labels: a gene enters the universe when
# it is expressed above `min_cpm` in at least `min_frac` of ALL samples in the matrix.
frozen_universe <- function(counts, symbols, min_cpm = 1, min_frac = 0.25) {
  cs <- colSums(counts)
  cpm <- t(t(counts) / cs) * 1e6
  keep <- rowMeans(cpm >= min_cpm) >= min_frac
  keep <- keep & !is.na(symbols) & nzchar(symbols)
  data.table(row = which(keep), symbol = symbols[keep])
}

# Within-sample relative rank of every universe gene (average ranks for ties), then the
# equal-weight mean rank of the program genes present. Returns NA when coverage fails.
rank_scores <- function(expr, universe_symbols, program_genes, min_frac = 0.80, min_n = 10) {
  present <- intersect(program_genes, universe_symbols)
  n_def <- length(program_genes)
  if (length(present) / n_def < min_frac || length(present) < min(min_n, n_def))
    return(list(score = rep(NA_real_, ncol(expr)), n_used = length(present), n_defined = n_def, evaluable = FALSE))
  idx <- match(present, universe_symbols)
  n <- length(universe_symbols)
  sc <- apply(expr, 2, function(x) mean(rank(x, ties.method = "average")[idx]) / n)
  list(score = as.numeric(sc), n_used = length(present), n_defined = n_def, evaluable = TRUE)
}

# Collapse a gene x sample matrix to unique symbols (sum of counts) restricted to the universe.
collapse_to_symbols <- function(counts, symbols, universe) {
  m <- counts[universe$row, , drop = FALSE]
  s <- universe$symbol
  if (anyDuplicated(s)) {
    m <- rowsum(as.matrix(m), s)
    s <- rownames(m)
  } else rownames(m) <- s
  list(expr = m, symbols = rownames(m))
}

score_matrix <- function(counts, symbols, programs, min_cpm = 1, min_frac_samples = 0.25,
                         cov_min_frac = 0.80, cov_min_n = 10) {
  u <- frozen_universe(counts, symbols, min_cpm, min_frac_samples)
  cs <- collapse_to_symbols(counts, symbols, u)
  # score on log-CPM (rank invariant to monotone transforms, but keeps values interpretable)
  lib <- colSums(cs$expr); lib[lib == 0] <- 1
  lcpm <- log2(t(t(cs$expr) / lib) * 1e6 + 1)
  res <- lapply(programs, function(g) rank_scores(lcpm, cs$symbols, g, cov_min_frac, cov_min_n))
  scores <- do.call(cbind, lapply(res, `[[`, "score"))
  colnames(scores) <- names(programs); rownames(scores) <- colnames(counts)
  cov <- rbindlist(lapply(names(res), function(k) data.table(program_id = k, n_defined = res[[k]]$n_defined,
                                                             n_used = res[[k]]$n_used, evaluable = res[[k]]$evaluable)))
  list(scores = as.data.table(scores, keep.rownames = "sample"), coverage = cov, universe_size = nrow(cs$expr))
}

# ---- paired patient-level inference -------------------------------------------------
# delta_i = score(PM) - score(primary); each patient contributes exactly one delta.
paired_summary <- function(delta, conf = 0.95, n_boot = 5000, seed = 1) {
  d <- delta[!is.na(delta)]
  n <- length(d)
  if (n < 3) return(data.table(n = n, mean = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_))
  tt <- stats::t.test(d, conf.level = conf)
  if (n_boot >= 100) {
    set.seed(seed)
    bs <- replicate(n_boot, mean(sample(d, n, replace = TRUE)))
    bq <- stats::quantile(bs, c((1 - conf) / 2, 1 - (1 - conf) / 2), names = FALSE)
  } else bq <- c(NA_real_, NA_real_)
  n_pos <- sum(d > 0); n_neg <- sum(d < 0)
  sgn <- stats::binom.test(n_pos, n_pos + n_neg, 0.5)
  wil <- suppressWarnings(stats::wilcox.test(d, conf.int = FALSE))
  data.table(n = n, mean = mean(d), ci_low = tt$conf.int[1], ci_high = tt$conf.int[2], p_value = tt$p.value,
             boot_ci_low = bq[1], boot_ci_high = bq[2], median = stats::median(d),
             n_positive = n_pos, n_negative = n_neg, sign_test_p = sgn$p.value, wilcoxon_p = wil$p.value,
             sd = stats::sd(d), cohen_dz = mean(d) / stats::sd(d))
}

leave_one_out <- function(delta, ids) {
  d <- data.table(id = ids, delta = delta)[!is.na(delta)]
  rbindlist(lapply(seq_len(nrow(d)), function(i) {
    s <- paired_summary(d$delta[-i], n_boot = 0)
    data.table(left_out = d$id[i], mean = s$mean, ci_low = s$ci_low, ci_high = s$ci_high, p_value = s$p_value)
  }))
}
