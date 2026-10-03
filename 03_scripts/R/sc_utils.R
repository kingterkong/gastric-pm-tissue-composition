# Single-cell helpers: readers, HGNC symbol harmonisation, per-sample robust QC,
# marker-based lineage annotation, expression-based CNV scoring for malignant calls,
# and pseudobulk aggregation. Program genes are never used to define lineages or
# malignant status (lineage markers are fixed in LINEAGE_MARKERS below).

suppressPackageStartupMessages({ library(Matrix); library(data.table) })

LINEAGE_MARKERS <- list(
  epithelial  = c("EPCAM", "KRT8", "KRT18", "KRT19", "CLDN4", "CDH1", "ELF3"),
  fibroblast  = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "FBLN1"),
  mural       = c("RGS5", "PDGFRB", "NOTCH3", "MYH11", "CSPG4", "MCAM"),
  endothelial = c("PECAM1", "VWF", "CDH5", "PLVAP", "CLDN5", "EMCN"),
  mesothelial = c("MSLN", "UPK3B", "CALB2", "WT1", "ITLN1", "PRG4", "LRRN4"),
  myeloid     = c("LYZ", "CD68", "CD14", "C1QA", "AIF1", "FCER1G", "CSF1R"),
  t_nk        = c("CD3E", "CD3D", "CD2", "NKG7", "GNLY", "CD247"),
  b_cell      = c("MS4A1", "CD79A", "CD19", "BANK1"),
  plasma      = c("MZB1", "JCHAIN", "SDC1", "TNFRSF17"),
  mast        = c("TPSAB1", "CPA3", "KIT", "HPGDS")
)

# ---- symbol harmonisation ------------------------------------------------------------
hgnc_maps <- function(hgnc_file) {
  h <- fread(hgnc_file, sep = "\t", quote = "", colClasses = "character",
             select = c("symbol", "status", "prev_symbol", "alias_symbol", "ensembl_gene_id"))
  h <- h[status == "Approved"]
  split_col <- function(col) { d <- h[nzchar(get(col)), .(symbol, x = get(col))]; d[, .(x = unlist(strsplit(x, "\\|"))), by = symbol] }
  prev <- split_col("prev_symbol"); alias <- split_col("alias_symbol")
  prev <- prev[!x %in% h$symbol]; prev <- prev[, if (.N == 1) .SD, by = x]            # unique only
  alias <- alias[!x %in% h$symbol & !x %in% prev$x]; alias <- alias[, if (.N == 1) .SD, by = x]
  list(approved = h$symbol, ens = h[nzchar(ensembl_gene_id), .(ensembl_gene_id, symbol)], prev = prev, alias = alias)
}
harmonise_symbols <- function(symbols, ensembl = NULL, maps) {
  out <- rep(NA_character_, length(symbols)); how <- rep("unmapped", length(symbols))
  if (!is.null(ensembl)) {
    i <- match(ensembl, maps$ens$ensembl_gene_id); ok <- !is.na(i)
    out[ok] <- maps$ens$symbol[i[ok]]; how[ok] <- "ensembl_id"
  }
  j <- is.na(out) & symbols %in% maps$approved; out[j] <- symbols[j]; how[j] <- "approved_symbol"
  k <- is.na(out); p <- match(symbols[k], maps$prev$x); out[k][!is.na(p)] <- maps$prev$symbol[p[!is.na(p)]]
  how[k][!is.na(p)] <- "previous_symbol"
  k <- is.na(out); a <- match(symbols[k], maps$alias$x); out[k][!is.na(a)] <- maps$alias$symbol[a[!is.na(a)]]
  how[k][!is.na(a)] <- "alias_unique"
  k <- is.na(out); out[k] <- symbols[k]
  data.table(original = symbols, ensembl = if (is.null(ensembl)) NA_character_ else ensembl, symbol = out, how = how)
}
# collapse duplicated harmonised symbols by summing counts
collapse_rows <- function(m, sym) {
  if (!anyDuplicated(sym)) { rownames(m) <- sym; return(m) }
  f <- factor(sym, levels = unique(sym))
  agg <- Matrix::sparse.model.matrix(~ 0 + f, transpose = TRUE)
  r <- agg %*% m
  rownames(r) <- levels(f)
  as(r, "CsparseMatrix")
}

# ---- readers -------------------------------------------------------------------------
read_10x_triplet <- function(mtx, features, barcodes) {
  m <- Matrix::readMM(mtx)
  f <- fread(features, header = FALSE)
  b <- fread(barcodes, header = FALSE)$V1
  m <- as(m, "CsparseMatrix")
  list(counts = m, ensembl = f$V1, symbol = if (ncol(f) >= 2) f$V2 else f$V1, barcodes = b)
}
read_dense_csv_gz <- function(f) {
  d <- fread(cmd = sprintf("gzip -dc %s", shQuote(f)), header = TRUE)
  genes <- d[[1]]; d[, 1 := NULL]
  m <- as(as.matrix(d), "CsparseMatrix")
  rownames(m) <- genes
  list(counts = m, symbol = genes, ensembl = NULL, barcodes = colnames(d))
}

# ---- QC --------------------------------------------------------------------------------
qc_cells <- function(m, sample_id, mt_genes, min_genes = 200, max_mt_hard = 25, nmad = 3) {
  n_counts <- Matrix::colSums(m); n_genes <- Matrix::colSums(m > 0)
  mt <- if (length(mt_genes)) Matrix::colSums(m[intersect(mt_genes, rownames(m)), , drop = FALSE]) / pmax(n_counts, 1) * 100 else rep(0, ncol(m))
  lc <- log1p(n_counts); lg <- log1p(n_genes)
  mad_out <- function(x) abs(x - median(x)) > nmad * mad(x)
  mt_thr <- min(max_mt_hard, median(mt) + nmad * mad(mt))
  mt_thr <- max(mt_thr, 5)
  keep <- n_genes >= min_genes & mt <= mt_thr & !mad_out(lc) & !mad_out(lg)
  data.table(sample_id = sample_id, barcode = colnames(m), n_counts = n_counts, n_genes = n_genes, pct_mt = mt,
             mt_threshold = mt_thr, qc_pass = keep)
}

# ---- normalisation and lineage annotation ------------------------------------------
lognorm <- function(m, sf = 1e4) {
  cs <- Matrix::colSums(m); cs[cs == 0] <- 1
  x <- m %*% Matrix::Diagonal(x = sf / cs)
  x@x <- log1p(x@x); dimnames(x) <- dimnames(m); x
}
marker_scores <- function(lx, markers = LINEAGE_MARKERS) {
  sapply(markers, function(g) {
    g <- intersect(g, rownames(lx))
    if (length(g) < 2) return(rep(NA_real_, ncol(lx)))
    Matrix::colMeans(lx[g, , drop = FALSE])
  })
}
# cluster-level lineage calls: per-sample Seurat graph clustering on HVG PCA; each cluster
# receives the lineage with the highest mean marker score if it exceeds the runner-up by a
# margin; otherwise "uncertain". Per-cell marker scores are retained for sensitivity.
annotate_lineages <- function(counts, lx, n_hvg = 2000, n_pcs = 30, k = 20, resolution = 1, margin = 0.15, seed = 1) {
  suppressPackageStartupMessages(library(Seurat))
  set.seed(seed)
  so <- CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
  so <- NormalizeData(so, verbose = FALSE)
  so <- FindVariableFeatures(so, nfeatures = n_hvg, verbose = FALSE)
  so <- ScaleData(so, verbose = FALSE)
  npc <- min(n_pcs, ncol(so) - 2)
  so <- RunPCA(so, npcs = npc, verbose = FALSE, seed.use = seed)
  so <- FindNeighbors(so, dims = seq_len(npc), k.param = min(k, ncol(so) - 1), verbose = FALSE)
  so <- FindClusters(so, resolution = resolution, random.seed = seed, verbose = FALSE)
  cl <- as.integer(as.character(so$seurat_clusters))
  sc <- marker_scores(lx)
  cl_means <- rowsum(sc, cl, na.rm = TRUE) / as.vector(table(cl))
  call <- apply(cl_means, 1, function(r) { o <- order(r, decreasing = TRUE); if (r[o[1]] - r[o[2]] >= margin) colnames(cl_means)[o[1]] else "uncertain" })
  cell_best <- colnames(sc)[max.col(replace(sc, is.na(sc), -Inf), ties.method = "first")]
  data.table(barcode = colnames(lx), cluster = cl, lineage_cluster = call[as.character(cl)], lineage_cell_best = cell_best, sc)
}

# ---- CNV scoring (expression-based; inferCNV/Tirosh-style reference bounds) ------------
# lx: log-normalised genes x cells for one sample; ref_group: character vector (NA for
# non-reference cells) naming the reference lineage of each reference cell; gene_pos:
# data.table(symbol, chr, start). Expression is centred on the mean of reference-group
# means; deviations inside the range of reference-group means are set to zero before
# smoothing, so lineage-specific expression in references does not register as CNV.
cnv_scores <- function(lx, ref_group, gene_pos, window = 101, n_genes = 7000, min_ref_group = 10) {
  ref <- !is.na(ref_group)
  keep_groups <- names(which(table(ref_group[ref]) >= min_ref_group))
  ref <- ref & ref_group %in% keep_groups
  gp <- gene_pos[symbol %in% rownames(lx) & chr %in% as.character(1:22)]
  gp <- gp[!duplicated(symbol)]
  mu <- Matrix::rowMeans(lx[gp$symbol, , drop = FALSE])
  gp <- gp[order(-mu)][seq_len(min(n_genes, nrow(gp)))]
  gp[, chr_n := as.integer(chr)]; setorder(gp, chr_n, start)
  X <- as.matrix(lx[gp$symbol, , drop = FALSE])
  gmeans <- sapply(keep_groups, function(g) rowMeans(X[, ref & ref_group == g, drop = FALSE]))
  if (is.null(dim(gmeans))) gmeans <- matrix(gmeans, ncol = 1)
  centre <- rowMeans(gmeans); lo <- apply(gmeans, 1, min) - centre; hi <- apply(gmeans, 1, max) - centre
  X <- X - centre
  X <- ifelse(X > hi, X - hi, ifelse(X < lo, X - lo, 0))
  X[X > 3] <- 3; X[X < -3] <- -3
  sm <- do.call(rbind, lapply(split(seq_len(nrow(X)), gp$chr_n), function(ix) {
    if (length(ix) < window) return(NULL)
    apply(X[ix, , drop = FALSE], 2, function(col) as.numeric(stats::filter(col, rep(1 / window, window), sides = 2)))
  }))
  sm <- sm[stats::complete.cases(sm), , drop = FALSE]
  signal <- colMeans(sm^2)
  thr <- stats::quantile(signal[ref], 0.95)
  nonref <- which(!ref)
  top <- nonref[signal[nonref] > stats::quantile(signal[nonref], 0.75)]
  prof <- if (length(top) >= 10) rowMeans(sm[, top, drop = FALSE]) else rowMeans(sm[, nonref, drop = FALSE])
  corr <- suppressWarnings(as.numeric(stats::cor(sm, prof)))
  data.table(barcode = colnames(lx), cnv_signal = signal, cnv_corr = corr, ref_signal_q95 = as.numeric(thr),
             n_genes_cnv = nrow(sm), n_ref_cells = sum(ref), ref_groups = paste(keep_groups, collapse = ","))
}

pseudobulk <- function(m, groups) {
  f <- factor(groups)
  agg <- Matrix::sparse.model.matrix(~ 0 + f)
  r <- m %*% agg
  colnames(r) <- levels(f)
  r
}
