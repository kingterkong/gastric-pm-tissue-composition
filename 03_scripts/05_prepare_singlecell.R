# Stage 5: single-cell preparation (program-agnostic; no site comparisons).
# For each source: read counts, harmonise symbols to HGNC, per-sample robust QC, per-sample
# marker-based lineage annotation (fixed lineage markers only), expression-based CNV
# scoring of epithelial cells against same-sample non-epithelial reference cells, and
# pseudobulk counts per sample x cell class. Outputs go to 02_data_processed/singlecell.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
source(pm_src("03_scripts", "R", "sc_utils.R"))
cfg <- pm_config()

hgnc_f <- pm_raw("other_sources", "annotation", "hgnc_complete_set_20260925.txt")
ann_f <- pm_in("02_data_processed", "bulk", "annotation_grch37_r75_hgnc.csv")
cw_f <- pm_in("00_admin", "sample_crosswalk.csv")
src <- list(
  GSE183904 = pm_raw("other_sources", "GSE183904", "GSE183904_RAW.tar"),
  GSE163558 = pm_raw("other_sources", "GSE163558", "GSE163558_RAW.tar"),
  GSE308231 = pm_raw("other_sources", "GSE308231", "GSE308231_RAW.tar"),
  GSE239676 = c(pm_raw("other_sources", "GSE239676", "GSE239676_count_matrix.mtx.gz"),
                pm_raw("other_sources", "GSE239676", "GSE239676_barcodes.tsv.gz"),
                pm_raw("other_sources", "GSE239676", "GSE239676_features.tsv.gz"),
                pm_raw("other_sources", "GSE239676", "GSE239676_meta.tsv.gz"))
)
step_begin("05_prepare_singlecell", inputs = c(hgnc_f, ann_f, cw_f, unlist(src)), seed = cfg$project$seed + 5)

maps <- hgnc_maps(hgnc_f)
ann <- fread(ann_f)
gene_pos <- unique(ann[!is.na(symbol_for_scoring), .(symbol = symbol_for_scoring, chr, start)])
cw <- fread(cw_f, colClasses = "character")
outdir <- pm_out("02_data_processed", "singlecell", "x"); outdir <- dirname(outdir)

REF_LINEAGES <- c("fibroblast", "mural", "endothelial", "myeloid", "t_nk", "b_cell", "plasma", "mast")

process_sample <- function(m, sym, ens, sample_id, dataset_id, seed) {
  hm <- harmonise_symbols(sym, ens, maps)
  m <- collapse_rows(m, hm$symbol)
  mt <- grep("^MT-", rownames(m), value = TRUE)
  qc <- qc_cells(m, sample_id, mt)
  m <- m[, qc$qc_pass, drop = FALSE]
  if (ncol(m) < 30) return(list(qc = qc, cells = NULL, counts = NULL, map = hm))
  # doublets (scDblFinder, per sample); doublets are excluded and counted in the QC table
  set.seed(seed)
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = m))
  sce <- suppressWarnings(suppressMessages(scDblFinder::scDblFinder(sce, verbose = FALSE)))
  dbl <- sce$scDblFinder.class == "doublet"
  qc[qc_pass == TRUE, doublet := dbl]
  qc[is.na(doublet), doublet := FALSE]
  rm(sce)
  m <- m[, !dbl, drop = FALSE]
  lx <- lognorm(m)
  la <- annotate_lineages(m, lx, seed = seed)
  la[, lineage := lineage_cluster]
  la[, lineage_high_conf := lineage_cluster == lineage_cell_best]
  epi <- la$lineage == "epithelial"
  ref_group <- ifelse(la$lineage %in% REF_LINEAGES, la$lineage, NA_character_)
  cnv <- if (sum(epi) >= 10 && sum(!is.na(ref_group)) >= 30) cnv_scores(lx, ref_group, gene_pos) else NULL
  cells <- cbind(data.table(dataset_id = dataset_id, sample_id = sample_id), la)
  if (!is.null(cnv)) cells <- merge(cells, cnv, by = "barcode", all.x = TRUE, sort = FALSE) else
    cells[, `:=`(cnv_signal = NA_real_, cnv_corr = NA_real_, ref_signal_q95 = NA_real_, n_genes_cnv = NA_integer_, n_ref_cells = NA_integer_, ref_groups = NA_character_)]
  cells <- merge(cells, qc[qc_pass == TRUE & doublet == FALSE, .(barcode, n_counts, n_genes, pct_mt)], by = "barcode", sort = FALSE)
  list(qc = qc, cells = cells, counts = m[, cells$barcode, drop = FALSE], map = hm)
}

finalise_dataset <- function(res, dataset_id) {
  res <- Filter(function(r) !is.null(r$cells), res)
  cells <- rbindlist(lapply(res, `[[`, "cells"), fill = TRUE)
  qc <- rbindlist(lapply(res, `[[`, "qc"))
  genes <- Reduce(union, lapply(res, function(r) rownames(r$counts)))
  pad <- function(m) { miss <- setdiff(genes, rownames(m)); if (length(miss)) m <- rbind(m, Matrix(0, length(miss), ncol(m), sparse = TRUE, dimnames = list(miss, colnames(m)))); m[genes, , drop = FALSE] }
  counts <- do.call(cbind, lapply(res, function(r) { m <- pad(r$counts); colnames(m) <- paste0(r$cells$sample_id[1], "|", colnames(m)); m }))
  cells[, cell_id := paste0(sample_id, "|", barcode)]
  stopifnot(identical(colnames(counts), cells$cell_id))
  gene_measured <- rbindlist(lapply(res, function(r) data.table(sample_id = r$cells$sample_id[1], symbol = rownames(r$counts))))
  maps_used <- unique(rbindlist(lapply(res, `[[`, "map")))
  list(cells = cells, qc = qc, counts = counts, gene_measured = gene_measured, symbol_map = maps_used)
}

# malignant classification (frozen rule; independent of any program genes)
classify_epithelium <- function(cells, cw_ds) {
  cells <- merge(cells, cw_ds[, .(sample_id, canonical_patient_id, lesion_role, specimen_context)], by = "sample_id", all.x = TRUE, sort = FALSE)
  tumour_like <- cells$lesion_role %in% c("primary_tumor", "solid_peritoneal_metastasis", "malignant_ascites", "liver_metastasis",
                                          "ovary_metastasis", "nodal_metastasis", "peritoneal_washings")
  cells[, cell_class := lineage]
  cells[lineage == "epithelial" & !tumour_like, cell_class := "epithelial_normal_sample"]
  cells[lineage == "epithelial" & tumour_like & !is.na(cnv_signal) & cnv_signal > ref_signal_q95 & cnv_corr > 0.2, cell_class := "malignant_cnv_high"]
  cells[lineage == "epithelial" & tumour_like & cell_class == "epithelial", cell_class := "tumour_epithelial_cnv_low_or_na"]
  cells
}

save_ds <- function(ds, dataset_id) {
  o <- list()
  o$counts <- pm_out("02_data_processed", "singlecell", paste0(dataset_id, "_counts_qc.rds")); save_rds_atomic(ds$counts, o$counts)
  o$cells <- pm_out("02_data_processed", "singlecell", paste0(dataset_id, "_cells.csv.gz")); write_csv_atomic(ds$cells, o$cells)
  o$qc <- pm_out("02_data_processed", "singlecell", paste0(dataset_id, "_qc_all_barcodes.csv.gz")); write_csv_atomic(ds$qc, o$qc)
  o$gm <- pm_out("02_data_processed", "singlecell", paste0(dataset_id, "_genes_measured.csv.gz")); write_csv_atomic(ds$gene_measured, o$gm)
  o$map <- pm_out("02_data_processed", "singlecell", paste0(dataset_id, "_symbol_map.csv.gz")); write_csv_atomic(ds$symbol_map, o$map)
  # pseudobulk per sample x cell_class (malignant, lineages)
  grp <- paste(ds$cells$sample_id, ds$cells$cell_class, sep = "||")
  pb <- pseudobulk(ds$counts, grp)
  pb_meta <- ds$cells[, .(n_cells = .N, total_umi = sum(n_counts)), by = .(sample_id, cell_class, canonical_patient_id, lesion_role, specimen_context)]
  pb_meta[, group := paste(sample_id, cell_class, sep = "||")]
  pb_meta <- pb_meta[match(colnames(pb), group)]
  o$pb <- pm_out("02_data_processed", "singlecell", paste0(dataset_id, "_pseudobulk.rds")); save_rds_atomic(list(counts = pb, meta = pb_meta), o$pb)
  unlist(o)
}

outs <- character()
tmp_root <- tempfile("sc_"); dir.create(tmp_root)
on.exit(unlink(tmp_root, recursive = TRUE), add = TRUE)
seed0 <- cfg$project$seed

# ---- GSE183904 -------------------------------------------------------------------
d <- file.path(tmp_root, "GSE183904"); dir.create(d)
system2("tar", c("-xf", shQuote(src$GSE183904), "-C", shQuote(d)))
cw183 <- cw[dataset_id == "GSE183904"][, sample_id := original_sample_id]
files <- list.files(d, pattern = "\\.csv\\.gz$", full.names = TRUE)
res <- lapply(seq_along(files), function(i) {
  f <- files[i]; sid <- sub("^GSM[0-9]+_(sample[0-9]+)\\.csv\\.gz$", "\\1", basename(f))
  r <- read_dense_csv_gz(f)
  message("GSE183904 ", sid, ": ", ncol(r$counts), " barcodes")
  out <- process_sample(r$counts, r$symbol, NULL, sid, "GSE183904", seed0 + i)
  rm(r); gc(verbose = FALSE); out
})
ds <- finalise_dataset(res, "GSE183904"); rm(res); gc(verbose = FALSE)
ds$cells <- classify_epithelium(ds$cells, cw183)
outs <- c(outs, save_ds(ds, "GSE183904")); rm(ds); gc(verbose = FALSE)
unlink(d, recursive = TRUE)

# ---- 10x triplet sources (GSE163558, GSE308231) ------------------------------------
triplet_source <- function(tarf, dataset_id, sample_of) {
  d <- file.path(tmp_root, dataset_id); dir.create(d)
  system2("tar", c("-xf", shQuote(tarf), "-C", shQuote(d)))
  mtx <- list.files(d, pattern = "matrix\\.mtx\\.gz$", full.names = TRUE)
  res <- lapply(seq_along(mtx), function(i) {
    pre <- sub("matrix\\.mtx\\.gz$", "", mtx[i])
    r <- read_10x_triplet(mtx[i], paste0(pre, "features.tsv.gz"), paste0(pre, "barcodes.tsv.gz"))
    colnames(r$counts) <- r$barcodes
    sid <- sample_of(basename(pre))
    message(dataset_id, " ", sid, ": ", ncol(r$counts), " barcodes")
    process_sample(r$counts, r$symbol, r$ensembl, sid, dataset_id, seed0 + 100 + i)
  })
  unlink(d, recursive = TRUE)
  finalise_dataset(res, dataset_id)
}
cw163 <- cw[dataset_id == "GSE163558"][, sample_id := gsm_or_run_id]
ds <- triplet_source(src$GSE163558, "GSE163558", function(pre) sub("^(GSM[0-9]+)_.*$", "\\1", pre))
ds$cells <- classify_epithelium(ds$cells, cw163)
outs <- c(outs, save_ds(ds, "GSE163558")); rm(ds); gc(verbose = FALSE)

cw308 <- cw[dataset_id == "GSE308231"][, sample_id := gsm_or_run_id]
ds <- triplet_source(src$GSE308231, "GSE308231", function(pre) sub("^(GSM[0-9]+)_.*$", "\\1", pre))
ds$cells <- classify_epithelium(ds$cells, cw308)
outs <- c(outs, save_ds(ds, "GSE308231")); rm(ds); gc(verbose = FALSE)

# ---- GSE239676: stream-filter to primary + ascites of patients with both --------
cw239 <- cw[dataset_id == "GSE239676"]
keep_samples <- cw239[grepl("_primary_ascites$", primary_pm_pair_id), expression_column]
bc <- fread(src$GSE239676[2], header = FALSE)$V1
mt239 <- fread(src$GSE239676[4])
stopifnot(length(bc) == nrow(mt239))
keep_idx <- which(mt239$Sample %in% keep_samples)
idx_file <- file.path(tmp_root, "keep_cols.txt"); writeLines(as.character(keep_idx), idx_file)
flt <- file.path(tmp_root, "gse239676_filtered.triplets")
awk <- sprintf("gzip -dc %s | awk 'NR==FNR{k[$1]=1; next} /^%%/{next} !h{h=1; next} ($2 in k){print $1, $2, $3}' %s - > %s",
               shQuote(src$GSE239676[1]), shQuote(idx_file), shQuote(flt))
st <- system(awk)
if (st != 0) stop("GSE239676 stream filter failed")
tr <- fread(flt, header = FALSE, col.names = c("i", "j", "x"))
feat <- fread(src$GSE239676[3], header = FALSE)$V1
newj <- match(tr$j, keep_idx)
M <- sparseMatrix(i = tr$i, j = newj, x = tr$x, dims = c(length(feat), length(keep_idx)))
rm(tr); gc(verbose = FALSE)
colnames(M) <- bc[keep_idx]
smp <- mt239$Sample[keep_idx]
res <- lapply(seq_along(keep_samples), function(i) {
  s <- keep_samples[i]
  message("GSE239676 ", s, ": ", sum(smp == s), " cells")
  process_sample(M[, smp == s, drop = FALSE], feat, NULL, s, "GSE239676", seed0 + 200 + i)
})
rm(M); gc(verbose = FALSE)
ds <- finalise_dataset(res, "GSE239676"); rm(res)
cw239[, sample_id := expression_column]
ds$cells <- classify_epithelium(ds$cells, cw239)
outs <- c(outs, save_ds(ds, "GSE239676")); rm(ds); gc(verbose = FALSE)

# ---- summary --------------------------------------------------------------------
summ <- rbindlist(lapply(c("GSE183904", "GSE163558", "GSE308231", "GSE239676"), function(g) {
  cl <- fread(pm_in("02_data_processed", "singlecell", paste0(g, "_cells.csv.gz")))
  cl[, .N, by = .(dataset_id, sample_id, lesion_role, cell_class)]
}))
o_s <- pm_out("04_results", "audit", "singlecell_cell_class_counts.csv"); write_csv_atomic(summ, o_s)
step_end(outputs = c(outs, o_s))
