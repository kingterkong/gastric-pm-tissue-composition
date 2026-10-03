# Stage 4: prepare bulk expression inputs (no site comparisons).
# - Ensembl GRCh37 r75 gene/transcript annotation and HGNC mapping (kept with original IDs)
# - GSE314812: RSEM expected counts -> DGEList-ready matrix for the audited specimens
# - GSE237876: Salmon transcript quant.sf -> gene-level counts (sum of NumReads, tximport
#   countsFromAbundance = "no" equivalent) and gene TPM (sum of transcript TPM)
# - matrix-type audit table for both sources (dimensions, integer fraction, column sums)
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
suppressPackageStartupMessages({ library(Matrix) })

gtf_f <- pm_raw("other_sources", "annotation", "Homo_sapiens.GRCh37.75.gtf.gz")
hgnc_f <- pm_raw("other_sources", "annotation", "hgnc_complete_set_20260925.txt")
c314_f <- pm_raw("GSE314812", "geo", "GSE314812_rsem.merged.gene_counts.tsv.gz")
tar237 <- pm_raw("GSE237876", "GSE237876_RAW.tar")
cw_f <- pm_in("00_admin", "sample_crosswalk.csv")
step_begin("04_prepare_bulk", inputs = c(gtf_f, hgnc_f, c314_f, tar237, cw_f))

# ---- annotation ----------------------------------------------------------------
gtf <- fread(cmd = sprintf("gzip -dc %s | grep -v '^#'", shQuote(gtf_f)), sep = "\t", header = FALSE,
             col.names = c("chr", "src", "feature", "start", "end", "score", "strand", "frame", "attr"))
attr_get <- function(a, key) sub(sprintf('.*%s "([^"]+)".*', key), "\\1", a)
genes <- gtf[feature == "gene"]
genes[, `:=`(gene_id = attr_get(attr, "gene_id"), gene_name = attr_get(attr, "gene_name"), gene_biotype = attr_get(attr, "gene_biotype"))]
ex <- gtf[feature == "exon"]
ex[, `:=`(gene_id = attr_get(attr, "gene_id"), transcript_id = attr_get(attr, "transcript_id"))]
# union exon length per gene
setkey(ex, gene_id, start)
ulen <- ex[, {
  s <- start; e <- end; o <- order(s); s <- s[o]; e <- e[o]
  tot <- 0; cs <- s[1]; ce <- e[1]
  if (length(s) > 1) for (i in 2:length(s)) { if (s[i] <= ce + 1) ce <- max(ce, e[i]) else { tot <- tot + ce - cs + 1; cs <- s[i]; ce <- e[i] } }
  .(union_exon_length = tot + ce - cs + 1)
}, by = gene_id]
tx2gene <- unique(ex[, .(transcript_id, gene_id)])
ann <- merge(genes[, .(gene_id, gene_name, gene_biotype, chr, start, end, strand)], ulen, by = "gene_id", all.x = TRUE)

hgnc <- fread(hgnc_f, sep = "\t", quote = "", colClasses = "character", select = c("hgnc_id", "symbol", "status", "prev_symbol", "alias_symbol", "ensembl_gene_id", "locus_group"))
hgnc <- hgnc[status == "Approved"]
ann <- merge(ann, hgnc[nzchar(ensembl_gene_id), .(ensembl_gene_id, hgnc_symbol = symbol, hgnc_id, locus_group)],
             by.x = "gene_id", by.y = "ensembl_gene_id", all.x = TRUE)
ann[, symbol_for_scoring := fifelse(!is.na(hgnc_symbol), hgnc_symbol, gene_name)]
ann[, symbol_source := fifelse(!is.na(hgnc_symbol), "HGNC approved (via Ensembl ID)", "Ensembl r75 gene_name")]
o_ann <- pm_out("02_data_processed", "bulk", "annotation_grch37_r75_hgnc.csv"); write_csv_atomic(ann, o_ann)
o_tx <- pm_out("02_data_processed", "bulk", "tx2gene_grch37_r75.csv"); write_csv_atomic(tx2gene, o_tx)

cw <- fread(cw_f, colClasses = "character")

audit_matrix <- function(m, label, unit_note) {
  v <- as.numeric(m)
  data.table(source = label, n_features = nrow(m), n_samples = ncol(m), n_na = sum(is.na(v)), n_negative = sum(v < 0, na.rm = TRUE),
             frac_zero = mean(v == 0), frac_noninteger = mean(abs(v - round(v)) > 1e-8), max_value = max(v),
             colsum_min = min(colSums(m)), colsum_median = median(colSums(m)), colsum_max = max(colSums(m)),
             duplicated_feature_ids = sum(duplicated(rownames(m))), duplicated_columns = sum(duplicated(colnames(m))),
             unit_judgement = unit_note)
}

# ---- GSE314812 ------------------------------------------------------------------
x <- fread(c314_f)
stopifnot(identical(names(x)[1:2], c("gene_id", "transcript_id(s)")))
m314 <- as.matrix(x[, -(1:2)]); rownames(m314) <- x$gene_id
x314 <- cw[dataset_id == "GSE314812"]
stopifnot(setequal(colnames(m314), x314$expression_column))
m314 <- m314[, x314$expression_column]
aud <- list(audit_matrix(m314, "GSE314812", "RSEM expected counts (gene-level; non-integer from multi-mapping allocation; column sums ~2-3e7 not 1e6)"))
missing_ann <- setdiff(rownames(m314), ann$gene_id)
save_rds_atomic(list(counts = m314, genes = ann[match(rownames(m314), gene_id)], samples = x314,
                     note = "values are RSEM expected counts; no effective lengths published"),
                pm_out("02_data_processed", "bulk", "GSE314812_counts.rds"))

# ---- GSE237876 ------------------------------------------------------------------
tmpd <- tempfile("gse237876_"); dir.create(tmpd)
on.exit(unlink(tmpd, recursive = TRUE), add = TRUE)
system2("tar", c("-xf", shQuote(tar237), "-C", shQuote(tmpd)))
x237 <- cw[dataset_id == "GSE237876"]
q <- lapply(x237$expression_column, function(f) {
  d <- fread(file.path(tmpd, f), select = c("Name", "Length", "EffectiveLength", "TPM", "NumReads"))
  d
})
tx_ids <- q[[1]]$Name
stopifnot(all(vapply(q, function(d) identical(d$Name, tx_ids), TRUE)))
tx_map <- tx2gene[match(tx_ids, transcript_id)]
unmapped_tx <- sum(is.na(tx_map$gene_id))
g_id <- tx_map$gene_id
keep <- !is.na(g_id)
counts_tx <- vapply(q, function(d) d$NumReads, numeric(length(tx_ids)))
tpm_tx <- vapply(q, function(d) d$TPM, numeric(length(tx_ids)))
colnames(counts_tx) <- colnames(tpm_tx) <- x237$original_sample_id
agg <- function(mat) { r <- rowsum(mat[keep, , drop = FALSE], g_id[keep]); r }
counts_g <- agg(counts_tx); tpm_g <- agg(tpm_tx)
aud[[2]] <- audit_matrix(counts_tx, "GSE237876 (transcript NumReads)", "Salmon estimated counts (transcript-level)")
aud[[3]] <- audit_matrix(counts_g, "GSE237876 (gene counts, summed)", "gene-level estimated counts (sum of transcript NumReads)")
save_rds_atomic(list(counts = counts_g, tpm = tpm_g, genes = ann[match(rownames(counts_g), gene_id)], samples = x237,
                     tx_mapping = data.table(n_transcripts = length(tx_ids), n_unmapped_to_r75 = unmapped_tx,
                                             unmapped_tpm_share_median = median(colSums(tpm_tx[!keep, , drop = FALSE]) / 1e6))),
                pm_out("02_data_processed", "bulk", "GSE237876_gene_counts.rds"))

aud <- rbindlist(aud, fill = TRUE)
aud[, annotation_unmapped_genes := c(length(missing_ann), unmapped_tx, NA)]
o_aud <- pm_out("04_results", "audit", "matrix_type_audit.csv"); write_csv_atomic(aud, o_aud)
print(aud)
step_end(outputs = c(o_ann, o_tx, pm_out("02_data_processed", "bulk", "GSE314812_counts.rds"),
                     pm_out("02_data_processed", "bulk", "GSE237876_gene_counts.rds"), o_aud))
