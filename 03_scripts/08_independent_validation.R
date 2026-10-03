# Stage 8: independent-source checks of the frozen primary endpoint. Runs the identical
# frozen estimand once in GSE237876 (independent paired bulk, 5 patients) and reports the
# separately-defined fluid-context estimand in GSE239676. Nothing is tuned here.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
source(pm_src("03_scripts", "R", "score_utils.R"))
cfg <- pm_config()

prog_f <- pm_in("02_data_processed", "frozen_programs", "program_definitions.rds")
lock_f <- pm_in("00_admin", "protocol_lock.json")
b237_f <- pm_in("02_data_processed", "bulk", "GSE237876_gene_counts.rds")
step_begin("08_independent_validation", inputs = c(prog_f, lock_f, b237_f), seed = cfg$project$seed + 8, uses_protocol = TRUE)

progs <- readRDS(prog_f)
gene_sets <- lapply(progs, `[[`, "genes")
outs <- character()

# ---- GSE237876: independent paired bulk (frozen estimand, one run) -------------------
b <- readRDS(b237_f)
smp <- as.data.table(b$samples)
sm <- score_matrix(b$counts, b$genes$symbol_for_scoring, gene_sets)
sc <- sm$scores; setnames(sc, "sample", "original_sample_id")
sc <- merge(smp[, .(original_sample_id, canonical_patient_id, lesion_role, region_id, therapy_before_sampling,
                    collection_time_or_order, purity_published, eligible_validation)], sc, by = "original_sample_id")
o1 <- pm_out("04_results", "validation", "GSE237876_program_scores_by_specimen.csv"); write_csv_atomic(sc, o1)
o1b <- pm_out("04_results", "validation", "GSE237876_program_coverage.csv")
write_csv_atomic(cbind(sm$coverage, universe_size = sm$universe_size), o1b); outs <- c(outs, o1, o1b)

prog_ids <- names(gene_sets)
val <- sc[eligible_validation == TRUE & lesion_role %in% c("primary_tumor", "solid_peritoneal_metastasis")]
# pre-specified multi-region rule: average eligible regions within patient x site first
agg <- val[, lapply(.SD, mean), by = .(canonical_patient_id, lesion_role), .SDcols = prog_ids]
nreg <- val[, .(n_specimens = .N), by = .(canonical_patient_id, lesion_role)]
agg <- merge(agg, nreg, by = c("canonical_patient_id", "lesion_role"))
o2 <- pm_out("04_results", "validation", "GSE237876_patient_site_values.csv"); write_csv_atomic(agg, o2); outs <- c(outs, o2)

w <- dcast(agg, canonical_patient_id ~ lesion_role, value.var = prog_ids)
D <- data.table(canonical_patient_id = w$canonical_patient_id)
for (p in prog_ids) D[, (p) := w[[paste0(p, "_solid_peritoneal_metastasis")]] - w[[paste0(p, "_primary_tumor")]]]
o3 <- pm_out("04_results", "validation", "GSE237876_patient_deltas.csv"); write_csv_atomic(D, o3); outs <- c(outs, o3)

ev <- sm$coverage[evaluable == TRUE, program_id]
val_res <- rbindlist(lapply(prog_ids, function(p) {
  if (!p %in% ev) return(data.table(program_id = p, cohort = "GSE237876", n = NA_integer_,
                                    note = "not evaluable (coverage rule)"))
  s <- paired_summary(D[[p]], n_boot = 5000, seed = cfg$project$seed)
  cbind(data.table(program_id = p, cohort = "GSE237876"), s)
}), fill = TRUE)
# identity-conflict sensitivity: GCM05's resected metastasis is labelled peritoneal in the
# sample table but liver in both clinical-course figures
D4 <- D[canonical_patient_id != "GCM05"]
val_res_excl <- rbindlist(lapply(intersect(prog_ids, ev), function(p) {
  cbind(data.table(program_id = p, cohort = "GSE237876 excluding GCM05"), paired_summary(D4[[p]], n_boot = 5000, seed = cfg$project$seed))
}), fill = TRUE)
val_all <- rbind(val_res, val_res_excl, fill = TRUE)
val_all[, max_attainable_two_sided_p := ifelse(!is.na(n) & n > 0, 2^(-(n - 1)), NA_real_)]
val_all[, precision_note := "n=5 (or 4): the smallest attainable two-sided exact P exceeds 0.05; directional consistency only"]
o4 <- pm_out("04_results", "validation", "GSE237876_frozen_endpoint.csv"); write_csv_atomic(val_all, o4); outs <- c(outs, o4)

# ---- GSE239676: fluid-context estimand, reported separately -------------------------
pb_f <- pm_in("02_data_processed", "singlecell", "GSE239676_pseudobulk.rds")
if (file.exists(pb_f)) {
  pb <- readRDS(pb_f)
  meta <- as.data.table(pb$meta)
  MIN_CELLS <- 20
  keep <- meta$cell_class == "malignant_cnv_high" & meta$n_cells >= MIN_CELLS &
    meta$lesion_role %in% c("primary_tumor", "malignant_ascites")
  if (sum(keep) >= 4) {
    cnts <- as.matrix(pb$counts[, keep, drop = FALSE])
    mk <- meta[keep]
    smf <- score_matrix(cnts, rownames(cnts), gene_sets)
    s2 <- cbind(mk[, .(sample_id, canonical_patient_id, lesion_role, n_cells, total_umi)], smf$scores[, -1])
    o5 <- pm_out("04_results", "validation", "GSE239676_malignant_pseudobulk_scores.csv"); write_csv_atomic(s2, o5)
    o5b <- pm_out("04_results", "validation", "GSE239676_program_coverage.csv")
    write_csv_atomic(cbind(smf$coverage, universe_size = smf$universe_size), o5b); outs <- c(outs, o5, o5b)
    both <- s2[, .N, by = .(canonical_patient_id, lesion_role)][, .N, by = canonical_patient_id][N == 2, canonical_patient_id]
    if (length(both) >= 3) {
      w2 <- dcast(s2[canonical_patient_id %in% both], canonical_patient_id ~ lesion_role, value.var = names(gene_sets))
      evf <- smf$coverage[evaluable == TRUE, program_id]
      fl <- rbindlist(lapply(intersect(prog_ids, evf), function(p) {
        d <- w2[[paste0(p, "_malignant_ascites")]] - w2[[paste0(p, "_primary_tumor")]]
        cbind(data.table(program_id = p, cohort = "GSE239676 malignant cells (ascites - primary)"),
              paired_summary(d, n_boot = 5000, seed = cfg$project$seed))
      }), fill = TRUE)
      fl[, estimand_note := "SEPARATE estimand: free-floating malignant cells in ascites versus primary tumour; NOT a solid peritoneal metastasis replication"]
      o6 <- pm_out("04_results", "validation", "GSE239676_fluid_context_estimand.csv"); write_csv_atomic(fl, o6); outs <- c(outs, o6)
      print(fl[program_id == "HRC_CORE"])
    }
  }
}

cat("\n=== GSE237876 frozen endpoint (independent paired bulk) ===\n")
print(val_all[program_id %in% c("HRC_CORE", "CTRL_EMT", "CTRL_STROMA", "NEG_YAP_STROMAL"),
              .(program_id, cohort, n, mean = round(mean, 4), ci_low = round(ci_low, 4), ci_high = round(ci_high, 4),
                n_positive, n_negative, p_value = signif(p_value, 3))])
step_end(outputs = outs)
