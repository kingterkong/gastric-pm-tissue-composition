# Stage 9: cell source and composition. Reports the three levels separately (plan 16.5):
# (1) in which cell types the frozen program is measurable, (2) how cell-type composition
# differs, (3) whether the program value changes inside a fixed cell type. Coupling between
# two programs is NOT a claim of this study (see protocol section 1); only the pre-specified
# host-compartment relationship is examined.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
source(pm_src("03_scripts", "R", "score_utils.R"))
cfg <- pm_config()

prog_f <- pm_in("02_data_processed", "frozen_programs", "program_definitions.rds")
lock_f <- pm_in("00_admin", "protocol_lock.json")
step_begin("09_composition_and_coupling", inputs = c(prog_f, lock_f), seed = cfg$project$seed + 9, uses_protocol = TRUE)
progs <- readRDS(prog_f); gene_sets <- lapply(progs, `[[`, "genes")
outs <- character()
MIN_CELLS <- 20

ds_all <- c("GSE183904", "GSE163558", "GSE308231", "GSE239676")
ds <- ds_all[file.exists(vapply(ds_all, function(g) pm_in("02_data_processed", "singlecell", paste0(g, "_pseudobulk.rds")), ""))]

# ---- (1) where is the program measurable, and in which lineages is it expressed ------
# Calibration sources only (development side): GSE183904 gastric donors + GSE163558.
calib <- intersect(ds, c("GSE183904", "GSE163558"))
lineage_rows <- list(); comp_rows <- list(); within_rows <- list(); score_rows <- list()
for (g in ds) {
  pb <- readRDS(pm_in("02_data_processed", "singlecell", paste0(g, "_pseudobulk.rds")))
  meta <- as.data.table(pb$meta); cnts <- pb$counts
  keep <- meta$n_cells >= MIN_CELLS
  if (!any(keep)) next
  sm <- score_matrix(as.matrix(cnts[, keep, drop = FALSE]), rownames(cnts), gene_sets)
  sc <- cbind(meta[keep], sm$scores[, -1])
  sc[, dataset_id := g]
  score_rows[[g]] <- sc
  lineage_rows[[g]] <- cbind(data.table(dataset_id = g), sm$coverage)
  # (2) composition: fraction of QC-passing cells per class within each sample
  cells <- fread(pm_in("02_data_processed", "singlecell", paste0(g, "_cells.csv.gz")))
  comp <- cells[, .N, by = .(dataset_id, sample_id, canonical_patient_id, lesion_role, specimen_context, cell_class)]
  comp[, frac := N / sum(N), by = .(sample_id)]
  comp_rows[[g]] <- comp
}
scores_all <- rbindlist(score_rows, fill = TRUE)
o1 <- pm_out("04_results", "cell_source", "pseudobulk_program_scores.csv"); write_csv_atomic(scores_all, o1)
o2 <- pm_out("04_results", "cell_source", "pseudobulk_program_coverage.csv"); write_csv_atomic(rbindlist(lineage_rows, fill = TRUE), o2)
comp_all <- rbindlist(comp_rows, fill = TRUE)
o3 <- pm_out("04_results", "cell_source", "cell_class_composition.csv"); write_csv_atomic(comp_all, o3)
outs <- c(outs, o1, o2, o3)

# Lineage profile of the primary program in calibration data only (development side).
if (length(calib)) {
  lin <- scores_all[dataset_id %in% calib & lesion_role %in% c("primary_tumor", "normal_adjacent_gastric")]
  prof <- lin[, .(n_groups = .N, n_cells = sum(n_cells), median_HRC = median(HRC_CORE, na.rm = TRUE),
                  median_CTRL_STROMA = median(CTRL_STROMA, na.rm = TRUE)), by = .(dataset_id, cell_class)][order(dataset_id, -median_HRC)]
  o4 <- pm_out("04_results", "cell_source", "calibration_lineage_profile.csv"); write_csv_atomic(prof, o4); outs <- c(outs, o4)
  cat("\n=== primary program by cell class (calibration sources; development side) ===\n"); print(prof)
}

# ---- (3) within fixed cell type: malignant cells in solid PM vs primary libraries ----
# Descriptive only. GSE308231 donor identity is unknown, so libraries are the unit and no
# paired inference is possible; GSE183904 peritoneal donors are held out from calibration.
mal <- scores_all[cell_class == "malignant_cnv_high" &
                    lesion_role %in% c("primary_tumor", "solid_peritoneal_metastasis")]
if (nrow(mal)) {
  o5 <- pm_out("04_results", "cell_source", "malignant_pseudobulk_solid_sites.csv"); write_csv_atomic(mal, o5); outs <- c(outs, o5)
  summ <- mal[, .(n_libraries = .N, n_cells = sum(n_cells), median_HRC = median(HRC_CORE, na.rm = TRUE),
                  min_HRC = min(HRC_CORE, na.rm = TRUE), max_HRC = max(HRC_CORE, na.rm = TRUE)),
              by = .(dataset_id, lesion_role)]
  summ[, unit_note := "library-level description; donors unverified in GSE308231; not a patient-level test"]
  o6 <- pm_out("04_results", "cell_source", "malignant_library_summary.csv"); write_csv_atomic(summ, o6); outs <- c(outs, o6)
  cat("\n=== malignant-cell pseudobulk in solid specimens (descriptive) ===\n"); print(summ)
}

# ---- host-compartment admixture check in the bulk target cohort ----------------------
D <- fread(pm_in("04_results", "primary", "GSE314812_patient_deltas.csv"))
comp_cols <- grep("^COMPART_", names(D), value = TRUE)
not_evaluable <- comp_cols[vapply(comp_cols, function(k) all(is.na(D[[k]])), TRUE)]
comp_cols <- setdiff(comp_cols, not_evaluable)
if (length(not_evaluable)) message("Compartment sets not evaluable under the frozen coverage rule: ", paste(not_evaluable, collapse = ", "))
adm <- rbindlist(lapply(comp_cols, function(k) {
  s <- paired_summary(D[[k]], n_boot = 5000, seed = cfg$project$seed)
  cbind(data.table(compartment = k), s)
}), fill = TRUE)
adm[, q_value := p.adjust(p_value, method = "BH")]
if (length(not_evaluable)) adm <- rbind(adm, data.table(compartment = not_evaluable, n = NA_integer_,
  note = "not evaluable on this platform under the frozen coverage rule"), fill = TRUE)
o7 <- pm_out("04_results", "cell_source", "host_compartment_deltas.csv"); write_csv_atomic(adm, o7); outs <- c(outs, o7)

# relationship between the primary program delta and the compartment deltas
rel <- rbindlist(lapply(comp_cols, function(k) {
  ct <- suppressWarnings(stats::cor.test(D$HRC_CORE, D[[k]], method = "spearman", exact = FALSE))
  data.table(compartment = k, spearman_rho = unname(ct$estimate), p_value = ct$p.value, n = sum(!is.na(D$HRC_CORE) & !is.na(D[[k]])))
}))
rel[, q_value := p.adjust(p_value, method = "BH")]
rel[, interpretation_limit := "association between deltas; does not establish that composition causes the program change"]
o8 <- pm_out("04_results", "cell_source", "program_vs_compartment_association.csv"); write_csv_atomic(rel, o8); outs <- c(outs, o8)
cat("\n=== host compartment deltas (PM - primary), GSE314812 ===\n")
print(adm[, .(compartment, n, mean = round(mean, 4), ci_low = round(ci_low, 4), ci_high = round(ci_high, 4), q_value = signif(q_value, 3))])
step_end(outputs = outs)
