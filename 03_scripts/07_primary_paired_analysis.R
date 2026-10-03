# Stage 7: the pre-specified primary analysis. Scores every GSE314812 specimen with the
# frozen programs, forms one delta per eligible patient, and runs the single primary test
# plus the pre-specified secondary families, adjusted models and sensitivity analyses.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
source(pm_src("03_scripts", "R", "score_utils.R"))
cfg <- pm_config()

prog_f <- pm_in("02_data_processed", "frozen_programs", "program_definitions.rds")
lock_f <- pm_in("00_admin", "protocol_lock.json")
bulk_f <- pm_in("02_data_processed", "bulk", "GSE314812_counts.rds")
step_begin("07_primary_paired_analysis", inputs = c(prog_f, lock_f, bulk_f), seed = cfg$project$seed + 7, uses_protocol = TRUE)

progs <- readRDS(prog_f)
gene_sets <- lapply(progs, `[[`, "genes")
b <- readRDS(bulk_f)
samples <- as.data.table(b$samples)
symbols <- b$genes$symbol_for_scoring

sm <- score_matrix(b$counts, symbols, gene_sets)
scores <- sm$scores
setnames(scores, "sample", "expression_column")
scores <- merge(samples[, .(expression_column, canonical_patient_id, visit, lesion_role, preservation, eligible_primary,
                            sex, age_if_public, stage, histology)], scores, by = "expression_column")
o_sc <- pm_out("04_results", "primary", "GSE314812_program_scores_by_specimen.csv"); write_csv_atomic(scores, o_sc)
o_cov <- pm_out("04_results", "primary", "GSE314812_program_coverage.csv")
write_csv_atomic(cbind(sm$coverage, universe_size = sm$universe_size), o_cov)

prog_ids <- names(gene_sets)
evaluable <- sm$coverage[evaluable == TRUE, program_id]
stopifnot("HRC_CORE" %in% evaluable)

# ---- patient-level deltas (one pair per patient; earliest eligible visit) ------------
pair <- scores[eligible_primary == TRUE]
stopifnot(pair[, .N, by = .(canonical_patient_id, lesion_role)][, all(N == 1)])
wide <- dcast(pair, canonical_patient_id + visit ~ lesion_role, value.var = c(prog_ids, "preservation", "expression_column"))
delta_of <- function(p) wide[[paste0(p, "_solid_peritoneal_metastasis")]] - wide[[paste0(p, "_primary_tumor")]]
D <- data.table(canonical_patient_id = wide$canonical_patient_id, visit = wide$visit)
for (p in prog_ids) D[, (p) := delta_of(p)]
D[, preservation_primary := wide$preservation_primary_tumor]
D[, preservation_pm := wide$preservation_solid_peritoneal_metastasis]
D[, preservation_change := fcase(preservation_primary == preservation_pm, "same",
                                 preservation_primary == "FFPE" & preservation_pm == "fresh-frozen", "FFPE_to_FF",
                                 default = "FF_to_FFPE")]
o_d <- pm_out("04_results", "primary", "GSE314812_patient_deltas.csv"); write_csv_atomic(D, o_d)

# ---- PRIMARY TEST -------------------------------------------------------------------
prim <- paired_summary(D$HRC_CORE, n_boot = 5000, seed = cfg$project$seed)
prim[, `:=`(program_id = "HRC_CORE", family = "primary", estimand = "mean within-patient delta (solid PM - primary), patient-weighted")]
setcolorder(prim, c("program_id", "family", "estimand"))
o_p <- pm_out("04_results", "primary", "primary_endpoint.csv"); write_csv_atomic(prim, o_p)

# ---- secondary families -------------------------------------------------------------
fam_of <- function(p) {
  r <- progs[[p]]$role
  fcase(r == "control", "secondary_A_control_programs", r == "negative_control", "secondary_A_control_programs",
        r == "composition_covariate", "secondary_B_host_composition", r == "sensitivity", "sensitivity_program_definition",
        default = "other")
}
sec <- rbindlist(lapply(setdiff(prog_ids, "HRC_CORE"), function(p) {
  if (!p %in% evaluable) return(data.table(program_id = p, family = fam_of(p), n = NA_integer_, mean = NA_real_,
                                           note = "not evaluable on this platform (coverage rule)"))
  s <- paired_summary(D[[p]], n_boot = 5000, seed = cfg$project$seed)
  cbind(data.table(program_id = p, family = fam_of(p)), s)
}), fill = TRUE)
sec[family %in% c("secondary_A_control_programs", "secondary_B_host_composition") & !is.na(p_value),
    q_value := p.adjust(p_value, method = "BH"), by = family]
o_s <- pm_out("04_results", "primary", "secondary_families.csv"); write_csv_atomic(sec, o_s)

# ---- pre-specified adjusted models (each with at most two extra covariates) ---------
fit <- function(formula, label) {
  m <- try(stats::lm(formula, data = D), silent = TRUE)
  if (inherits(m, "try-error")) return(data.table(model = label, term = NA_character_, estimate = NA_real_, note = "fit failed"))
  s <- summary(m)$coefficients
  r <- as.data.table(s, keep.rownames = "term")
  setnames(r, c("term", "estimate", "std_error", "t_value", "p_value"))
  ci <- stats::confint(m)
  r[, `:=`(ci_low = ci[, 1], ci_high = ci[, 2], model = label, n = stats::nobs(m),
           design_rank = qr(stats::model.matrix(m))$rank, n_coef = ncol(stats::model.matrix(m)))]
  r[]
}
adj <- rbindlist(list(
  fit(HRC_CORE ~ 1, "unadjusted (primary)"),
  fit(HRC_CORE ~ preservation_change, "preservation change"),
  # pre-specified as mesothelium + adipocyte; the mesothelium set is not evaluable on this
  # platform under the frozen coverage rule, so it runs with adipocyte alone (deviations DEV-002)
  fit(HRC_CORE ~ COMPART_ADIPOCYTE, "host compartment (adipocyte delta; mesothelium not evaluable)"),
  fit(HRC_CORE ~ COMPART_FIBROBLAST, "fibroblast compartment delta"),
  fit(HRC_CORE ~ CTRL_STROMA, "stromal abundance delta"),
  fit(HRC_CORE ~ CTRL_PROLIFERATION, "proliferation delta")
), fill = TRUE)
o_a <- pm_out("04_results", "primary", "adjusted_models.csv"); write_csv_atomic(adj, o_a)

# ---- sensitivity --------------------------------------------------------------------
loo <- leave_one_out(D$HRC_CORE, D$canonical_patient_id)
o_l <- pm_out("04_results", "sensitivity", "leave_one_patient_out.csv"); write_csv_atomic(loo, o_l)

same_pres <- D[preservation_change == "same"]
sens <- rbindlist(list(
  cbind(data.table(analysis = "primary (all eligible patients)"), paired_summary(D$HRC_CORE, n_boot = 5000, seed = cfg$project$seed)),
  cbind(data.table(analysis = "same-preservation pairs only"), paired_summary(same_pres$HRC_CORE, n_boot = 5000, seed = cfg$project$seed)),
  cbind(data.table(analysis = "epithelium-restricted program subset"), paired_summary(D$SENS_HRC_EPI_RESTRICTED, n_boot = 5000, seed = cfg$project$seed)),
  cbind(data.table(analysis = "alternative source definition (epiHR)"), paired_summary(D$SENS_HRC_EPIHR, n_boot = 5000, seed = cfg$project$seed))
), fill = TRUE)

# GC05: swap the pre-treatment pair for the post-treatment recurrence pair
alt <- scores[canonical_patient_id != "GC05" & eligible_primary == TRUE]
gc05b <- scores[canonical_patient_id == "GC05" & visit == "2" & lesion_role %in% c("primary_tumor", "solid_peritoneal_metastasis")]
if (nrow(gc05b) == 2) {
  a <- rbind(alt, gc05b)
  w2 <- dcast(a, canonical_patient_id ~ lesion_role, value.var = "HRC_CORE")
  d2 <- w2$solid_peritoneal_metastasis - w2$primary_tumor
  sens <- rbind(sens, cbind(data.table(analysis = "GC05 post-treatment pair instead of pre-treatment"),
                            paired_summary(d2, n_boot = 5000, seed = cfg$project$seed)), fill = TRUE)
}
o_se <- pm_out("04_results", "sensitivity", "primary_sensitivity.csv"); write_csv_atomic(sens, o_se)

cat("\n=== PRIMARY ENDPOINT (HRC_CORE, GSE314812) ===\n"); print(prim)
cat("\n=== secondary families ===\n"); print(sec[, .(program_id, family, n, mean = round(mean, 4), ci_low = round(ci_low, 4), ci_high = round(ci_high, 4), p_value = signif(p_value, 3), q_value = signif(q_value, 3))])
cat("\n=== sensitivity ===\n"); print(sens[, .(analysis, n, mean = round(mean, 4), ci_low = round(ci_low, 4), ci_high = round(ci_high, 4), p_value = signif(p_value, 3))])
cat("\n=== adjusted ===\n"); print(adj[, .(model, term, estimate = round(estimate, 4), ci_low = round(ci_low, 4), ci_high = round(ci_high, 4), p_value = signif(p_value, 3), n, design_rank, n_coef)])
step_end(outputs = c(o_sc, o_cov, o_d, o_p, o_s, o_a, o_l, o_se))
