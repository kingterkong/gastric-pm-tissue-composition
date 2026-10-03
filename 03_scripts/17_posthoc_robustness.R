#!/usr/bin/env Rscript
# Post-hoc robustness analyses added after the frozen primary analysis. These analyses
# do not alter the primary endpoint: they test whether its direction depends on the
# scoring algorithm, whether an independent epithelial marker set behaves similarly,
# and whether the frozen score tracks externally reported histological tumour purity.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
source(pm_src("03_scripts", "R", "score_utils.R"))
suppressPackageStartupMessages(library(singscore))
cfg <- pm_config()

prog_f <- pm_in("02_data_processed", "frozen_programs", "program_definitions.rds")
bulk_f <- pm_in("02_data_processed", "bulk", "GSE314812_counts.rds")
score_f <- pm_in("04_results", "validation", "GSE237876_program_scores_by_specimen.csv")
primary_f <- pm_in("04_results", "primary", "primary_endpoint.csv")
step_begin("17_posthoc_robustness", inputs = c(prog_f, bulk_f, score_f, primary_f),
           seed = cfg$project$seed + 17, uses_protocol = FALSE)

boot_dz <- function(x, n_boot = 5000, seed = 1L) {
  x <- x[is.finite(x)]
  set.seed(seed)
  bs <- replicate(n_boot, {
    z <- sample(x, length(x), replace = TRUE)
    if (sd(z) == 0) NA_real_ else mean(z) / sd(z)
  })
  unname(quantile(bs, c(0.025, 0.975), na.rm = TRUE))
}

summarise_delta <- function(delta, method, target, scale_note, seed) {
  s <- paired_summary(delta, n_boot = 5000, seed = seed)
  b <- boot_dz(delta, seed = seed + 1000L)
  cbind(data.table(method = method, target = target, scale_note = scale_note,
                   analysis_status = "post hoc", dz_ci_low = b[1], dz_ci_high = b[2]), s)
}

# ----- Primary cohort: four algorithms on the identical 100-gene programme ----------
progs <- readRDS(prog_f)
hrc <- progs[["HRC_CORE"]]$genes
b <- readRDS(bulk_f)
u <- frozen_universe(b$counts, b$genes$symbol_for_scoring, min_cpm = 1, min_frac = 0.25)
cs <- collapse_to_symbols(b$counts, b$genes$symbol_for_scoring, u)
lib <- colSums(cs$expr); lib[lib == 0] <- 1
lcpm <- log2(t(t(cs$expr) / lib) * 1e6 + 1)
hrc_present <- intersect(hrc, rownames(lcpm))
stopifnot(length(hrc_present) / length(hrc) >= 0.80)

score_rank <- rank_scores(lcpm, rownames(lcpm), hrc, min_frac = 0.80, min_n = 10)$score
score_mean <- colMeans(lcpm[hrc_present, , drop = FALSE])
gene_sd <- apply(lcpm[hrc_present, , drop = FALSE], 1, sd)
z <- sweep(lcpm[hrc_present, , drop = FALSE], 1,
           rowMeans(lcpm[hrc_present, , drop = FALSE]), "-")
z <- sweep(z, 1, gene_sd, "/")
z[!is.finite(z)] <- 0
score_z <- colMeans(z)
ranked <- rankGenes(lcpm)
score_sing <- simpleScore(ranked, upSet = hrc_present, centerScore = TRUE,
                          knownDirection = TRUE)$TotalScore

# Independent, canonical epithelial marker set. None is part of HRC_CORE.
epi <- c("EPCAM", "KRT8", "KRT18", "KRT19")
stopifnot(length(intersect(epi, hrc)) == 0L, all(epi %in% rownames(lcpm)))
score_epi_rank <- rank_scores(lcpm, rownames(lcpm), epi, min_frac = 1, min_n = 4)$score
score_epi_mean <- colMeans(lcpm[epi, , drop = FALSE])

score_dt <- data.table(expression_column = colnames(lcpm),
                       hrc_mean_rank = score_rank,
                       hrc_singscore = as.numeric(score_sing),
                       hrc_mean_logcpm = score_mean,
                       hrc_mean_gene_z = score_z,
                       epithelial_marker_mean_rank = score_epi_rank,
                       epithelial_marker_mean_logcpm = score_epi_mean)
meta <- as.data.table(b$samples)[, .(expression_column, canonical_patient_id, visit,
                                     lesion_role, eligible_primary)]
score_dt <- merge(meta, score_dt, by = "expression_column")
pair <- score_dt[eligible_primary == "TRUE" | eligible_primary == TRUE]
stopifnot(pair[, .N, by = .(canonical_patient_id, lesion_role)][, all(N == 1)])
methods <- c("hrc_mean_rank", "hrc_singscore", "hrc_mean_logcpm", "hrc_mean_gene_z",
             "epithelial_marker_mean_rank", "epithelial_marker_mean_logcpm")
wide <- dcast(pair, canonical_patient_id ~ lesion_role, value.var = methods)
delta <- data.table(canonical_patient_id = wide$canonical_patient_id)
for (m in methods) {
  delta[, (m) := wide[[paste0(m, "_solid_peritoneal_metastasis")]] -
                       wide[[paste0(m, "_primary_tumor")]]]
}

labels <- c(
  hrc_mean_rank = "HRC: frozen mean-rank score",
  hrc_singscore = "HRC: singscore",
  hrc_mean_logcpm = "HRC: mean log-CPM",
  hrc_mean_gene_z = "HRC: mean gene-wise z score",
  epithelial_marker_mean_rank = "Independent epithelial markers: mean-rank score",
  epithelial_marker_mean_logcpm = "Independent epithelial markers: mean log-CPM"
)
targets <- c(rep("100-gene HRC programme", 4), rep("EPCAM/KRT8/KRT18/KRT19", 2))
scales <- c("relative rank", "centered singscore", "log2(CPM + 1)",
            "gene-wise standardized expression", "relative rank", "log2(CPM + 1)")
robust <- rbindlist(lapply(seq_along(methods), function(i) {
  summarise_delta(delta[[methods[i]]], labels[[methods[i]]], targets[i], scales[i],
                  cfg$project$seed + 17L + i)
}), fill = TRUE)
robust[, marker_overlap_with_hrc := ifelse(grepl("Independent", method), 0L, NA_integer_)]

# The first row must reproduce the frozen primary result to numerical precision.
frozen <- fread(primary_f)
stopifnot(abs(robust[1, mean] - frozen[1, mean]) < 1e-12,
          abs(robust[1, p_value] - frozen[1, p_value]) < 1e-12)

o1 <- pm_out("04_results", "sensitivity", "scoring_method_robustness.csv")
write_csv_atomic(robust[target == "100-gene HRC programme"], o1)
o2 <- pm_out("04_results", "sensitivity", "independent_epithelial_marker_check.csv")
write_csv_atomic(robust[target == "EPCAM/KRT8/KRT18/KRT19"], o2)
o3 <- pm_out("04_results", "sensitivity", "posthoc_patient_deltas.csv")
write_csv_atomic(delta, o3)

# ----- External histological-purity check in GSE237876 -------------------------------
ext <- fread(score_f)
ext[, purity := as.numeric(purity_published)]
ext <- ext[is.finite(purity) & is.finite(HRC_CORE)]
ext[, `:=`(purity_centered = purity - mean(purity),
           hrc_centered = HRC_CORE - mean(HRC_CORE)), by = canonical_patient_id]
stopifnot(nrow(ext) == 66L, uniqueN(ext$canonical_patient_id) == 14L)

raw_cor <- cor.test(ext$purity, ext$HRC_CORE, method = "spearman", exact = FALSE)
within_cor <- cor.test(ext$purity_centered, ext$hrc_centered, method = "spearman", exact = FALSE)
fe <- lm(HRC_CORE ~ purity + factor(canonical_patient_id), data = ext)
fe_tab <- summary(fe)$coefficients["purity", ]
fe_ci <- confint(fe)["purity", ]

# Patient-cluster bootstrap of the within-patient (fixed-effects) slope.
patients <- unique(ext$canonical_patient_id)
cross <- ext[, .(xx = sum(purity_centered^2), xy = sum(purity_centered * hrc_centered)),
             by = canonical_patient_id]
set.seed(cfg$project$seed + 1701L)
bs_slope <- replicate(5000, {
  take <- sample(patients, length(patients), replace = TRUE)
  d <- cross[match(take, canonical_patient_id)]
  sum(d$xy) / sum(d$xx)
})
bs_ci <- unname(quantile(bs_slope, c(0.025, 0.975), na.rm = TRUE))

# Patient-site aggregation is a sensitivity analysis for multiple regions.
site <- ext[, .(purity = mean(purity), HRC_CORE = mean(HRC_CORE)),
            by = .(canonical_patient_id, lesion_role)]
site[, `:=`(purity_centered = purity - mean(purity),
            hrc_centered = HRC_CORE - mean(HRC_CORE)), by = canonical_patient_id]
site_fe <- lm(HRC_CORE ~ purity + factor(canonical_patient_id), data = site)
site_ci <- confint(site_fe)["purity", ]

# The five eligible primary-PM pairs are shown transparently; this small contrast does
# not establish a consistent purity shift between sites.
eligible <- ext[eligible_validation == TRUE &
                lesion_role %in% c("primary_tumor", "solid_peritoneal_metastasis")]
eligible_site <- eligible[, .(purity = mean(purity), HRC_CORE = mean(HRC_CORE)),
                          by = .(canonical_patient_id, lesion_role)]
pw <- dcast(eligible_site, canonical_patient_id ~ lesion_role, value.var = c("purity", "HRC_CORE"))
purity_delta <- pw$purity_solid_peritoneal_metastasis - pw$purity_primary_tumor
hrc_delta <- pw$HRC_CORE_solid_peritoneal_metastasis - pw$HRC_CORE_primary_tumor
purity_pair <- paired_summary(purity_delta, n_boot = 5000, seed = cfg$project$seed + 1702L)
pair_cor <- cor.test(purity_delta, hrc_delta, method = "spearman", exact = FALSE)

purity_results <- rbindlist(list(
  data.table(analysis = "all specimens: raw Spearman correlation", n_specimens = nrow(ext),
             n_patients = uniqueN(ext$canonical_patient_id), estimate = unname(raw_cor$estimate),
             ci_low = NA_real_, ci_high = NA_real_, p_value = raw_cor$p.value,
             estimand = "Spearman rho; specimens not treated as independent for the main inference"),
  data.table(analysis = "all specimens: within-patient centered Spearman correlation", n_specimens = nrow(ext),
             n_patients = uniqueN(ext$canonical_patient_id), estimate = unname(within_cor$estimate),
             ci_low = NA_real_, ci_high = NA_real_, p_value = within_cor$p.value,
             estimand = "Spearman rho on patient-centered values; descriptive"),
  data.table(analysis = "all specimens: patient fixed-effects slope", n_specimens = nrow(ext),
             n_patients = uniqueN(ext$canonical_patient_id), estimate = unname(coef(fe)["purity"]),
             ci_low = fe_ci[1], ci_high = fe_ci[2], p_value = fe_tab[4],
             estimand = "HRC-score change per unit published histological purity; parametric CI and P"),
  data.table(analysis = "all specimens: patient fixed-effects slope, cluster bootstrap CI", n_specimens = nrow(ext),
             n_patients = uniqueN(ext$canonical_patient_id), estimate = unname(coef(fe)["purity"]),
             ci_low = bs_ci[1], ci_high = bs_ci[2], p_value = mean(bs_slope <= 0) * 2,
             estimand = "HRC-score change per unit purity; 5000 patient-cluster bootstrap resamples"),
  data.table(analysis = "patient-site aggregated sensitivity: fixed-effects slope", n_specimens = nrow(site),
             n_patients = uniqueN(site$canonical_patient_id), estimate = unname(coef(site_fe)["purity"]),
             ci_low = site_ci[1], ci_high = site_ci[2],
             p_value = summary(site_fe)$coefficients["purity", 4],
             estimand = "HRC-score change per unit purity after averaging regions within patient and site"),
  data.table(analysis = "eligible primary-PM pairs: mean purity difference", n_specimens = nrow(eligible_site),
             n_patients = length(purity_delta), estimate = purity_pair$mean,
             ci_low = purity_pair$ci_low, ci_high = purity_pair$ci_high, p_value = purity_pair$p_value,
             estimand = "mean within-patient purity difference (solid PM minus primary)"),
  data.table(analysis = "eligible primary-PM pairs: correlation of purity and HRC differences", n_specimens = nrow(eligible_site),
             n_patients = length(purity_delta), estimate = unname(pair_cor$estimate),
             ci_low = NA_real_, ci_high = NA_real_, p_value = pair_cor$p.value,
             estimand = "Spearman rho between paired purity and HRC differences")
), fill = TRUE)
purity_results[, analysis_status := "post hoc"]

o4 <- pm_out("04_results", "validation", "GSE237876_published_purity_association.csv")
write_csv_atomic(purity_results, o4)
o5 <- pm_out("04_results", "validation", "GSE237876_purity_plot_values.csv")
write_csv_atomic(ext[, .(original_sample_id, canonical_patient_id, lesion_role, purity,
                         HRC_CORE, purity_centered, hrc_centered)], o5)
o6 <- pm_out("04_results", "validation", "GSE237876_paired_purity_deltas.csv")
write_csv_atomic(data.table(canonical_patient_id = pw$canonical_patient_id,
                            purity_delta = purity_delta, hrc_delta = hrc_delta), o6)

cat("\n=== Alternative scoring and independent epithelial markers ===\n")
print(robust[, .(method, n, mean = round(mean, 4), ci_low = round(ci_low, 4),
                 ci_high = round(ci_high, 4), p_value = signif(p_value, 4),
                 n_negative, cohen_dz = round(cohen_dz, 3))])
cat("\n=== External published histological purity ===\n")
print(purity_results[, .(analysis, n_specimens, n_patients, estimate = round(estimate, 4),
                         ci_low = round(ci_low, 4), ci_high = round(ci_high, 4),
                         p_value = signif(p_value, 4))])
step_end(outputs = c(o1, o2, o3, o4, o5, o6))
