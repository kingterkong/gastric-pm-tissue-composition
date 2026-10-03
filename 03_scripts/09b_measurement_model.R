# Stage 9b: measurement-model calibration (plan section 17.2). Using calibration donors
# only, builds pseudo-bulk mixtures with KNOWN malignant-cell fractions and measures how
# far the frozen program score moves per unit change in that fraction. This converts the
# observed bulk change into a statement about how much purity difference would be needed
# to produce it, without claiming any biology.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
source(pm_src("03_scripts", "R", "score_utils.R"))
cfg <- pm_config()

prog_f <- pm_in("02_data_processed", "frozen_programs", "program_definitions.rds")
cells_f <- pm_in("02_data_processed", "singlecell", "GSE183904_cells.csv.gz")
counts_f <- pm_in("02_data_processed", "singlecell", "GSE183904_counts_qc.rds")
prim_f <- pm_in("04_results", "primary", "primary_endpoint.csv")
step_begin("09b_measurement_model", inputs = c(prog_f, cells_f, counts_f, prim_f),
           seed = cfg$project$seed + 96, uses_protocol = TRUE)

progs <- readRDS(prog_f); gene_sets <- lapply(progs, `[[`, "genes")
cells <- fread(cells_f)
counts <- readRDS(counts_f)
stopifnot(identical(colnames(counts), cells$cell_id))

# calibration donors only: gastric primary tumours of GSE183904 (development side)
use <- cells[lesion_role == "primary_tumor"]
EPI_CLASSES <- c("malignant_cnv_high", "tumour_epithelial_cnv_low_or_na", "epithelial_normal_sample")
use[, is_epithelial := cell_class %in% EPI_CLASSES]
elig <- use[, .(n_mal = sum(cell_class == "malignant_cnv_high"),
                n_epi = sum(is_epithelial), n_nonepi = sum(!is_epithelial),
                n_other = sum(cell_class != "malignant_cnv_high")), by = sample_id][n_mal >= 100 & n_other >= 200]
FRACS <- c(0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60)
N_CELLS <- 300; N_REP <- 10
set.seed(cfg$project$seed + 96)

build_mixtures <- function(axis) {
  rbindlist(lapply(elig$sample_id, function(s) {
    if (axis == "malignant_fraction") {
      inpool <- use[sample_id == s & cell_class == "malignant_cnv_high", cell_id]
      outpool <- use[sample_id == s & cell_class != "malignant_cnv_high", cell_id]
    } else {
      inpool <- use[sample_id == s & is_epithelial == TRUE, cell_id]
      outpool <- use[sample_id == s & is_epithelial == FALSE, cell_id]
    }
    rbindlist(lapply(FRACS, function(f) {
      ni <- round(N_CELLS * f); no <- N_CELLS - ni
      if (ni > length(inpool) || no > length(outpool)) return(NULL)
      rbindlist(lapply(seq_len(N_REP), function(r) {
        idx <- c(sample(inpool, ni), sample(outpool, no))
        data.table(sample_id = s, axis = axis, fraction = f, rep = r, col = list(match(idx, colnames(counts))))
      }))
    }))
  }))
}
mix <- rbind(build_mixtures("malignant_fraction"), build_mixtures("epithelial_fraction"))
if (!nrow(mix)) stop("no calibration donor met the mixture requirements")

pseudo <- do.call(cbind, lapply(seq_len(nrow(mix)), function(i) Matrix::rowSums(counts[, mix$col[[i]], drop = FALSE])))
colnames(pseudo) <- sprintf("%s_%s_f%03d_r%02d", mix$sample_id, substr(mix$axis, 1, 3), round(mix$fraction * 100), mix$rep)
sm <- score_matrix(as.matrix(pseudo), rownames(counts), gene_sets)
sc <- cbind(mix[, .(sample_id, axis, fraction, rep)], sm$scores[, -1])
o1 <- pm_out("04_results", "sensitivity", "pseudobulk_mixture_scores.csv"); write_csv_atomic(sc, o1)

# slope of each evaluable program score on the known malignant fraction (donor-adjusted)
ev <- sm$coverage[evaluable == TRUE, program_id]
slopes <- rbindlist(lapply(unique(sc$axis), function(ax) rbindlist(lapply(ev, function(p) {
  d <- sc[axis == ax]
  m <- stats::lm(stats::as.formula(paste0("`", p, "` ~ fraction + sample_id")), data = d)
  s <- summary(m)$coefficients["fraction", ]
  ci <- stats::confint(m)["fraction", ]
  data.table(axis = ax, program_id = p, slope_per_unit_fraction = unname(s[1]), ci_low = ci[1], ci_high = ci[2],
             p_value = unname(s[4]), score_change_per_10pct_fraction = unname(s[1]) * 0.10,
             n_mixtures = nrow(d), n_donors = uniqueN(d$sample_id))
}))))
o2 <- pm_out("04_results", "sensitivity", "program_score_vs_malignant_fraction.csv"); write_csv_atomic(slopes, o2)

# how large a purity difference would reproduce the observed bulk change?
prim <- fread(prim_f)
interp <- rbindlist(lapply(unique(slopes$axis), function(ax) {
  hs <- slopes[axis == ax & program_id == "HRC_CORE"]
  data.table(axis = ax,
    observed_bulk_delta = prim$mean[1], observed_ci_low = prim$ci_low[1], observed_ci_high = prim$ci_high[1],
    slope_per_unit_fraction = hs$slope_per_unit_fraction[1],
    implied_fraction_difference = prim$mean[1] / hs$slope_per_unit_fraction[1],
    implied_ci_low = prim$ci_low[1] / hs$slope_per_unit_fraction[1],
    implied_ci_high = prim$ci_high[1] / hs$slope_per_unit_fraction[1],
    interpretation = paste("A difference of this size in the stated cell fraction between the paired specimens would",
                           "reproduce the observed tissue-level change with no change in the program inside cells.",
                           "This is a measurement-model statement, not evidence that such a difference exists;",
                           "it bounds what the bulk contrast can distinguish."),
    calibration_source = "GSE183904 gastric primary tumours (development side); 300-cell mixtures, 10 replicates per fraction")
}))
o3 <- pm_out("04_results", "sensitivity", "purity_equivalence.csv"); write_csv_atomic(interp, o3)

# ---- does ONE compositional shift explain every program at once? --------------------
# For each evaluable program, convert its observed bulk delta into the epithelial-fraction
# difference that would reproduce it. If the programs agree on a single value, one
# compositional shift is a sufficient explanation for all of them.
sec <- fread(pm_in("04_results", "primary", "secondary_families.csv"))
obs <- rbind(prim[, .(program_id, mean, ci_low, ci_high, p_value)],
             sec[!is.na(mean), .(program_id, mean, ci_low, ci_high, p_value)])
epi_slopes <- slopes[axis == "epithelial_fraction", .(program_id, slope = slope_per_unit_fraction)]
cons <- merge(obs, epi_slopes, by = "program_id")
cons <- cons[abs(slope) > 0.05]   # a program whose score barely tracks composition cannot calibrate it
cons[, implied_epithelial_fraction_difference := mean / slope]
cons[, implied_ci_low := pmin(ci_low / slope, ci_high / slope)]
cons[, implied_ci_high := pmax(ci_low / slope, ci_high / slope)]
setorder(cons, implied_epithelial_fraction_difference)
cons[, common_explanation_note := paste("if a single difference in epithelial content explained every program,",
                                        "these implied values would agree")]
o4 <- pm_out("04_results", "sensitivity", "single_composition_consistency.csv"); write_csv_atomic(cons, o4)
cat("\n=== implied epithelial-fraction difference per program (one common shift?) ===\n")
print(cons[, .(program_id, observed = round(mean, 4), slope = round(slope, 3),
               implied = round(implied_epithelial_fraction_difference, 3),
               implied_lo = round(implied_ci_low, 3), implied_hi = round(implied_ci_high, 3))])

cat("\n=== score change per unit cell fraction (calibration mixtures) ===\n")
print(dcast(slopes, program_id ~ axis, value.var = "slope_per_unit_fraction")[order(-epithelial_fraction)])
cat("\n=== cell-fraction difference that would reproduce the observed bulk change ===\n")
print(interp[, .(axis, observed_bulk_delta = round(observed_bulk_delta, 4),
                 slope = round(slope_per_unit_fraction, 4),
                 implied_fraction_difference = round(implied_fraction_difference, 3),
                 implied_ci_low = round(implied_ci_low, 3), implied_ci_high = round(implied_ci_high, 3))])
step_end(outputs = c(o1, o2, o3, o4))
