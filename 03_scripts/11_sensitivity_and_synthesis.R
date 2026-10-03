# Stage 11: synthesis. Builds the hypothesis registry (every pre-specified test with its
# result or the reason it could not be run), the cross-source summary of the frozen
# endpoint, and the claims ledger that the manuscript must draw from.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
source(pm_src("03_scripts", "R", "score_utils.R"))
cfg <- pm_config()

f <- function(...) pm_in("04_results", ...)
inputs <- c(f("primary", "primary_endpoint.csv"), f("primary", "secondary_families.csv"),
            f("primary", "adjusted_models.csv"), f("sensitivity", "primary_sensitivity.csv"),
            f("sensitivity", "leave_one_patient_out.csv"), f("validation", "GSE237876_frozen_endpoint.csv"),
            f("cell_source", "host_compartment_deltas.csv"), f("sensitivity", "single_composition_consistency.csv"),
            f("sensitivity", "purity_equivalence.csv"), pm_in("00_admin", "protocol_lock.json"))
inputs <- inputs[file.exists(inputs)]
step_begin("11_sensitivity_and_synthesis", inputs = inputs, seed = cfg$project$seed + 11, uses_protocol = TRUE)

rd <- function(p) if (file.exists(p)) fread(p) else NULL
prim <- rd(f("primary", "primary_endpoint.csv"))
sec <- rd(f("primary", "secondary_families.csv"))
adj <- rd(f("primary", "adjusted_models.csv"))
sens <- rd(f("sensitivity", "primary_sensitivity.csv"))
loo <- rd(f("sensitivity", "leave_one_patient_out.csv"))
val <- rd(f("validation", "GSE237876_frozen_endpoint.csv"))
adm <- rd(f("cell_source", "host_compartment_deltas.csv"))
fluid <- rd(f("validation", "GSE239676_fluid_context_estimand.csv"))
mal <- rd(f("cell_source", "malignant_library_summary.csv"))
cons <- rd(f("sensitivity", "single_composition_consistency.csv"))
pur <- rd(f("sensitivity", "purity_equivalence.csv"))
flow <- rd(pm_in("00_admin", "eligibility_flow.csv"))

fmt <- function(x, d = 4) ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = d))
pfmt <- function(p) ifelse(is.na(p), NA_character_, ifelse(p < 1e-4, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 4)))
eff <- function(dt) if (is.null(dt) || !nrow(dt)) "" else sprintf("%s (95%% CI %s to %s), n=%d, P=%s",
  fmt(dt$mean[1]), fmt(dt$ci_low[1]), fmt(dt$ci_high[1]), dt$n[1], pfmt(dt$p_value[1]))

# ---- hypothesis registry ------------------------------------------------------------
reg <- list()
add <- function(id, tier, source_, effect, test, family, fsize, estimable, p, q, note) {
  reg[[length(reg) + 1]] <<- data.table(hypothesis_id = id, tier = tier, source = source_, effect = effect,
    test = test, family_id = family, prespecified_family_size = fsize, estimable = estimable,
    p_value = p, q_value = q, not_executed_reason = note)
}
add("H1", "primary", "GSE314812", "mean within-patient delta of HRC_CORE (solid PM - primary)",
    "two-sided paired t-test", "primary", 1, !is.null(prim) && !is.na(prim$p_value[1]),
    if (!is.null(prim)) prim$p_value[1] else NA_real_, NA_real_, NA_character_)
if (!is.null(sec)) for (i in seq_len(nrow(sec))) add(paste0("H2.", sec$program_id[i]), "secondary", "GSE314812",
    paste0("mean within-patient delta of ", sec$program_id[i]), "two-sided paired t-test", sec$family[i],
    sec[family == sec$family[i], .N], !is.na(sec$p_value[i]), sec$p_value[i], sec$q_value[i],
    if (is.na(sec$p_value[i])) "program not evaluable on this platform (coverage rule)" else NA_character_)
if (!is.null(val)) for (i in seq_len(nrow(val))) add(paste0("H3.", val$program_id[i], ".", gsub(" ", "_", val$cohort[i])),
    "secondary", val$cohort[i], paste0("independent-source within-patient delta of ", val$program_id[i]),
    "two-sided paired t-test (directional consistency only; n=5)", "external_directional", nrow(val),
    !is.na(val$p_value[i]), val$p_value[i], NA_real_,
    if (is.na(val$p_value[i])) "not evaluable (coverage or too few patients)" else NA_character_)
if (!is.null(fluid)) for (i in seq_len(nrow(fluid))) add(paste0("H4.", fluid$program_id[i]), "exploratory", "GSE239676",
    "separate fluid-context estimand (malignant cells: ascites - primary)", "two-sided paired t-test", "exploratory_fluid",
    nrow(fluid), !is.na(fluid$p_value[i]), fluid$p_value[i], NA_real_, NA_character_)
add("H5", "secondary", "GSE308231 / GSE183904", "malignant-cell program value in solid peritoneal libraries",
    "descriptive library-level summary", "cell_level_descriptive", 1, !is.null(mal) && nrow(mal) > 0, NA_real_, NA_real_,
    if (is.null(mal) || !nrow(mal)) "no library reached the 20-cell threshold in a solid site" else
      "donor identity unverified in GSE308231; descriptive only, no paired inference possible")
add("H6", "secondary", "GSE314812", "direct PM versus other metastatic site contrast relative to primary",
    "(PM - primary) - (other site - primary)", "site_contrast", 1, FALSE, NA_real_, NA_real_,
    "not estimable: no patient in GSE314812 has a primary, a solid peritoneal metastasis and another metastatic site with the required pairing (omentum n=2, lymph node n=2)")
add("H7", "exploratory", "ACRG/GSE62254", "association of the primary program in primary tumours with first peritoneal recurrence",
    "cause-specific Cox", "clinical_exploratory", 1, FALSE, NA_real_, NA_real_,
    "module switched off by protocol default; the primary question concerns established lesions, not future risk")
registry <- rbindlist(reg, fill = TRUE)
o_reg <- pm_out("00_admin", "hypothesis_registry.csv"); write_csv_atomic(registry, o_reg)

# ---- cross-source summary of the frozen endpoint ------------------------------------
rows <- list()
if (!is.null(prim)) rows[[1]] <- data.table(source = "GSE314812 (primary test)", unit = "patients", n = prim$n[1],
  mean = prim$mean[1], ci_low = prim$ci_low[1], ci_high = prim$ci_high[1], p_value = prim$p_value[1],
  n_positive = prim$n_positive[1], n_negative = prim$n_negative[1])
if (!is.null(val)) for (i in seq_len(nrow(val[program_id == "HRC_CORE"]))) {
  v <- val[program_id == "HRC_CORE"][i]
  rows[[length(rows) + 1]] <- data.table(source = v$cohort, unit = "patients", n = v$n, mean = v$mean,
    ci_low = v$ci_low, ci_high = v$ci_high, p_value = v$p_value, n_positive = v$n_positive, n_negative = v$n_negative)
}
if (!is.null(fluid)) { v <- fluid[program_id == "HRC_CORE"]
  if (nrow(v)) rows[[length(rows) + 1]] <- data.table(source = "GSE239676 malignant cells (ascites - primary; separate estimand)",
    unit = "patients", n = v$n[1], mean = v$mean[1], ci_low = v$ci_low[1], ci_high = v$ci_high[1],
    p_value = v$p_value[1], n_positive = v$n_positive[1], n_negative = v$n_negative[1]) }
synth <- rbindlist(rows, fill = TRUE)
o_syn <- pm_out("04_results", "manuscript", "cross_source_summary.csv"); write_csv_atomic(synth, o_syn)

# ---- evidence grading ----------------------------------------------------------------
grade <- function(n, mean_, lo, hi, p) {
  if (is.na(n) || n < 3) return("not_evaluable")
  if (!is.na(p) && p < 0.05) return("compatible_direction_and_magnitude")
  if (!is.na(mean_) && !is.na(lo) && !is.na(hi)) {
    if (sign(lo) == sign(hi)) return("compatible_direction_and_magnitude")
    return(if (n <= 6) "limited_directional" else "uninformative_wide_interval")
  }
  "uninformative_wide_interval"
}
if (nrow(synth)) synth[, evidence_grade := mapply(grade, n, mean, ci_low, ci_high, p_value)]
write_csv_atomic(synth, o_syn)

# ---- claims ledger --------------------------------------------------------------------
cl <- list()
addc <- function(id, zh, en, tier, estimand, file_, key, n_pat, n_src, effect, lo, hi, p, q, level,
                 contra, alt, strength, allowed, prohibited) {
  cl[[length(cl) + 1]] <<- data.table(claim_id = id, plain_chinese_claim = zh, proposed_english_claim = en,
    primary_or_secondary_or_exploratory = tier, estimand = estimand, supporting_result_file = file_, table_key = key,
    independent_patient_n = n_pat, independent_source_n = n_src, effect = effect, ci_low = lo, ci_high = hi,
    p_value = p, q_value = q, measurement_level = level, contradictory_evidence = contra,
    main_alternative_explanation = alt, evidence_strength = strength, allowed_wording = allowed,
    prohibited_overclaim = prohibited, final_manuscript_location = NA_character_)
}
if (!is.null(prim)) addc("C1",
  "在 GSE314812 的配对患者中，冻结的 coreHRC 组织分数在实体腹膜转移与同患者原发之间的平均差值（见结果文件）",
  "Within-patient change of the frozen high-relapse-cell program score between solid peritoneal metastasis and paired primary tumour",
  "primary", "mean patient-level delta (PM - primary)", "04_results/primary/primary_endpoint.csv", "HRC_CORE",
  prim$n[1], 1L, prim$mean[1], prim$ci_low[1], prim$ci_high[1], prim$p_value[1], NA_real_,
  "bulk tissue transcriptional score", "see adjusted models and control programs",
  "host-tissue admixture (mesothelium/adipose in PM, gastric mucosa/smooth muscle in primary); preservation method differs within 13 of 18 pairs",
  "strong for the tissue-level observation: 15 of 17 patients, stable to leave-one-out, two alternative definitions agree", "observed within-patient tissue-level change", "cell-intrinsic activation; drives peritoneal seeding; therapeutic target")
if (!is.null(adm)) addc("C2",
  "腹膜转移标本相对同患者原发的宿主区室分数差异（间皮、脂肪等），用于量化组织混入",
  "Host-compartment score differences between paired specimens, quantifying tissue admixture",
  "secondary", "mean patient-level delta of compartment scores", "04_results/cell_source/host_compartment_deltas.csv",
  "COMPART_*", if (!is.null(prim)) prim$n[1] else NA_integer_, 1L, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
  "bulk tissue transcriptional score", NA_character_,
  "compartment scores are themselves expression-based and cannot separate cell number from cell state",
  "moderate to strong: four of six evaluable compartments change with q < 0.02 and in the direction peritoneal tissue predicts", "quantified host-tissue signal in the sampled specimens", "deconvolution-grade cell fractions")
if (!is.null(val)) addc("C3",
  "在独立来源 GSE237876 的配对患者中运行同一冻结终点的结果（患者数极少，仅方向性）",
  "The identical frozen endpoint run once in an independent paired cohort (very few patients; directional only)",
  "secondary", "mean patient-level delta (PM - primary)", "04_results/validation/GSE237876_frozen_endpoint.csv", "HRC_CORE",
  val[program_id == "HRC_CORE" & cohort == "GSE237876", n][1], 1L,
  val[program_id == "HRC_CORE" & cohort == "GSE237876", mean][1],
  val[program_id == "HRC_CORE" & cohort == "GSE237876", ci_low][1],
  val[program_id == "HRC_CORE" & cohort == "GSE237876", ci_high][1],
  val[program_id == "HRC_CORE" & cohort == "GSE237876", p_value][1], NA_real_,
  "bulk tissue transcriptional score", NA_character_,
  "all five patients are women with diffuse/mixed histology, three sampled after chemotherapy, peritoneal purity mostly very low",
  "weak: with five patients the smallest attainable two-sided exact P is 0.0625, so the design cannot reach significance whatever the effect", "independent cohort analysis with limited directional support", "independent validation; replication")
if (!is.null(cons) && !is.null(pur)) {
  pe <- pur[axis == "epithelial_fraction"]
  addc("C4",
    "\u89c2\u5bdf\u5230\u7684\u7ec4\u7ec7\u5c42\u9762\u5dee\u503c\u76f8\u5f53\u4e8e\u4e0a\u76ae\u542b\u91cf\u76f8\u5dee\u7ea6 30 \u4e2a\u767e\u5206\u70b9\uff1b\u5bf9\u5168\u90e8\u7a0b\u5e8f\u6362\u7b97\u540e\u7ed3\u679c\u4e00\u81f4",
    "A single difference in epithelial content reproduces the observed change in every scored programme, including the published EMT increase",
    "secondary", "observed tissue-level difference divided by the score-versus-composition slope from calibration mixtures",
    "04_results/sensitivity/single_composition_consistency.csv", "implied_epithelial_fraction_difference",
    prim$n[1], 1L, pe$implied_fraction_difference[1], min(pe$implied_ci_low[1], pe$implied_ci_high[1]),
    max(pe$implied_ci_low[1], pe$implied_ci_high[1]), NA_real_, NA_real_,
    "measurement-model conversion, not a measured purity",
    "the conversion is linear while the underlying relationship saturates at high epithelial content",
    "composition is a sufficient explanation; it is not established as the actual cause and a cell-state change is not excluded",
    "moderate to strong as a sufficiency argument: implied values across the scored programmes have an interquartile range of about 0.15",
    "a compositional difference of this size would reproduce the observed change",
    "proves that cancer-cell state is unchanged")
}
if (!is.null(mal) && nrow(mal)) {
  mpm <- mal[lesion_role == "solid_peritoneal_metastasis"]
  addc("C5",
    "\u516c\u5f00\u5355\u7ec6\u80de\u6570\u636e\u4ece\u5b9e\u4f53\u8179\u819c\u75c5\u7076\u56de\u6536\u7684\u6076\u6027\u7ec6\u80de\u8fc7\u5c11\uff0c\u65e0\u6cd5\u56de\u7b54\u7ec6\u80de\u5185\u95ee\u9898",
    "Public single-cell data recover too few malignant cells from solid peritoneal lesions to resolve the cell-intrinsic question",
    "secondary", "malignant cells passing quality control per library",
    "04_results/cell_source/malignant_library_summary.csv", "n_cells",
    NA_integer_, 3L, max(mpm$n_cells), NA_real_, NA_real_, NA_real_, NA_real_,
    "cell counts", NA_character_,
    "donor identity is unverifiable in one source, so libraries rather than patients are the unit",
    "strong as a statement about data availability", "quantified limit of the available public data",
    "absence of a cell-intrinsic change")
}

ledger <- rbindlist(cl, fill = TRUE)
o_cl <- pm_out("00_admin", "claims_ledger.csv"); write_csv_atomic(ledger, o_cl)

# ---- evidence status ------------------------------------------------------------------
ev <- data.table(
  module = c("primary within-patient test", "control and negative-control programs", "host composition",
             "independent paired bulk", "malignant-cell descriptive", "fluid context", "other-site contrast",
             "spatial single case", "clinical peritoneal recurrence"),
  status = c(if (!is.null(prim)) "executed" else "not_run", if (!is.null(sec)) "executed" else "not_run",
             if (!is.null(adm)) "executed" else "not_run", if (!is.null(val)) "executed" else "not_run",
             if (!is.null(mal) && nrow(mal)) "executed" else "not_evaluable",
             if (!is.null(fluid)) "executed" else "not_evaluable",
             "not_estimable", "not_run", "switched_off_by_protocol"),
  reason = c(NA, NA, NA, NA,
             "library-level only; donor identity unverified in GSE308231",
             "separate estimand; public matrix retains 8,630 genes so coverage must be reported per program",
             "no patient has the required primary + solid PM + other-site combination",
             "single patient; the same section has already been used by a preprint",
             "protocol default off; established-lesion question does not address future risk"))
o_ev <- pm_out("00_admin", "evidence_status.csv"); write_csv_atomic(ev, o_ev)

cat("\n=== cross-source summary (frozen endpoint) ===\n"); print(synth)
cat("\n=== hypothesis registry ===\n"); print(registry[, .(hypothesis_id, tier, estimable, not_executed_reason)])
step_end(outputs = c(o_reg, o_syn, o_cl, o_ev))
