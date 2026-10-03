# Stage 13: assemble the manuscript. Every study number is read from a result file, written
# to key_numbers.json and number_trace.csv, and substituted into the template by key. No
# number is typed into the prose.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
cfg <- pm_config()
R <- function(...) pm_in("04_results", ...)

tmpl_f <- pm_src("06_manuscript", "source", "manuscript_en_template.md")
abs_zh_f <- pm_src("06_manuscript", "source", "abstract_zh_template.md")
supp_f <- pm_src("06_manuscript", "source", "supplementary_template.md")
cover_f <- pm_src("06_manuscript", "source", "cover_letter_template.md")
inputs <- c(R("primary", "primary_endpoint.csv"), R("primary", "secondary_families.csv"),
            R("primary", "adjusted_models.csv"), R("primary", "GSE314812_program_coverage.csv"),
            R("sensitivity", "primary_sensitivity.csv"), R("sensitivity", "leave_one_patient_out.csv"),
            R("sensitivity", "single_composition_consistency.csv"), R("sensitivity", "purity_equivalence.csv"),
            R("sensitivity", "program_score_vs_malignant_fraction.csv"),
            R("validation", "GSE237876_frozen_endpoint.csv"), R("cell_source", "host_compartment_deltas.csv"),
            R("cell_source", "calibration_lineage_profile.csv"), R("cell_source", "malignant_library_summary.csv"),
            R("audit", "source_counts.csv"), R("audit", "matrix_type_audit.csv"),
            pm_in("00_admin", "eligibility_flow.csv"), tmpl_f, abs_zh_f, supp_f, cover_f,
            pm_in("07_references", "references.bib"))
step_begin("13_build_manuscript", inputs = inputs, seed = cfg$project$seed + 13, uses_protocol = TRUE)

prim <- fread(R("primary", "primary_endpoint.csv"))
sec <- fread(R("primary", "secondary_families.csv"))
adj <- fread(R("primary", "adjusted_models.csv"))
cov <- fread(R("primary", "GSE314812_program_coverage.csv"))
sens <- fread(R("sensitivity", "primary_sensitivity.csv"))
loo <- fread(R("sensitivity", "leave_one_patient_out.csv"))
cons <- fread(R("sensitivity", "single_composition_consistency.csv"))
pur <- fread(R("sensitivity", "purity_equivalence.csv"))
slp <- fread(R("sensitivity", "program_score_vs_malignant_fraction.csv"))
val <- fread(R("validation", "GSE237876_frozen_endpoint.csv"))
adm <- fread(R("cell_source", "host_compartment_deltas.csv"))
lin <- fread(R("cell_source", "calibration_lineage_profile.csv"))
mal <- fread(R("cell_source", "malignant_library_summary.csv"))
srcc <- fread(R("audit", "source_counts.csv"))
flow <- fread(pm_in("00_admin", "eligibility_flow.csv"))
mtx <- fread(R("audit", "matrix_type_audit.csv"))

trace <- list()
N <- function(key, value, file_, object, denom = NA_character_, unit = NA_character_, digits = 3) {
  disp <- if (is.character(value)) value else if (is.na(value)) "NA" else
    if (abs(value) >= 1 && value == round(value)) formatC(value, format = "d", big.mark = ",") else
      { v <- formatC(value, format = "f", digits = digits); if (grepl("^-0\\.0*$", v)) sub("^-", "", v) else v }
  trace[[length(trace) + 1]] <<- data.table(number_key = key, display_text = disp,
    raw_value = if (is.character(value)) NA_real_ else as.numeric(value),
    result_file = file_, object_key = object, denominator = denom, unit = unit,
    script = "03_scripts/13_build_manuscript.R")
  disp
}
P <- function(p) if (is.na(p)) "NA" else if (p < 1e-4) formatC(p, format = "e", digits = 1) else formatC(p, format = "f", digits = 4)
fmt3 <- function(x) { v <- formatC(x, format = "f", digits = 3); ifelse(v == "-0.000", "0.000", v) }
ci <- function(r) sprintf("%s to %s", fmt3(r$ci_low[1]), fmt3(r$ci_high[1]))
getp <- function(id) sec[program_id == id]
getc <- function(id) adm[compartment == id]

K <- list()
K$n_primary_patients <- N("n_primary_patients", prim$n[1], "04_results/primary/primary_endpoint.csv", "n", "patients", "patients")
K$n_specimens_314812 <- N("n_specimens_314812", srcc[dataset_id == "GSE314812", n_specimens], "04_results/audit/source_counts.csv", "n_specimens", "specimens", "specimens")
K$n_patients_314812 <- N("n_patients_314812", srcc[dataset_id == "GSE314812", n_patients_identifiable], "04_results/audit/source_counts.csv", "n_patients_identifiable", "patients", "patients")
K$n_pair_units <- N("n_pair_units", flow[dataset_id == "GSE314812" & grepl("pair units", step), n], "00_admin/eligibility_flow.csv", "pair units", "pair units", "pair units")
K$primary_mean <- N("primary_mean", prim$mean[1], "04_results/primary/primary_endpoint.csv", "mean", NA, "score units")
K$primary_ci <- N("primary_ci", ci(prim), "04_results/primary/primary_endpoint.csv", "ci_low..ci_high")
K$primary_p <- N("primary_p", P(prim$p_value[1]), "04_results/primary/primary_endpoint.csv", "p_value")
K$primary_neg <- N("primary_neg", prim$n_negative[1], "04_results/primary/primary_endpoint.csv", "n_negative", "patients", "patients")
K$primary_dz <- N("primary_dz", prim$cohen_dz[1], "04_results/primary/primary_endpoint.csv", "cohen_dz", NA, "standardised")
K$primary_boot_ci <- N("primary_boot_ci", sprintf("%s to %s", fmt3(prim$boot_ci_low[1]), fmt3(prim$boot_ci_high[1])), "04_results/primary/primary_endpoint.csv", "boot_ci")
K$primary_sign_p <- N("primary_sign_p", P(prim$sign_test_p[1]), "04_results/primary/primary_endpoint.csv", "sign_test_p")
K$primary_wilcox_p <- N("primary_wilcox_p", P(prim$wilcoxon_p[1]), "04_results/primary/primary_endpoint.csv", "wilcoxon_p")
K$loo_min <- N("loo_min", min(loo$mean), "04_results/sensitivity/leave_one_patient_out.csv", "min(mean)")
K$loo_max <- N("loo_max", max(loo$mean), "04_results/sensitivity/leave_one_patient_out.csv", "max(mean)")
K$loo_maxp <- N("loo_maxp", P(max(loo$p_value)), "04_results/sensitivity/leave_one_patient_out.csv", "max(p_value)")
K$universe <- N("universe", cov$universe_size[1], "04_results/primary/GSE314812_program_coverage.csv", "universe_size", "genes", "genes")
K$hrc_measured <- N("hrc_measured", cov[program_id == "HRC_CORE", n_used], "04_results/primary/GSE314812_program_coverage.csv", "n_used", "of 100", "genes")

for (nm in c("CTRL_EMT", "CTRL_STROMA", "CTRL_PROLIFERATION", "CTRL_IEG_STRESS", "NEG_YAP_STROMAL",
             "SENS_HRC_EPI_RESTRICTED", "SENS_HRC_EPIHR")) {
  r <- getp(nm); low <- tolower(nm)
  K[[paste0(low, "_mean")]] <- N(paste0(low, "_mean"), r$mean[1], "04_results/primary/secondary_families.csv", nm)
  K[[paste0(low, "_ci")]] <- N(paste0(low, "_ci"), ci(r), "04_results/primary/secondary_families.csv", nm)
  K[[paste0(low, "_p")]] <- N(paste0(low, "_p"), P(r$p_value[1]), "04_results/primary/secondary_families.csv", nm)
  K[[paste0(low, "_q")]] <- N(paste0(low, "_q"), P(r$q_value[1]), "04_results/primary/secondary_families.csv", nm)
}
for (nm in c("COMPART_GASTRIC_EPITHELIUM", "COMPART_FIBROBLAST", "COMPART_ADIPOCYTE", "COMPART_ENDOTHELIUM",
             "COMPART_IMMUNE", "COMPART_SMOOTH_MUSCLE")) {
  r <- getc(nm); low <- tolower(nm)
  K[[paste0(low, "_mean")]] <- N(paste0(low, "_mean"), r$mean[1], "04_results/cell_source/host_compartment_deltas.csv", nm)
  K[[paste0(low, "_ci")]] <- N(paste0(low, "_ci"), ci(r), "04_results/cell_source/host_compartment_deltas.csv", nm)
  K[[paste0(low, "_q")]] <- N(paste0(low, "_q"), P(r$q_value[1]), "04_results/cell_source/host_compartment_deltas.csv", nm)
}
# the number that matters is coverage in the SCORING universe, not in the platform annotation
K$meso_measured <- N("meso_measured", cov[program_id == "COMPART_MESOTHELIUM", n_used],
                     "04_results/primary/GSE314812_program_coverage.csv", "n_used", "of 40", "genes")
K$meso_defined <- N("meso_defined", cov[program_id == "COMPART_MESOTHELIUM", n_defined],
                    "04_results/primary/GSE314812_program_coverage.csv", "n_defined", NA, "genes")

am <- function(lbl) adj[model == lbl & term == "(Intercept)"]
K$adj_pres_mean <- N("adj_pres_mean", am("preservation change")$estimate[1], "04_results/primary/adjusted_models.csv", "preservation intercept")
K$adj_pres_p <- N("adj_pres_p", P(am("preservation change")$p_value[1]), "04_results/primary/adjusted_models.csv", "preservation intercept p")
K$adj_fib_mean <- N("adj_fib_mean", am("fibroblast compartment delta")$estimate[1], "04_results/primary/adjusted_models.csv", "fibroblast intercept")
K$adj_fib_ci <- N("adj_fib_ci", ci(am("fibroblast compartment delta")), "04_results/primary/adjusted_models.csv", "fibroblast intercept ci")
K$adj_fib_p <- N("adj_fib_p", P(am("fibroblast compartment delta")$p_value[1]), "04_results/primary/adjusted_models.csv", "fibroblast intercept p")
K$adj_stroma_mean <- N("adj_stroma_mean", am("stromal abundance delta")$estimate[1], "04_results/primary/adjusted_models.csv", "stroma intercept")
K$adj_stroma_ci <- N("adj_stroma_ci", ci(am("stromal abundance delta")), "04_results/primary/adjusted_models.csv", "stroma intercept ci")
K$adj_stroma_p <- N("adj_stroma_p", P(am("stromal abundance delta")$p_value[1]), "04_results/primary/adjusted_models.csv", "stroma intercept p")
K$adj_prolif_mean <- N("adj_prolif_mean", am("proliferation delta")$estimate[1], "04_results/primary/adjusted_models.csv", "proliferation intercept")
K$adj_prolif_p <- N("adj_prolif_p", P(am("proliferation delta")$p_value[1]), "04_results/primary/adjusted_models.csv", "proliferation intercept p")

sm <- function(a) sens[analysis == a]
K$sens_samepres_n <- N("sens_samepres_n", sm("same-preservation pairs only")$n[1], "04_results/sensitivity/primary_sensitivity.csv", "n", "patients", "patients")
K$sens_samepres_mean <- N("sens_samepres_mean", sm("same-preservation pairs only")$mean[1], "04_results/sensitivity/primary_sensitivity.csv", "mean")
K$sens_samepres_ci <- N("sens_samepres_ci", ci(sm("same-preservation pairs only")), "04_results/sensitivity/primary_sensitivity.csv", "ci")
K$sens_samepres_p <- N("sens_samepres_p", P(sm("same-preservation pairs only")$p_value[1]), "04_results/sensitivity/primary_sensitivity.csv", "p")
K$sens_gc05_mean <- N("sens_gc05_mean", sm("GC05 post-treatment pair instead of pre-treatment")$mean[1], "04_results/sensitivity/primary_sensitivity.csv", "mean")
K$sens_gc05_p <- N("sens_gc05_p", P(sm("GC05 post-treatment pair instead of pre-treatment")$p_value[1]), "04_results/sensitivity/primary_sensitivity.csv", "p")

v <- val[program_id == "HRC_CORE" & cohort == "GSE237876"]
v4 <- val[program_id == "HRC_CORE" & cohort == "GSE237876 excluding GCM05"]
K$ext_n <- N("ext_n", v$n[1], "04_results/validation/GSE237876_frozen_endpoint.csv", "n", "patients", "patients")
K$ext_mean <- N("ext_mean", v$mean[1], "04_results/validation/GSE237876_frozen_endpoint.csv", "mean")
K$ext_ci <- N("ext_ci", ci(v), "04_results/validation/GSE237876_frozen_endpoint.csv", "ci")
K$ext_p <- N("ext_p", P(v$p_value[1]), "04_results/validation/GSE237876_frozen_endpoint.csv", "p")
K$ext_neg <- N("ext_neg", v$n_negative[1], "04_results/validation/GSE237876_frozen_endpoint.csv", "n_negative", "of 5", "patients")
K$ext_minp <- N("ext_minp", formatC(2^-(v$n[1] - 1), format = "f", digits = 4), "04_results/validation/GSE237876_frozen_endpoint.csv", "max_attainable_two_sided_p")
K$ext4_n <- N("ext4_n", v4$n[1], "04_results/validation/GSE237876_frozen_endpoint.csv", "n excl GCM05", "patients", "patients")
K$ext4_mean <- N("ext4_mean", v4$mean[1], "04_results/validation/GSE237876_frozen_endpoint.csv", "mean excl GCM05")
K$ext4_ci <- N("ext4_ci", ci(v4), "04_results/validation/GSE237876_frozen_endpoint.csv", "ci excl GCM05")

pe <- pur[axis == "epithelial_fraction"]; pm_ <- pur[axis == "malignant_fraction"]
K$slope_epi <- N("slope_epi", pe$slope_per_unit_fraction[1], "04_results/sensitivity/purity_equivalence.csv", "slope epithelial")
K$implied_epi <- N("implied_epi", pe$implied_fraction_difference[1], "04_results/sensitivity/purity_equivalence.csv", "implied epithelial")
K$implied_epi_ci <- N("implied_epi_ci", sprintf("%s to %s", fmt3(min(pe$implied_ci_low[1], pe$implied_ci_high[1])), fmt3(max(pe$implied_ci_low[1], pe$implied_ci_high[1]))), "04_results/sensitivity/purity_equivalence.csv", "implied ci")
K$slope_mal <- N("slope_mal", pm_$slope_per_unit_fraction[1], "04_results/sensitivity/purity_equivalence.csv", "slope malignant")
K$implied_mal <- N("implied_mal", pm_$implied_fraction_difference[1], "04_results/sensitivity/purity_equivalence.csv", "implied malignant")
K$cons_median <- N("cons_median", median(cons$implied_epithelial_fraction_difference), "04_results/sensitivity/single_composition_consistency.csv", "median implied")
K$cons_q1 <- N("cons_q1", quantile(cons$implied_epithelial_fraction_difference, 0.25), "04_results/sensitivity/single_composition_consistency.csv", "q1 implied")
K$cons_q3 <- N("cons_q3", quantile(cons$implied_epithelial_fraction_difference, 0.75), "04_results/sensitivity/single_composition_consistency.csv", "q3 implied")
K$cons_n <- N("cons_n", nrow(cons), "04_results/sensitivity/single_composition_consistency.csv", "n programs", "programs", "programs")
K$n_mix_donors <- N("n_mix_donors", slp[axis == "epithelial_fraction" & program_id == "HRC_CORE", n_donors], "04_results/sensitivity/program_score_vs_malignant_fraction.csv", "n_donors", "donors", "donors")

l183 <- lin[dataset_id == "GSE183904"]
K$lin_epi_min <- N("lin_epi_min", min(l183[grepl("epithelial|malignant", cell_class), median_HRC]), "04_results/cell_source/calibration_lineage_profile.csv", "epithelial min")
K$lin_epi_max <- N("lin_epi_max", max(l183[grepl("epithelial|malignant", cell_class), median_HRC]), "04_results/cell_source/calibration_lineage_profile.csv", "epithelial max")
K$lin_non_min <- N("lin_non_min", min(l183[!grepl("epithelial|malignant", cell_class), median_HRC]), "04_results/cell_source/calibration_lineage_profile.csv", "non-epithelial min")
K$lin_non_max <- N("lin_non_max", max(l183[!grepl("epithelial|malignant", cell_class), median_HRC]), "04_results/cell_source/calibration_lineage_profile.csv", "non-epithelial max")
mpm <- mal[lesion_role == "solid_peritoneal_metastasis"]
K$mal_pm_cells_max <- N("mal_pm_cells_max", max(mpm$n_cells), "04_results/cell_source/malignant_library_summary.csv", "max malignant cells in PM library", "cells", "cells")
K$mal_pm_libs <- N("mal_pm_libs", sum(mpm$n_libraries), "04_results/cell_source/malignant_library_summary.csv", "PM libraries", "libraries", "libraries")
K$mal_pm_range <- N("mal_pm_range", sprintf("%s to %s", formatC(min(mpm$median_HRC), format = "f", digits = 3), formatC(max(mpm$median_HRC), format = "f", digits = 3)), "04_results/cell_source/malignant_library_summary.csv", "PM malignant score range")
mpr <- mal[lesion_role == "primary_tumor" & !is.na(median_HRC)]
K$mal_pt_range <- N("mal_pt_range", sprintf("%s to %s", formatC(min(mpr$median_HRC), format = "f", digits = 3), formatC(max(mpr$median_HRC), format = "f", digits = 3)), "04_results/cell_source/malignant_library_summary.csv", "primary malignant score range")

K$n_237876_specimens <- N("n_237876_specimens", srcc[dataset_id == "GSE237876", n_specimens], "04_results/audit/source_counts.csv", "n_specimens", "specimens", "specimens")
K$n_237876_patients <- N("n_237876_patients", srcc[dataset_id == "GSE237876", n_patients_identifiable], "04_results/audit/source_counts.csv", "n_patients", "patients", "patients")
K$n_237876_pm <- N("n_237876_pm", srcc[dataset_id == "GSE237876", n_solid_pm], "04_results/audit/source_counts.csv", "n_solid_pm", "specimens", "specimens")
K$ffpe_314812 <- N("ffpe_314812", srcc[dataset_id == "GSE314812", n_ffpe], "04_results/audit/source_counts.csv", "n_ffpe", "specimens", "specimens")
K$mtx_genes_314812 <- N("mtx_genes_314812", mtx[source == "GSE314812", n_features], "04_results/audit/matrix_type_audit.csv", "n_features", "genes", "genes")

numbers_f <- pm_out("04_results", "manuscript", "key_numbers.json"); write_json_atomic(K, numbers_f)
trace_dt <- rbindlist(trace)
trace_f <- pm_out("04_results", "manuscript", "number_trace.csv"); write_csv_atomic(trace_dt, trace_f)

# ---- fill templates -------------------------------------------------------------------
fill <- function(path) {
  txt <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  for (k in names(K)) txt <- gsub(paste0("{{", k, "}}"), K[[k]], txt, fixed = TRUE)
  left <- regmatches(txt, gregexpr("\\{\\{[a-z0-9_]+\\}\\}", txt))[[1]]
  if (length(left)) stop("unfilled placeholders: ", paste(unique(left), collapse = ", "))
  txt
}
man_md <- pm_out("06_manuscript", "journal_submission", "manuscript_en.md")
write_text_atomic(strsplit(fill(tmpl_f), "\n")[[1]], man_md)
abs_zh <- pm_out("06_manuscript", "journal_submission", "abstract_zh.md")
write_text_atomic(strsplit(fill(abs_zh_f), "\n")[[1]], abs_zh)
supp_md <- pm_out("06_manuscript", "journal_submission", "supplementary_methods_and_legends.md")
write_text_atomic(strsplit(fill(supp_f), "\n")[[1]], supp_md)
cover_md <- pm_out("06_manuscript", "journal_submission", "cover_letter_en.md")
write_text_atomic(strsplit(fill(cover_f), "\n")[[1]], cover_md)

# abstract_en extracted from the manuscript's abstract section
ml <- readLines(man_md, warn = FALSE)
a0 <- grep("^## Abstract", ml)[1]; a1 <- grep("^## Introduction", ml)[1]
abs_en <- pm_out("06_manuscript", "journal_submission", "abstract_en.md")
write_text_atomic(ml[a0:(a1 - 1)], abs_en)

wc <- function(f) {
  l <- readLines(f, warn = FALSE)
  body <- l[!grepl("^\\s*[|!\\[]", l)]
  length(unlist(strsplit(paste(body, collapse = " "), "\\s+")))
}
cat(sprintf("manuscript words (approx): %d\nabstract words: %d\n", wc(man_md), wc(abs_en)))

# ---- render docx / pdf -----------------------------------------------------------------
pandoc <- cfg$tools$pandoc; soffice <- cfg$tools$soffice
rendered <- character()
if (file.exists(pandoc)) {
  bib <- pm_in("07_references", "references.bib")
  csl <- pm_src("06_manuscript", "source", "vancouver.csl")
  for (src in c(man_md, abs_en, abs_zh, supp_md, cover_md)) {
    docx <- sub("\\.md$", ".docx", src)
    args <- c(shQuote(src), "-o", shQuote(docx), "--from", "markdown+pipe_tables")
    if (grepl("manuscript_en", src) && file.exists(bib)) {
      args <- c(args, "--citeproc", "--bibliography", shQuote(bib))
      if (file.exists(csl)) args <- c(args, "--csl", shQuote(csl))
    }
    system2(pandoc, args, stdout = TRUE, stderr = TRUE)
    if (file.exists(docx)) rendered <- c(rendered, docx)
  }
}
if (file.exists(soffice) && length(rendered)) {
  outdir <- dirname(rendered[1])
  system2(soffice, c("--headless", "--convert-to", "pdf", "--outdir", shQuote(outdir), shQuote(rendered)),
          stdout = TRUE, stderr = TRUE)
  rendered <- c(rendered, sub("\\.docx$", ".pdf", rendered))
}
rendered <- rendered[file.exists(rendered)]
cat("rendered:\n"); print(basename(rendered))
step_end(outputs = c(numbers_f, trace_f, man_md, abs_en, abs_zh, supp_md, cover_md, rendered))
