# Stage 12: main figures and tables, generated only from the result files. Colours are the
# validated colourblind-safe pair/triplet; every panel names the unit (patients, libraries,
# cells) and no number is re-computed here.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(openxlsx) })
cfg <- pm_config()

R <- function(...) pm_in("04_results", ...)
inputs <- c(R("primary", "primary_endpoint.csv"), R("primary", "secondary_families.csv"),
            R("primary", "GSE314812_patient_deltas.csv"), R("primary", "GSE314812_program_scores_by_specimen.csv"),
            R("primary", "adjusted_models.csv"), R("sensitivity", "primary_sensitivity.csv"),
            R("sensitivity", "leave_one_patient_out.csv"), R("sensitivity", "single_composition_consistency.csv"),
            R("sensitivity", "program_score_vs_malignant_fraction.csv"), R("sensitivity", "pseudobulk_mixture_scores.csv"),
            R("validation", "GSE237876_frozen_endpoint.csv"), R("cell_source", "host_compartment_deltas.csv"),
            R("cell_source", "calibration_lineage_profile.csv"), R("cell_source", "malignant_library_summary.csv"),
            R("manuscript", "cross_source_summary.csv"), pm_in("00_admin", "eligibility_flow.csv"),
            pm_in("04_results", "audit", "source_counts.csv"))
inputs <- inputs[file.exists(inputs)]
step_begin("12_generate_figures_tables", inputs = inputs, seed = cfg$project$seed + 12, uses_protocol = TRUE)

# ---- palette (validated: see 00_admin/figure_palette_validation.txt) -----------------
COL_PRIMARY <- "#0072B2"; COL_PM <- "#D55E00"; COL_THIRD <- "#009E73"
INK <- "#1a1a1a"; INK2 <- "#4d4d4d"; GRID <- "#e2e2e0"
theme_pm <- function(base = 8) {
  theme_minimal(base_size = base, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = GRID, linewidth = 0.25),
          axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK),
          plot.title = element_text(colour = INK, face = "bold", size = base + 1),
          plot.subtitle = element_text(colour = INK2, size = base - 0.5),
          plot.caption = element_text(colour = INK2, size = base - 1.5, hjust = 0),
          legend.position = "top", legend.title = element_blank(), legend.text = element_text(colour = INK),
          strip.text = element_text(colour = INK, face = "bold"))
}
sav <- function(p, name, w, h) {
  pdf_f <- pm_out("05_figures", "main", paste0(name, ".pdf"))
  png_f <- pm_out("05_figures", "main", paste0(name, ".png"))
  ggsave(pdf_f, p, width = w, height = h, units = "in", device = "pdf")
  ggsave(png_f, p, width = w, height = h, units = "in", dpi = 300)
  c(pdf_f, png_f)
}
outs <- character()

D <- fread(R("primary", "GSE314812_patient_deltas.csv"))
sc <- fread(R("primary", "GSE314812_program_scores_by_specimen.csv"))
prim <- fread(R("primary", "primary_endpoint.csv"))
sec <- fread(R("primary", "secondary_families.csv"))
adj <- fread(R("primary", "adjusted_models.csv"))
sens <- fread(R("sensitivity", "primary_sensitivity.csv"))
loo <- fread(R("sensitivity", "leave_one_patient_out.csv"))
cons <- fread(R("sensitivity", "single_composition_consistency.csv"))
mixsc <- fread(R("sensitivity", "pseudobulk_mixture_scores.csv"))
val <- fread(R("validation", "GSE237876_frozen_endpoint.csv"))
adm <- fread(R("cell_source", "host_compartment_deltas.csv"))
lin <- fread(R("cell_source", "calibration_lineage_profile.csv"))
mal <- fread(R("cell_source", "malignant_library_summary.csv"))
flow <- fread(pm_in("00_admin", "eligibility_flow.csv"))
meso_cov <- fread(R("primary", "GSE314812_program_coverage.csv"))[program_id == "COMPART_MESOTHELIUM"]
srcc <- fread(pm_in("04_results", "audit", "source_counts.csv"))

pretty <- function(x) {
  x <- gsub("^COMPART_", "", x); x <- gsub("^CTRL_", "", x); x <- gsub("^SENS_", "", x); x <- gsub("^NEG_", "", x)
  x <- gsub("_", " ", tolower(x))
  x <- gsub("hrc core", "relapse program (primary endpoint)", x)
  x <- gsub("hrc epi restricted", "relapse program, epithelium-restricted", x)
  x <- gsub("hrc epihr", "relapse program, alternative definition", x)
  x <- gsub("yap stromal", "stroma-restricted YAP targets (negative control)", x)
  x <- gsub("ieg stress", "immediate-early stress", x)
  x <- gsub("emt", "EMT", x); x <- gsub("g2m", "G2/M", x)
  x
}

# ---- Figure 1: which patients and specimens carry which evidence --------------------
f1a <- flow[dataset_id == "GSE314812"]
f1a[, step := factor(step, levels = rev(step))]
p1a <- ggplot(f1a, aes(n, step)) +
  geom_col(fill = COL_PRIMARY, width = 0.62) +
  geom_text(aes(label = paste0(n, " ", unit, ifelse(n == 1, "", "s"))), hjust = -0.12, size = 2.5, colour = INK) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
  labs(title = "a  Primary cohort (GSE314812)", x = NULL, y = NULL) + theme_pm() +
  theme(panel.grid.major.y = element_blank())
f1b <- flow[dataset_id == "GSE237876"]
f1b[, step := factor(step, levels = rev(step))]
p1b <- ggplot(f1b, aes(n, step)) +
  geom_col(fill = COL_PM, width = 0.62) +
  geom_text(aes(label = paste0(n, " ", unit, ifelse(n == 1, "", "s"))), hjust = -0.12, size = 2.5, colour = INK) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.32))) +
  labs(title = "b  Independent paired cohort (GSE237876)", x = NULL, y = NULL) + theme_pm() +
  theme(panel.grid.major.y = element_blank())
cw_all <- fread(pm_in("00_admin", "sample_crosswalk.csv"))
n_blood <- cw_all[dataset_id != "GSE234129" & specimen_context == "blood", .(blood = .N), by = dataset_id]
sp0 <- srcc[, .(dataset_id, primary = n_primary, `solid PM` = n_solid_pm, ascites = n_ascites,
                lavage = n_lavage, `other metastasis` = n_other_met, normal = n_normal)]
sp0 <- merge(sp0, n_blood, by = "dataset_id", all.x = TRUE)
sp0[is.na(blood), blood := 0]
sp <- melt(sp0, id.vars = "dataset_id", variable.name = "type", value.name = "n")[n > 0]
sp[, type := factor(type, levels = c("primary", "solid PM", "ascites", "lavage", "other metastasis", "normal", "blood"))]
p1c <- ggplot(sp, aes(n, dataset_id, fill = type)) +
  geom_col(width = 0.68) +
  scale_fill_manual(values = c(primary = COL_PRIMARY, `solid PM` = COL_PM, ascites = COL_THIRD,
                               lavage = "#7a7a7a", `other metastasis` = "#9a6fb0", normal = "#c9c9c7",
                               blood = "#8c6d4f"), drop = FALSE) +
  labs(title = "c  Specimens by type in every source used", x = "specimens", y = NULL) + theme_pm() +
  theme(panel.grid.major.y = element_blank())
fig1 <- (p1a / p1b / p1c) + plot_layout(heights = c(1.1, 1.05, 1.5))
outs <- c(outs, sav(fig1, "Figure1_study_design", 7.0, 7.4))

# ---- Figure 2: the primary endpoint --------------------------------------------------
pair_long <- sc[eligible_primary == TRUE, .(canonical_patient_id, lesion_role, HRC_CORE)]
pair_long[, site := ifelse(lesion_role == "primary_tumor", "primary tumour", "solid peritoneal metastasis")]
ord <- D[order(HRC_CORE), canonical_patient_id]
pair_long[, canonical_patient_id := factor(canonical_patient_id, levels = ord)]
p2a <- ggplot(pair_long, aes(site, HRC_CORE, group = canonical_patient_id)) +
  geom_line(colour = "#b0b0ae", linewidth = 0.4) +
  geom_point(aes(colour = site), size = 1.9) +
  scale_colour_manual(values = c(`primary tumour` = COL_PRIMARY, `solid peritoneal metastasis` = COL_PM)) +
  labs(title = "a  Paired specimens", subtitle = sprintf("%d patients, one pair each", prim$n[1]),
       y = "relapse-program score", x = NULL) + theme_pm() + theme(legend.position = "none")
D2 <- copy(D); D2[, canonical_patient_id := factor(canonical_patient_id, levels = ord)]
D2[, dir := ifelse(HRC_CORE < 0, "lower in metastasis", "higher in metastasis")]
p2b <- ggplot(D2, aes(HRC_CORE, canonical_patient_id, fill = dir)) +
  geom_col(width = 0.66) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  scale_fill_manual(values = c(`lower in metastasis` = COL_PM, `higher in metastasis` = COL_PRIMARY)) +
  labs(title = "b  Within-patient difference", x = "metastasis − primary", y = NULL) + theme_pm() +
  theme(panel.grid.major.y = element_blank())
p2c <- ggplot(prim, aes(x = mean, y = "mean")) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = COL_PM, linewidth = 0.7) +
  geom_point(size = 2.6, colour = COL_PM) +
  geom_text(aes(label = sprintf("%.3f (95%% CI %.3f to %.3f), P = %.4f, %d of %d patients lower",
                                mean, ci_low, ci_high, p_value, n_negative, n)),
            vjust = -1.5, size = 2.5, colour = INK) +
  scale_y_discrete(expand = expansion(add = c(0.4, 0.9))) +
  labs(title = "c  Primary endpoint", x = "mean within-patient difference", y = NULL) + theme_pm()
fig2 <- (p2a | p2b) / p2c + plot_layout(heights = c(2.6, 1))
outs <- c(outs, sav(fig2, "Figure2_primary_endpoint", 7.0, 5.6))

# ---- Figure 3: what else moved, and in which direction ------------------------------
allp <- rbind(prim[, .(program_id, mean, ci_low, ci_high, family = "primary endpoint")],
              sec[!is.na(mean), .(program_id, mean, ci_low, ci_high, family)], fill = TRUE)
allp[, grp := fcase(family == "primary endpoint", "Primary endpoint",
                    grepl("secondary_A", family), "Control and negative-control programs",
                    grepl("secondary_B", family), "Host-tissue compartments",
                    default = "Alternative definitions of the primary program")]
allp[, grp := factor(grp, levels = c("Primary endpoint", "Alternative definitions of the primary program",
                                     "Host-tissue compartments", "Control and negative-control programs"))]
allp[, label := pretty(program_id)]
allp[, label := factor(label, levels = allp[order(mean), label])]
allp[, dirc := ifelse(mean < 0, "lower in metastasis", "higher in metastasis")]
rng <- range(c(allp$ci_low, allp$ci_high), na.rm = TRUE)
panel_of <- function(g, show_x) {
  d <- allp[grp == g]
  ggplot(d, aes(mean, label, colour = dirc)) +
    geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
    geom_linerange(aes(xmin = ci_low, xmax = ci_high), linewidth = 0.6) +
    geom_point(size = 1.9) +
    scale_colour_manual(values = c(`higher in metastasis` = COL_PRIMARY, `lower in metastasis` = COL_PM),
                        limits = c("higher in metastasis", "lower in metastasis"), drop = FALSE) +
    scale_x_continuous(limits = rng) +
    labs(title = g, x = if (show_x) "mean within-patient difference" else NULL, y = NULL) +
    theme_pm() +
    theme(panel.grid.major.y = element_blank(), plot.title = element_text(size = 8, face = "bold"),
          legend.position = "top",
          axis.text.x = if (show_x) element_text() else element_blank())
}
lv <- levels(allp$grp)
panels <- lapply(seq_along(lv), function(i) panel_of(lv[i], i == length(lv)))
p3 <- Reduce(`/`, panels) +
  plot_layout(heights = allp[, .N, by = grp][match(lv, grp), N], guides = "collect") +
  plot_annotation(title = "Every program moves as tissue composition predicts",
                  subtitle = sprintf("mean within-patient difference, %d patients; bars are 95%% confidence intervals", prim$n[1]),
                  caption = sprintf("Mesothelium is absent: it did not meet the frozen coverage rule on this platform (%d of %d genes measured).",
                         meso_cov$n_used[1], meso_cov$n_defined[1]),
                  theme = theme_pm())
outs <- c(outs, sav(p3, "Figure3_program_and_compartment_changes", 6.8, 5.6))

# ---- Figure 4: the program is an epithelial-content readout --------------------------
lin4 <- lin[dataset_id == "GSE183904"][order(median_HRC)]
nice_class <- function(x) {
  x <- gsub("malignant_cnv_high", "malignant (CNV high)", x)
  x <- gsub("tumour_epithelial_cnv_low_or_na", "tumour epithelial (CNV low or undetermined)", x)
  x <- gsub("epithelial_normal_sample", "epithelial (normal-tissue sample)", x)
  x <- gsub("^t_nk$", "T and NK", x); x <- gsub("^b_cell$", "B cell", x)
  gsub("_", " ", x)
}
lin4[, cell_class := factor(nice_class(cell_class), levels = nice_class(cell_class))]
lin4[, kind := ifelse(grepl("malignant|epithelial|tumour", cell_class), "epithelial", "non-epithelial")]
p4a <- ggplot(lin4, aes(median_HRC, cell_class, fill = kind)) +
  geom_col(width = 0.66) +
  scale_fill_manual(values = c(epithelial = COL_PRIMARY, `non-epithelial` = "#9a9a98")) +
  labs(title = "a  The program scores epithelial identity", x = "median score (calibration donors)", y = NULL,
       subtitle = "GSE183904 gastric primary tumours, development side") +
  theme_pm() + theme(panel.grid.major.y = element_blank())
mx <- mixsc[axis == "epithelial_fraction", .(score = mean(HRC_CORE)), by = .(sample_id, fraction)]
p4b <- ggplot(mx, aes(fraction, score, group = sample_id)) +
  geom_line(colour = "#b0b0ae", linewidth = 0.4) + geom_point(colour = COL_PRIMARY, size = 1.1) +
  labs(title = "b  Known mixtures track epithelial content",
       subtitle = "300-cell pseudo-bulk mixtures, 8 calibration donors",
       x = "epithelial-cell fraction (known)", y = "relapse-program score") + theme_pm()
cons4 <- copy(cons)
cons4[, label := pretty(program_id)]
cons4[, label := factor(label, levels = cons4[order(implied_epithelial_fraction_difference), label])]
p4c <- ggplot(cons4, aes(implied_epithelial_fraction_difference, label)) +
  geom_vline(xintercept = median(cons$implied_epithelial_fraction_difference), colour = COL_PM,
             linetype = "22", linewidth = 0.5) +
  geom_linerange(aes(xmin = implied_ci_low, xmax = implied_ci_high), colour = INK2, linewidth = 0.5) +
  geom_point(size = 1.8, colour = COL_PRIMARY) +
  labs(title = "c  One compositional shift explains every program",
       subtitle = sprintf("epithelial-content difference each program implies; dashed line is the median (%.2f)",
                          median(cons$implied_epithelial_fraction_difference)),
       x = "implied difference in epithelial fraction", y = NULL,
       caption = paste("Programmes cluster near a single value. Gastric-mucosa markers are the expected outlier:\n",
                       "peritoneal tissue contains none, so they shift further than dilution alone predicts.\n",
                       "The relationship in b saturates at high epithelial content, so these conversions are approximations.")) +
  theme_pm() + theme(panel.grid.major.y = element_blank())
fig4 <- (p4a | p4b) / p4c + plot_layout(heights = c(1.5, 1.6))
outs <- c(outs, sav(fig4, "Figure4_composition_explains_the_change", 7.6, 6.8))

# ---- Figure 5: robustness, adjustment and the independent cohort --------------------
s5 <- copy(sens)[!is.na(mean)]
s5[, analysis := factor(analysis, levels = rev(analysis))]
p5a <- ggplot(s5, aes(mean, analysis)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = INK2, linewidth = 0.55) +
  geom_point(aes(size = n), colour = COL_PM) + scale_size_continuous(range = c(1.2, 2.6), guide = "none") +
  labs(title = "a  Sensitivity analyses", x = "mean within-patient difference", y = NULL,
       subtitle = "point size is the number of patients") +
  theme_pm() + theme(panel.grid.major.y = element_blank())
a5 <- adj[grepl("Intercept", term) | is.na(term)][!is.na(estimate)]
a5[, model := factor(model, levels = rev(unique(model)))]
p5b <- ggplot(a5, aes(estimate, model)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = INK2, linewidth = 0.55) +
  geom_point(size = 2, colour = COL_PM) +
  labs(title = "b  Difference remaining after adjustment", x = "adjusted mean difference (model intercept)", y = NULL,
       subtitle = "each model adds at most one composition term") +
  theme_pm() + theme(panel.grid.major.y = element_blank())
v5 <- fread(R("manuscript", "cross_source_summary.csv"))
v5[, source := factor(source, levels = rev(source))]
lab_x <- max(v5$ci_high, na.rm = TRUE) + 0.035
p5c <- ggplot(v5, aes(mean, source)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = INK2, linewidth = 0.55) +
  geom_point(size = 2.2, colour = COL_PM) +
  geom_text(aes(x = lab_x, label = sprintf("%d patients", n)), hjust = 0, size = 2.4, colour = INK2) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.34))) +
  labs(title = "c  The same frozen endpoint across sources", x = "mean within-patient difference", y = NULL,
       caption = paste("With five patients the smallest attainable two-sided exact P is 0.0625,\n",
                       "so the independent cohort can show direction only.")) +
  theme_pm() + theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 6.4))
fig5 <- (p5a / p5b / p5c) + plot_layout(heights = c(1.2, 1.1, 1.1))
outs <- c(outs, sav(fig5, "Figure5_robustness_and_independent_source", 7.2, 6.8))

# ---- Figure 6: what the available single-cell data can and cannot resolve -----------
m6 <- copy(mal)[!is.na(median_HRC)]
m6[, site := ifelse(lesion_role == "primary_tumor", "primary tumour", "solid peritoneal metastasis")]
bulk_rng <- range(pair_long$HRC_CORE, na.rm = TRUE)
p6a <- ggplot(m6, aes(median_HRC, dataset_id, colour = site)) +
  annotate("rect", xmin = bulk_rng[1], xmax = bulk_rng[2], ymin = -Inf, ymax = Inf,
           fill = "#f0f0ee", colour = NA) +
  geom_point(aes(size = n_cells), alpha = 0.9) +
  scale_colour_manual(values = c(`primary tumour` = COL_PRIMARY, `solid peritoneal metastasis` = COL_PM)) +
  scale_size_continuous(range = c(1.4, 5), name = "malignant cells") +
  scale_x_continuous(limits = c(min(bulk_rng[1], min(m6$median_HRC)) - 0.01, max(bulk_rng[2], max(m6$median_HRC)) + 0.01)) +
  labs(title = "a  Program value inside malignant cells", x = "median score across libraries", y = NULL,
       subtitle = paste0("libraries reaching 20 malignant cells; donors unverified in GSE308231.\n",
                         "Shaded band is the full range of paired tissue-level scores in the primary cohort.")) +
  theme_pm() + theme(panel.grid.major.y = element_blank())
cnt6 <- m6[, .(cells = sum(n_cells), libs = .N), by = .(dataset_id, site)]
p6b <- ggplot(cnt6, aes(cells, dataset_id, fill = site)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.68) +
  geom_text(aes(label = sprintf("%d cells / %d lib", cells, libs)), position = position_dodge(width = 0.72),
            hjust = -0.08, size = 2.2, colour = INK) +
  scale_fill_manual(values = c(`primary tumour` = COL_PRIMARY, `solid peritoneal metastasis` = COL_PM)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.45))) +
  labs(title = "b  How little malignant material the peritoneal libraries contain",
       x = "malignant cells passing quality control", y = NULL,
       caption = "This is why the cell-intrinsic question cannot be settled with the public data available.") +
  theme_pm() + theme(panel.grid.major.y = element_blank())
fig6 <- (p6a / p6b) + plot_layout(heights = c(1, 1))
outs <- c(outs, sav(fig6, "Figure6_cell_level_limits", 7.0, 5.0))

# ---- main tables ----------------------------------------------------------------------
t1 <- srcc[, .(Source = dataset_id, Specimens = n_specimens, `Patients (identifiable)` = n_patients_identifiable,
               Primary = n_primary, `Solid PM` = n_solid_pm, Ascites = n_ascites, Lavage = n_lavage,
               `Other metastasis` = n_other_met, Normal = n_normal,
               `Patients with primary + solid PM` = n_patients_primary_and_solid_pm,
               `Therapy known (patients)` = n_therapy_known, FFPE = n_ffpe, `Fresh-frozen` = n_fresh_frozen)]
t2 <- allp[, .(Program = label, Role = grp, `Mean difference` = round(mean, 4),
               `95% CI low` = round(ci_low, 4), `95% CI high` = round(ci_high, 4))]
pq <- rbind(sec[, .(Program = pretty(program_id), P = p_value, q = q_value)],
            prim[, .(Program = pretty(program_id), P = p_value, q = NA_real_)], fill = TRUE)
t2 <- merge(t2, unique(pq, by = "Program"), by = "Program", all.x = TRUE)
setorder(t2, Role, -`Mean difference`)
t3 <- cons[, .(Program = pretty(program_id), `Observed difference` = round(mean, 4),
               `Slope vs epithelial fraction` = round(slope, 3),
               `Implied epithelial-fraction difference` = round(implied_epithelial_fraction_difference, 3),
               `Implied CI low` = round(implied_ci_low, 3), `Implied CI high` = round(implied_ci_high, 3))]
wb <- createWorkbook()
for (nm in c("Table1_sources", "Table2_program_effects", "Table3_composition_equivalence")) addWorksheet(wb, nm)
writeData(wb, "Table1_sources", t1); writeData(wb, "Table2_program_effects", t2); writeData(wb, "Table3_composition_equivalence", t3)
o_tab <- pm_out("04_results", "manuscript", "main_tables.xlsx")
assert_writable(o_tab); saveWorkbook(wb, o_tab, overwrite = TRUE)
for (x in list(list(t1, "Table1_sources"), list(t2, "Table2_program_effects"), list(t3, "Table3_composition_equivalence")))
  write_csv_atomic(x[[1]], pm_out("04_results", "manuscript", paste0(x[[2]], ".csv")))
outs <- c(outs, o_tab, pm_in("04_results", "manuscript", paste0(c("Table1_sources", "Table2_program_effects", "Table3_composition_equivalence"), ".csv")))

cat("figures written:\n"); print(basename(grep("\\.pdf$", outs, value = TRUE)))
step_end(outputs = outs)
