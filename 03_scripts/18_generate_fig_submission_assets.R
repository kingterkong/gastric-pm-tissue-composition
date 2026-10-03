#!/usr/bin/env Rscript
# Build the strengthened robustness figure, a supplementary adjustment figure, and a
# supplementary workbook containing the post-hoc analyses. Existing frozen-analysis
# figures and result files are not modified.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(openxlsx) })
cfg <- pm_config()
R <- function(...) pm_in("04_results", ...)
inputs <- c(R("sensitivity", "primary_sensitivity.csv"),
            R("sensitivity", "scoring_method_robustness.csv"),
            R("sensitivity", "independent_epithelial_marker_check.csv"),
            R("validation", "GSE237876_published_purity_association.csv"),
            R("validation", "GSE237876_purity_plot_values.csv"),
            R("manuscript", "cross_source_summary.csv"),
            R("primary", "adjusted_models.csv"),
            pm_in("06_manuscript", "journal_submission", "supplementary_tables.xlsx"))
step_begin("18_generate_fig_submission_assets", inputs = inputs,
           seed = cfg$project$seed + 18, uses_protocol = FALSE)

COL_PRIMARY <- "#0072B2"; COL_PM <- "#D55E00"; COL_THIRD <- "#009E73"
INK <- "#1a1a1a"; INK2 <- "#4d4d4d"; GRID <- "#e2e2e0"
theme_pm <- function(base = 8.3) {
  theme_minimal(base_size = base, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = GRID, linewidth = 0.25),
          axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK),
          plot.title = element_text(colour = INK, face = "bold", size = base + 1),
          plot.subtitle = element_text(colour = INK2, size = base - 0.4),
          plot.caption = element_text(colour = INK2, size = base - 1.4, hjust = 0),
          legend.position = "top", legend.title = element_blank(),
          legend.text = element_text(colour = INK),
          strip.text = element_text(colour = INK, face = "bold"))
}
sav <- function(p, directory, name, w, h) {
  pdf_f <- pm_out("05_figures", directory, paste0(name, ".pdf"))
  png_f <- pm_out("05_figures", directory, paste0(name, ".png"))
  tif_f <- pm_out("05_figures", directory, paste0(name, ".tiff"))
  ggsave(pdf_f, p, width = w, height = h, units = "in", device = "pdf")
  ggsave(png_f, p, width = w, height = h, units = "in", dpi = 300)
  ggsave(tif_f, p, width = w, height = h, units = "in", dpi = 600,
         device = "tiff", compression = "lzw")
  c(pdf_f, png_f, tif_f)
}

sens <- fread(inputs[1])
alg <- rbind(fread(inputs[2]), fread(inputs[3]), fill = TRUE)
pur <- fread(inputs[4])
pv <- fread(inputs[5])
cross <- fread(inputs[6])
adj <- fread(inputs[7])

# a: analyses specified before seeing the primary result
sens[, label := fcase(
  analysis == "primary (all eligible patients)", "Primary analysis",
  analysis == "same-preservation pairs only", "Same preservation",
  analysis == "epithelium-restricted program subset", "Epithelium-restricted HRC subset",
  analysis == "alternative source definition (epiHR)", "Alternative programme definition",
  analysis == "GC05 post-treatment pair instead of pre-treatment", "Alternative eligible pair for GC05",
  default = analysis)]
sens[, label := factor(label, levels = rev(label))]
p5a <- ggplot(sens, aes(mean, label)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = INK2, linewidth = 0.55) +
  geom_point(aes(size = n), colour = COL_PM) +
  scale_size_continuous(range = c(1.3, 2.6), guide = "none") +
  labs(title = "a  Pre-specified sensitivity analyses",
       subtitle = "Mean paired difference (solid PM minus primary)", x = "score difference", y = NULL) +
  theme_pm() + theme(panel.grid.major.y = element_blank())

# b: algorithms have different units, so compare their paired standardized effects.
alg[, short := fcase(
  method == "HRC: frozen mean-rank score", "HRC: frozen mean rank",
  method == "HRC: singscore", "HRC: singscore",
  method == "HRC: mean log-CPM", "HRC: mean log-CPM",
  method == "HRC: mean gene-wise z score", "HRC: mean gene z score",
  method == "Independent epithelial markers: mean-rank score", "Epithelial markers: mean rank",
  method == "Independent epithelial markers: mean log-CPM", "Epithelial markers: mean log-CPM",
  default = method)]
alg[, set := ifelse(grepl("Independent", method), "Epithelial markers", "HRC programme")]
alg[, short := factor(short, levels = rev(short))]
p5b <- ggplot(alg, aes(cohen_dz, short, colour = set)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_linerange(aes(xmin = dz_ci_low, xmax = dz_ci_high), linewidth = 0.6) +
  geom_point(size = 2) +
  scale_colour_manual(values = c(`HRC programme` = COL_PM,
                                 `Epithelial markers` = COL_PRIMARY)) +
  labs(title = "b  Scoring and marker robustness",
       subtitle = "Paired standardized effect; 95% bootstrap CI",
       x = "Cohen's dz", y = NULL) +
  theme_pm() + theme(panel.grid.major.y = element_blank())

# c: published histological purity, centered within each patient.
fe <- pur[analysis == "all specimens: patient fixed-effects slope, cluster bootstrap CI"]
pv[, site_group := fcase(lesion_role == "primary_tumor", "primary tumour",
                         lesion_role == "solid_peritoneal_metastasis", "solid PM",
                         default = "other metastasis")]
p5c <- ggplot(pv, aes(purity_centered, hrc_centered, colour = site_group)) +
  geom_hline(yintercept = 0, colour = GRID, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = GRID, linewidth = 0.3) +
  geom_abline(intercept = 0, slope = fe$estimate, colour = INK, linewidth = 0.65) +
  geom_point(size = 1.45, alpha = 0.80) +
  scale_colour_manual(values = c(`primary tumour` = COL_PRIMARY, `solid PM` = COL_PM,
                                 `other metastasis` = "#7A7A78")) +
  labs(title = "c  Score tracks published histological purity",
       subtitle = sprintf("66 specimens, 14 patients\nwithin-patient slope %.3f (cluster-bootstrap 95%% CI %.3f to %.3f)",
                          fe$estimate, fe$ci_low, fe$ci_high),
       x = "purity deviation from patient mean", y = "HRC-score deviation from patient mean") +
  theme_pm()

# d: same frozen endpoint in both sources.
cross[, source2 := fcase(source == "GSE314812 (primary test)", "Primary cohort (n=17)",
                         source == "GSE237876", "External cohort (n=5)",
                         default = "External, disputed site excluded (n=4)")]
cross[, source2 := factor(source2, levels = rev(source2))]
p5d <- ggplot(cross, aes(mean, source2)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = INK2, linewidth = 0.55) +
  geom_point(size = 2.1, colour = COL_PM) +
  labs(title = "d  Independent paired cohort",
       subtitle = "Identical frozen endpoint\nlimited precision with five patients",
       x = "mean paired score difference", y = NULL) +
  theme_pm() + theme(panel.grid.major.y = element_blank())

fig5 <- (p5a | p5b) / (p5c | p5d) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(title = "The composition result is robust to scoring choices and external purity data",
                  theme = theme_pm())
outs <- sav(fig5, "main", "Figure5_robustness_external_purity", 7.6, 7.1)

# Supplementary Figure S1 retains every one-covariate adjustment from the frozen plan.
a <- adj[term == "(Intercept)" & !is.na(estimate)]
a[, label := fcase(
  model == "unadjusted (primary)", "Unadjusted",
  model == "preservation change", "Preservation change",
  grepl("adipocyte", model), "Adipocyte compartment",
  model == "fibroblast compartment delta", "Fibroblast compartment",
  model == "stromal abundance delta", "Stromal abundance",
  model == "proliferation delta", "Proliferation",
  default = model)]
a[, label := factor(label, levels = rev(label))]
ps1 <- ggplot(a, aes(estimate, label)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = INK2, linewidth = 0.65) +
  geom_point(size = 2.2, colour = COL_PM) +
  labs(title = "Adjustment for a single composition or technical covariate",
       subtitle = "Model intercept: HRC difference remaining after adjustment; 17 paired patients",
       x = "adjusted mean difference", y = NULL,
       caption = "Each model adds one term, except the categorical preservation-change model.") +
  theme_pm(9) + theme(panel.grid.major.y = element_blank())
outs <- c(outs, sav(ps1, "supplementary", "Supplementary_Figure_S1_adjustment", 6.8, 4.4))

# Add all post-hoc numerical results to a new copy of the supplementary workbook.
src_wb <- inputs[8]
if (file.exists(src_wb)) {
  wb <- loadWorkbook(src_wb)
} else {
  wb <- createWorkbook()
}
for (nm in c("S11_scoring_robustness", "S12_epithelial_markers", "S13_external_purity", "S14_paired_purity")) {
  if (nm %in% names(wb)) removeWorksheet(wb, nm)
  addWorksheet(wb, nm)
}
writeData(wb, "S11_scoring_robustness", fread(inputs[2]))
writeData(wb, "S12_epithelial_markers", fread(inputs[3]))
writeData(wb, "S13_external_purity", fread(inputs[4]))
writeData(wb, "S14_paired_purity", fread(R("validation", "GSE237876_paired_purity_deltas.csv")))
o_wb <- pm_out("04_results", "manuscript", "supplementary_tables_FIG.xlsx")
assert_writable(o_wb); saveWorkbook(wb, o_wb, overwrite = TRUE)
outs <- c(outs, o_wb)

step_end(outputs = outs)
