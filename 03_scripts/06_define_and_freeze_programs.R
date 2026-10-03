# Stage 6: derive every program from its primary source file, harmonise symbols to HGNC,
# characterise cell-type specificity from an independent reference (HPA v25), build the
# host-compartment covariate sets, compute per-platform coverage, and freeze everything
# with hashes. Nothing here looks at any target-cohort site contrast.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
source(pm_src("03_scripts", "R", "sc_utils.R"))
suppressPackageStartupMessages({ library(readxl); library(jsonlite) })
cfg <- pm_config()

P <- function(...) pm_raw("other_sources", "annotation", ...)
S <- function(...) pm_raw("literature", "programs", ...)
src <- c(
  hrc      = S("alt", "canellas2022_hrc", "CanellasSocias2022_Nature_MOESM4_SuppTable2_EpiHR_coreHRC.xlsx"),
  hallmark = S("alt", "msigdb_2026.1.Hs", "h.all.v2026.1.Hs.json"),
  tca      = S("alt", "msigdb_2026.1.Hs", "c4.3ca.v2026.1.Hs.json"),
  estimate = S("alt", "estimate_yoshihara2013", "Yoshihara2013_NatCommun_MOESM488_SuppData1_stromal_immune_signatures.xlsx"),
  yap      = S("yap_ecm", "genesets", "yap_ecm_programs_v1.json"),
  hpa      = P("hpa_rna_single_cell_type.tsv.zip"),
  hgnc     = P("hgnc_complete_set_20260925.txt")
)
ann_f <- pm_in("02_data_processed", "bulk", "annotation_grch37_r75_hgnc.csv")
step_begin("06_define_and_freeze_programs", inputs = c(src, ann_f), seed = cfg$project$seed + 6)

maps <- hgnc_maps(src[["hgnc"]])
harm <- function(g) {
  g <- unique(trimws(g)); g <- g[nzchar(g) & !is.na(g)]
  h <- harmonise_symbols(g, NULL, maps)
  h[, .(gene_original = original, gene_mapped = symbol, mapping = how)]
}

# ---- primary program: EMP1+ high-relapse-cell core program (Canellas-Socias 2022) ----
hrc_tab <- as.data.table(read_excel(src[["hrc"]], sheet = 1, skip = 3))
setnames(hrc_tab, 1:6, c("gene", "corr_SMC", "corr_KUL", "epiHR", "cluster", "core_HRC"))
hrc_core <- hrc_tab[as.numeric(core_HRC) == 1, gene]
hrc_epi <- hrc_tab[as.numeric(epiHR) == 1, gene]

# ---- MSigDB / 3CA sets ----
msig_get <- function(f, name) { j <- fromJSON(f); stopifnot(name %in% names(j)); j[[name]]$geneSymbols }
hallmark_emt <- msig_get(src[["hallmark"]], "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
hallmark_g2m <- msig_get(src[["hallmark"]], "HALLMARK_G2M_CHECKPOINT")
tca_json <- fromJSON(src[["tca"]])
stress_nm <- grep("STRESS$", names(tca_json), value = TRUE)
stress_nm <- stress_nm[!grepl("VITRO", stress_nm)][1]
tca_stress <- tca_json[[stress_nm]]$geneSymbols

# ---- ESTIMATE stromal signature ----
# Supplementary Data 1 is wide: one row per signature, genes spread across the columns.
est <- as.data.table(read_excel(src[["estimate"]], sheet = "signatures", col_names = FALSE))
est_row <- est[grepl("^Stromal", trimws(as.character(est[[1]])))]
stopifnot(nrow(est_row) == 1)
estimate_stromal <- unlist(est_row[, 3:ncol(est_row)], use.names = FALSE)
estimate_stromal <- trimws(estimate_stromal[!is.na(estimate_stromal) & nzchar(trimws(estimate_stromal))])

# ---- YAP consensus target set (for the pre-specified negative control) ----
yap_json <- fromJSON(src[["yap"]], simplifyVector = FALSE)$programs
yap_key <- grep("CONSENSUS_CORE_DERIVED", names(yap_json), value = TRUE)[1]
stopifnot(!is.na(yap_key))
yap_consensus <- unlist(yap_json[[yap_key]]$genes, use.names = FALSE)
stopifnot(length(yap_consensus) > 10)

# ---- independent cell-type reference: HPA v25 single cell type ----
hpa <- fread(cmd = sprintf("unzip -p %s", shQuote(src[["hpa"]])), col.names = c("ensembl", "gene", "cell_type", "ncpm"))
COMPART <- list(
  gastric_epithelium = c("foveolar cells", "gastric chief cells", "gastric progenitor cells", "mucous neck cells", "parietal cells"),
  intestinal_epithelium = c("enterocytes", "goblet cells", "paneth cells", "intestinal stem cells", "tuft cells"),
  mesothelium = "mesothelial cells",
  adipocyte = "adipocytes",
  fibroblast = c("fibroblasts", "fibro-adipogenic progenitors"),
  smooth_muscle = c("smooth muscle cells", "vascular smooth muscle cells", "pericytes"),
  endothelium = c("vascular endothelial cells", "lymphatic endothelial cells"),
  immune = c("macrophages", "monocytes", "t-cells", "b-cells", "nk-cells", "plasma cells", "granulocytes", "neutrophils",
             "dendritic cells", "cdc", "pdcs", "mast cells", "kupffer cells", "microglia", "langerhans cells")
)
COMPART_requested <- COMPART
COMPART <- lapply(COMPART, function(x) intersect(x, unique(hpa$cell_type)))
missing_ct <- unlist(mapply(function(req, got) setdiff(req, got), COMPART_requested, COMPART))
if (length(missing_ct)) message("HPA cell types not found (ignored): ", paste(missing_ct, collapse = ", "))
stopifnot(all(lengths(COMPART) > 0))
hpa_c <- rbindlist(lapply(names(COMPART), function(k) {
  if (!length(COMPART[[k]])) return(NULL)
  hpa[cell_type %in% COMPART[[k]], .(value = max(ncpm)), by = .(gene)][, compartment := k]
}))
hpa_w <- dcast(hpa_c, gene ~ compartment, value.var = "value", fill = 0)
epi_cols <- c("gastric_epithelium", "intestinal_epithelium")
non_epi <- setdiff(names(hpa_w), c("gene", epi_cols))
hpa_w[, epi_max := do.call(pmax, .SD), .SDcols = epi_cols]
hpa_w[, nonepi_max := do.call(pmax, .SD), .SDcols = non_epi]
hpa_w[, nonepi_which := non_epi[max.col(as.matrix(.SD), ties.method = "first")], .SDcols = non_epi]
hpa_w[, specificity := fifelse(pmax(epi_max, nonepi_max) < 5, "low_expression",
                        fifelse(epi_max >= 2 * nonepi_max, "epithelium_dominant",
                        fifelse(nonepi_max >= 2 * epi_max, paste0("nonepithelial:", nonepi_which), "shared")))]
o_hpa <- pm_out("02_data_processed", "frozen_programs", "hpa_v25_compartment_max_ncpm.csv.gz"); write_csv_atomic(hpa_w, o_hpa)

# ---- host-compartment covariate sets (independent of every program above) ----
prog_raw <- list(HRC_CORE = hrc_core, HRC_EPI = hrc_epi, EMT = hallmark_emt, G2M = hallmark_g2m,
                 STROMA = estimate_stromal, STRESS = tca_stress, YAP = yap_consensus)
prog_mapped <- lapply(prog_raw, function(g) harm(g)$gene_mapped)
excluded_from_compartments <- unique(unlist(prog_mapped))
N_COMPART_GENES <- 40
compartment_sets <- lapply(setdiff(names(COMPART), "intestinal_epithelium"), function(k) {
  d <- hpa_w[specificity != "low_expression" & !gene %in% excluded_from_compartments]
  others <- setdiff(names(COMPART), k)
  d <- d[get(k) >= 10]
  d[, other_max := do.call(pmax, .SD), .SDcols = others]
  d <- d[get(k) >= 3 * other_max][order(-get(k))]
  head(d$gene, N_COMPART_GENES)
})
names(compartment_sets) <- paste0("COMPART_", toupper(setdiff(names(COMPART), "intestinal_epithelium")))

# ---- pre-specified derived subsets ----
hrc_map <- harm(hrc_core)
hrc_spec <- merge(hrc_map, hpa_w[, .(gene = gene, specificity)], by.x = "gene_mapped", by.y = "gene", all.x = TRUE)
hrc_epi_restricted <- hrc_spec[specificity == "epithelium_dominant", gene_mapped]
yap_map <- harm(yap_consensus)
yap_spec <- merge(yap_map, hpa_w[, .(gene = gene, specificity)], by.x = "gene_mapped", by.y = "gene", all.x = TRUE)
yap_stromal <- yap_spec[grepl("^nonepithelial:(fibroblast|mesothelium|endothelium|smooth_muscle|adipocyte)", specificity), gene_mapped]

# ---- assemble frozen program table ----
defs <- list(
  list(id = "HRC_CORE", label = "EMP1+ high-relapse-cell core program (colorectal cancer origin)", role = "primary",
       genes = prog_mapped$HRC_CORE, raw = hrc_core, direction = 1,
       source = "Canellas-Socias et al. 2022 Nature, Supplementary Table 2 (core_HRC == 1)", doi = "10.1038/s41586-022-05402-9", pmid = "36352230",
       version = "Nature ESM MOESM4", overlap = "derived from human colorectal cancer bulk/scRNA and mouse models; no gastric cancer; none of this project's cohorts contributed"),
  list(id = "CTRL_EMT", label = "Pan-EMT (hallmark)", role = "control",
       genes = prog_mapped$EMT, raw = hallmark_emt, direction = 1,
       source = "MSigDB HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", doi = "10.1016/j.cels.2015.12.004", pmid = "26771021",
       version = "MSigDB 2026.1.Hs", overlap = "independent of all cohorts"),
  list(id = "CTRL_PROLIFERATION", label = "Proliferation (G2/M checkpoint)", role = "control",
       genes = prog_mapped$G2M, raw = hallmark_g2m, direction = 1,
       source = "MSigDB HALLMARK_G2M_CHECKPOINT", doi = "10.1016/j.cels.2015.12.004", pmid = "26771021",
       version = "MSigDB 2026.1.Hs", overlap = "independent of all cohorts"),
  list(id = "CTRL_STROMA", label = "Overall stromal abundance (ESTIMATE stromal signature)", role = "control",
       genes = prog_mapped$STROMA, raw = estimate_stromal, direction = 1,
       source = "Yoshihara et al. 2013 Nat Commun, Supplementary Data 1 (stromal signature)", doi = "10.1038/ncomms3612", pmid = "24113773",
       version = "Nature Communications ESM MOESM488", overlap = "derived from bulk tumour compendia; independent of this project's cohorts"),
  list(id = "CTRL_IEG_STRESS", label = "Immediate-early / stress response (ischaemia and dissociation confounder)", role = "control",
       genes = prog_mapped$STRESS, raw = tca_stress, direction = 1,
       source = paste0("MSigDB C4:3CA ", stress_nm, " (Gavish et al. 2023 malignant meta-program)"), doi = "10.1038/s41586-023-06130-4", pmid = "37258682",
       version = "MSigDB 2026.1.Hs (3CA meta_programs 2023-07-13)", overlap = "3CA contributing studies contain no gastric cancer cohort"),
  list(id = "NEG_YAP_STROMAL", label = "Stroma-restricted YAP/TAZ-TEAD target genes (admixture negative control)", role = "negative_control",
       genes = yap_stromal, raw = yap_stromal, direction = 1,
       source = "consensus YAP/TAZ-TEAD targets appearing in >=2 independent source signatures, restricted to genes that are non-epithelium-dominant in HPA v25",
       doi = "multiple (see program_candidates_yap_ecm.csv)", pmid = "NA", version = "derived 2026-09-25",
       overlap = "signature sources are cell-line based; HPA is an independent normal-tissue reference"),
  list(id = "SENS_HRC_EPI_RESTRICTED", label = "coreHRC restricted to epithelium-dominant genes (pre-specified sensitivity)", role = "sensitivity",
       genes = hrc_epi_restricted, raw = hrc_epi_restricted, direction = 1,
       source = "coreHRC genes with HPA v25 specificity == epithelium_dominant", doi = "10.1038/s41586-022-05402-9", pmid = "36352230",
       version = "derived 2026-09-25", overlap = "as coreHRC"),
  list(id = "SENS_HRC_EPIHR", label = "Epithelial high-relapse program (pre-specified alternative definition)", role = "sensitivity",
       genes = prog_mapped$HRC_EPI, raw = hrc_epi, direction = 1,
       source = "Canellas-Socias et al. 2022 Nature, Supplementary Table 2 (epiHR == 1)", doi = "10.1038/s41586-022-05402-9", pmid = "36352230",
       version = "Nature ESM MOESM4", overlap = "as coreHRC")
)
for (k in names(compartment_sets)) defs[[length(defs) + 1]] <- list(
  id = k, label = paste0("Host compartment marker set: ", tolower(sub("^COMPART_", "", k))), role = "composition_covariate",
  genes = compartment_sets[[k]], raw = compartment_sets[[k]], direction = 1,
  source = "HPA v25 single-cell-type nCPM: genes with compartment max >= 10 nCPM and >= 3x every other compartment, top 40, excluding all program genes",
  doi = "10.1126/science.aal3321", pmid = "28495876", version = "HPA rna_single_cell_type downloaded 2026-09-25",
  overlap = "normal-tissue reference; independent of all analysis cohorts")

# ---- coverage per platform ----
ann <- fread(ann_f)
bulk314 <- readRDS(pm_in("02_data_processed", "bulk", "GSE314812_counts.rds"))
bulk237 <- readRDS(pm_in("02_data_processed", "bulk", "GSE237876_gene_counts.rds"))
sym_of <- function(ids) unique(ann[gene_id %in% ids & !is.na(symbol_for_scoring), symbol_for_scoring])
platform_genes <- list(
  GSE314812 = sym_of(rownames(bulk314$counts)),
  GSE237876 = sym_of(rownames(bulk237$counts))
)
for (g in c("GSE183904", "GSE163558", "GSE308231", "GSE239676")) {
  f <- pm_in("02_data_processed", "singlecell", paste0(g, "_genes_measured.csv.gz"))
  if (file.exists(f)) platform_genes[[g]] <- unique(fread(f)$symbol)
}
COV_MIN_FRAC <- 0.80; COV_MIN_N <- 10
cov <- rbindlist(lapply(defs, function(d) rbindlist(lapply(names(platform_genes), function(p) {
  n_def <- length(d$genes); n_meas <- length(intersect(d$genes, platform_genes[[p]]))
  data.table(program_id = d$id, platform = p, n_defined = n_def, n_measured = n_meas,
             frac_measured = round(n_meas / n_def, 4),
             evaluable = (n_meas / n_def >= COV_MIN_FRAC) & (n_meas >= min(COV_MIN_N, n_def)))
}))))
o_cov <- pm_out("02_data_processed", "frozen_programs", "program_coverage_by_platform.csv"); write_csv_atomic(cov, o_cov)

# ---- write definitions ----
defs_long <- rbindlist(lapply(defs, function(d) {
  m <- harm(d$raw)
  m <- merge(m, hpa_w[, .(gene_mapped = gene, hpa_specificity = specificity)], by = "gene_mapped", all.x = TRUE, sort = FALSE)
  data.table(program_id = d$id, biological_label = d$label, role = d$role,
             source_doi_or_url = d$doi, source_pmid = d$pmid, source_version = d$version,
             source_dataset_overlap_status = d$overlap,
             gene_original = m$gene_original, gene_mapped = m$gene_mapped, mapping_route = m$mapping,
             hpa_specificity = fifelse(is.na(m$hpa_specificity), "not_in_HPA", m$hpa_specificity),
             direction = d$direction, fixed_weight = 1,
             inclusion_reason = "member of the source definition as published (or of the pre-specified derived subset)",
             exclusion_reason = NA_character_, scoring_method = "within-sample equal-weight mean gene rank on the frozen common universe",
             score_scale = "0-1 relative transcriptional rank position")
}), fill = TRUE)
defs_long[, frozen_at := "2026-09-25"]
o_defs <- pm_out("02_data_processed", "frozen_programs", "program_definitions.csv"); write_csv_atomic(defs_long, o_defs)

per_prog <- lapply(defs, function(d) {
  g <- sort(unique(harm(d$raw)$gene_mapped))
  list(program_id = d$id, biological_label = d$label, role = d$role, source = d$source, doi = d$doi, pmid = d$pmid,
       source_version = d$version, dataset_overlap_status = d$overlap, direction = d$direction, weights = "equal",
       n_genes = length(g), genes = g, scoring_method = "within-sample equal-weight mean gene rank",
       coverage_rule = sprintf(">=%.0f%% of defined genes measured and >=%d genes", 100 * COV_MIN_FRAC, COV_MIN_N),
       definition_sha256 = sha256_obj(list(id = d$id, genes = g, direction = d$direction)))
})
names(per_prog) <- vapply(defs, function(d) d$id, "")
o_json <- pm_out("02_data_processed", "frozen_programs", "program_definitions.json"); write_json_atomic(per_prog, o_json)
o_rds <- pm_out("02_data_processed", "frozen_programs", "program_definitions.rds"); save_rds_atomic(per_prog, o_rds)

summ <- rbindlist(lapply(per_prog, function(p) data.table(program_id = p$program_id, role = p$role, n_genes = p$n_genes,
                                                          definition_sha256 = substr(p$definition_sha256, 1, 12))))
ov <- CJ(a = summ$program_id, b = summ$program_id)[a < b]
ov[, shared := mapply(function(x, y) length(intersect(per_prog[[x]]$genes, per_prog[[y]]$genes)), a, b)]
o_ov <- pm_out("02_data_processed", "frozen_programs", "program_pairwise_gene_overlap.csv"); write_csv_atomic(ov[shared > 0][order(-shared)], o_ov)
# ---- protocol lock -------------------------------------------------------------------
man <- fread(pm_record("00_admin", "data_manifest.csv"), colClasses = "character")
lock <- list(
  protocol_version = "v1",
  frozen_at = "2026-09-25",
  is_preregistered = FALSE,
  registration_note = "not registered on any platform; the term preregistered is not used in the manuscript",
  branch = "B (single frozen malignant-epithelial program; within-patient change)",
  primary_question = paste("Does an independently defined, malignant-epithelium-biased dissemination/relapse program",
                           "change reproducibly between paired primary gastric cancer and solid peritoneal metastasis",
                           "within the same patient, and does any change survive adjustment for host-tissue composition?"),
  primary_program = "HRC_CORE",
  primary_estimand = "mean of within-patient delta = score(solid peritoneal metastasis) - score(primary tumour); one pair per patient; patient-weighted",
  primary_test = "two-sided paired t-test on patient deltas; 95% t confidence interval; patient-level bootstrap (5000) as robustness",
  primary_cohort = "GSE314812 eligible paired patients (GC05 pre-treatment visit only)",
  roles = list(primary_test = "GSE314812", external_bulk_directional = "GSE237876",
               calibration_development = c("GSE183904 gastric primary/normal donors", "GSE163558"),
               cell_level_descriptive = c("GSE308231", "GSE183904 peritoneal donors"),
               fluid_context_separate_estimand = "GSE239676", spatial_single_case = "GSE251950",
               clinical_exploratory_default_off = "GSE62254/ACRG"),
  multiplicity = list(primary_family = list(members = "HRC_CORE", correction = "none (single test)"),
                      secondary_A = list(members = c("CTRL_EMT", "CTRL_PROLIFERATION", "CTRL_STROMA", "CTRL_IEG_STRESS", "NEG_YAP_STROMAL"), correction = "BH FDR"),
                      secondary_B = list(members = grep("^COMPART_", names(per_prog), value = TRUE), correction = "BH FDR"),
                      exploratory = c("gene-level DE", "pathway enrichment", "fluid context", "spatial", "clinical")),
  coverage_rule = list(min_fraction_measured = COV_MIN_FRAC, min_genes = COV_MIN_N,
                       failure_action = "program marked not evaluable on that platform; no score computed"),
  scoring = list(method = "within-sample equal-weight mean gene rank", scale = "0-1 relative rank",
                 universe_rule = "CPM >= 1 in >= 25% of samples on that platform; defined without any site/group information"),
  adjusted_models_prespecified = c("preservation change", "mesothelium + adipocyte compartment deltas", "stromal abundance delta", "proliferation delta"),
  sensitivity_prespecified = c("patient bootstrap", "median and sign test", "Wilcoxon signed rank", "leave-one-patient-out",
                               "same-preservation pairs only", "GC05 post-treatment pair", "epithelium-restricted subset", "epiHR alternative definition"),
  success_definitions = list(
    compatible = "independent effect same direction with interval compatible with the primary estimate",
    limited_directional = "point estimate same direction but interval wide (the maximum achievable in GSE237876, n=5)",
    inconsistent = "direction or magnitude in substantive conflict",
    uninformative = "too few patients or too imprecise",
    not_evaluable = "coverage, identity or availability failure"),
  evidence_ceiling = list(
    prohibited_claims = c("drives", "mediates", "causes", "therapeutic target", "predicts recurrence", "treatment benefit"),
    conditional_claims = list(independent_validation = "only if patients and development are truly independent and the endpoint identical; GSE237876 (n=5) is an independent cohort analysis with limited directional support, not validation",
                              cell_intrinsic = "requires within-cell-type measured evidence",
                              peritoneal_specific = "requires a directly estimable other-site contrast"),
    default_framing = "tissue-level transcriptional change; not malignant-cell-intrinsic activity"),
  contacted_published_results = "see topic_decision.md section 6; no target-cohort site-stratified expression was computed before this lock",
  program_definition_hashes = lapply(per_prog, function(p) p$definition_sha256),
  input_hashes = setNames(as.list(man$sha256), man$local_path),
  script_hashes = setNames(as.list(vapply(sort(list.files(pm_src("03_scripts"), pattern = "\\.R$", recursive = TRUE, full.names = TRUE)), sha256_file, "")),
                           sort(sub(paste0("^", pm_project_root(), "/"), "", list.files(pm_src("03_scripts"), pattern = "\\.R$", recursive = TRUE, full.names = TRUE)))),
  curation_hashes = setNames(as.list(vapply(sort(list.files(pm_src("03_scripts", "curation"), full.names = TRUE)), sha256_file, "")),
                             sort(basename(list.files(pm_src("03_scripts", "curation"))))),
  protocol_text_sha256 = sha256_file(pm_src("00_admin", "protocol_v1.md")),
  crosswalk_sha256 = sha256_file(pm_record("00_admin", "sample_crosswalk.csv"))
)
o_lock <- pm_out("00_admin", "protocol_lock.json"); write_json_atomic(lock, o_lock)

print(summ)
print(cov[program_id %in% c("HRC_CORE", "SENS_HRC_EPI_RESTRICTED", "CTRL_EMT", "CTRL_STROMA")])
step_end(outputs = c(o_hpa, o_cov, o_defs, o_json, o_rds, o_ov, o_lock))
