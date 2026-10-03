# Stage 2: sample and independence audit. Builds one row per biological specimen /
# sequencing object for every candidate source from GEO SOFT records, paper supplementary
# tables and the small curated maps in 03_scripts/curation (each with a source locator).
# Writes 00_admin/sample_crosswalk.csv, overlap_register.csv, eligibility_flow.csv and
# 04_results/audit/{source_counts,invariants}.csv. No expression values are summarised
# by site; the only expression-file reads are headers/column names for matching.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
suppressPackageStartupMessages(library(readxl))

meta <- function(gse) fread(pm_in("02_data_processed", "metadata", paste0(gse, "_geo_samples.csv")), colClasses = "character")
cur <- function(f) fread(pm_src("03_scripts", "curation", f), colClasses = "character")
first_field <- function(x) trimws(sub(" \\|\\|.*$", "", x))

inputs <- c(
  pm_in("02_data_processed", "metadata", paste0(c("GSE314812", "GSE237876", "GSE183904", "GSE163558", "GSE308231", "GSE239676", "GSE234129", "GSE228598"), "_geo_samples.csv")),
  pm_raw("GSE314812", "paper_supplement", "42003_2026_9906_MOESM3_ESM_SupplementaryData1.xlsx"),
  pm_raw("GSE314812", "geo", "GSE314812_rsem.merged.gene_counts.tsv.gz"),
  pm_raw("GSE237876", "paper_kim2024", "CAC2-44-514-s002.xlsx"),
  pm_raw("GSE237876", "GSE237876_RAW.tar"),
  pm_raw("other_sources", "GSE183904", "GSE183904_family.soft.gz"),
  pm_raw("other_sources", "GSE183904", "paper", "suppl", "cd-21-0683_supp.xlsx"),
  pm_raw("other_sources", "GSE239676", "GSE239676_barcodes.tsv.gz"),
  pm_raw("other_sources", "GSE239676", "GSE239676_meta.tsv.gz"),
  pm_raw("other_sources", "GSE234129", "GSE234129_meta.tsv.gz"),
  pm_src("03_scripts", "curation", c("GSE163558_patient_map.csv", "GSE314812_visit_rules.csv"))
)
step_begin("02_audit_samples", inputs = inputs)

COLS <- c("dataset_id", "source_study_id", "publication_pmid", "gsm_or_run_id", "expression_column",
          "original_sample_id", "original_patient_id", "canonical_patient_id", "specimen_id", "region_id",
          "technical_replicate_id", "visit", "anatomical_site", "specimen_context", "lesion_role", "pm_status",
          "primary_pm_pair_id", "paired_identity_evidence", "collection_time_or_order", "therapy_before_sampling",
          "therapy_regimen", "preservation", "assay", "genome_build", "matrix_unit", "histology", "stage", "sex",
          "age_if_public", "purity_published", "metadata_source", "metadata_source_locator", "identity_confidence",
          "eligible_primary", "eligible_validation", "exclusion_reason", "site_raw")
fill <- function(dt) { for (c in setdiff(COLS, names(dt))) dt[, (c) := NA_character_]; dt[, ..COLS] }

# ---------------- GSE314812 (Charton 2026; primary within-patient target) -------------
g <- meta("GSE314812")
g[, library := sub("^Library name: ", "", first_field(description))]
g[, geo_label := sub("^.*patient ", "", title)]
vr <- cur("GSE314812_visit_rules.csv")
g <- merge(g, vr, by.x = "geo_label", by.y = "geo_patient_label", all.x = TRUE)
g[is.na(canonical_patient_id), `:=`(canonical_patient_id = geo_label, visit = "1",
                                    therapy_before_sampling = "unknown", collection_time_or_order = "unknown (single surgical visit reported)",
                                    identity_confidence = "high")]
sd1 <- as.data.table(read_excel(pm_raw("GSE314812", "paper_supplement", "42003_2026_9906_MOESM3_ESM_SupplementaryData1.xlsx"), skip = 2))
setnames(sd1, c("no", "patient", "gender", "age_range", "lauren", "stage", "msi", "ebv", "her2ihc", "who", "tcga"))
sd1 <- sd1[grepl("^GC[0-9]+$", patient)]
g <- merge(g, sd1[, .(patient, lauren, stage, age_range, gender)], by.x = "canonical_patient_id", by.y = "patient", all.x = TRUE)
site_map <- c("GC Primary tumor" = "stomach", "GC Peritoneal metastatic tumor" = "peritoneum", "NAT" = "stomach",
              "GC Omentum metastatic tumor" = "omentum", "GC Lymph node metastatic tumor" = "lymph_node")
role_map <- c("GC Primary tumor" = "primary_tumor", "GC Peritoneal metastatic tumor" = "solid_peritoneal_metastasis",
              "NAT" = "normal_adjacent_gastric", "GC Omentum metastatic tumor" = "omental_metastasis",
              "GC Lymph node metastatic tumor" = "nodal_metastasis")
hdr <- names(fread(cmd = sprintf("gzip -dc %s | head -2", shQuote(pm_raw("GSE314812", "geo", "GSE314812_rsem.merged.gene_counts.tsv.gz")))))
x314 <- g[, .(dataset_id = "GSE314812", source_study_id = "Charton2026_CommunBiol", publication_pmid = "41882239",
              gsm_or_run_id = gsm, expression_column = ifelse(library %in% hdr, library, NA_character_),
              original_sample_id = library, original_patient_id = geo_label, canonical_patient_id,
              specimen_id = paste0("GSE314812_", library), region_id = "single_region_reported", technical_replicate_id = "1",
              visit, anatomical_site = site_map[char_tissue], specimen_context = "solid_tissue",
              lesion_role = role_map[char_tissue], pm_status = "PM_patient_cohort",
              collection_time_or_order, therapy_before_sampling, therapy_regimen = "unknown",
              preservation = char_preservation_method,
              assay = "bulk whole-transcriptome RNA-seq (TruSeq RNA Access or exome capture; per-sample kit not reported)",
              genome_build = "GRCh37", matrix_unit = "RSEM expected counts (gene-level, nf-core star_rsem)",
              histology = paste0("Lauren=", lauren), stage = stage, sex = char_sex, age_if_public = char_age,
              metadata_source = "GEO SOFT; Charton 2026 Supplementary Data 1; curation/GSE314812_visit_rules.csv",
              metadata_source_locator = paste0("GSE314812_family.soft.gz:", gsm, "; SD1 row ", canonical_patient_id),
              identity_confidence, site_raw = char_tissue)]
x314[, has_pt := any(lesion_role == "primary_tumor"), by = .(canonical_patient_id, visit)]
x314[, has_pm := any(lesion_role == "solid_peritoneal_metastasis"), by = .(canonical_patient_id, visit)]
x314[has_pt & has_pm, primary_pm_pair_id := paste0("GSE314812_", canonical_patient_id, "_v", visit)]
x314[!is.na(primary_pm_pair_id), paired_identity_evidence := "same GEO patient label and SD1 patient row (clinical table key)"]
# primary eligibility: earliest eligible visit (pre-treatment preferred), primary or solid PM only
x314[, first_pair_visit := if (any(!is.na(primary_pm_pair_id))) min(as.integer(visit[!is.na(primary_pm_pair_id)])) else NA_integer_, by = canonical_patient_id]
x314[, eligible_primary := !is.na(primary_pm_pair_id) & as.integer(visit) == first_pair_visit &
       lesion_role %in% c("primary_tumor", "solid_peritoneal_metastasis")]
x314[, exclusion_reason := fifelse(eligible_primary, "",
        fifelse(lesion_role %in% c("normal_adjacent_gastric"), "normal tissue (descriptive only)",
        fifelse(lesion_role %in% c("omental_metastasis", "nodal_metastasis"), "other metastatic site (descriptive / pre-specified extension)",
        fifelse(!is.na(primary_pm_pair_id), "later visit of a patient already contributing an earlier pair (post-treatment recurrence; sensitivity only)",
        fifelse(lesion_role == "solid_peritoneal_metastasis", "PM without same-patient primary", "primary without same-patient solid PM")))))]
x314[, eligible_validation := FALSE]
x314[, c("has_pt", "has_pm", "first_pair_visit") := NULL]

# ---------------- GSE237876 (Kim 2024; independent bulk paired source) -----------------
k <- meta("GSE237876")
k[, file := basename(supplementary_file_1)]
k[, pat := sub("^.*_GCMET_([0-9]+)_.*$", "\\1", file)]
k[, tag := sub("^.*_GCMET_[0-9]+_(.*)_exp\\.txt\\.gz$", "\\1", file)]
k[, sample_id := paste0("GCM", substr(pat, 2, 3), "_", tag)]
s2 <- as.data.table(read_excel(pm_raw("GSE237876", "paper_kim2024", "CAC2-44-514-s002.xlsx"), sheet = "Sup. Table S2", skip = 1))
setnames(s2, 1:7, c("patient", "sample_id", "tumor_type", "timing", "chemo", "lauren", "purity"))
s2 <- s2[grepl("^GCM", patient)]
s1 <- as.data.table(read_excel(pm_raw("GSE237876", "paper_kim2024", "CAC2-44-514-s002.xlsx"), sheet = "Sup. Table S1", skip = 2))
setnames(s1, 1:8, c("patient", "age", "sex", "tnm", "lauren_s1", "histology", "msi", "ebv"))
s1 <- s1[grepl("^GCM", patient)]
k <- merge(k, s2, by = "sample_id", all.x = TRUE)
k <- merge(k, s1[, .(patient, age, sex, tnm, histology)], by = "patient", all.x = TRUE)
tar_list <- system2("tar", c("-tf", shQuote(pm_raw("GSE237876", "GSE237876_RAW.tar"))), stdout = TRUE)
site237 <- function(tag) fifelse(grepl("^T", tag), "stomach", fifelse(grepl("^perit", tag), "peritoneum",
                         fifelse(grepl("^liver", tag), "liver", fifelse(grepl("^ovary", tag), "ovary", fifelse(grepl("^lung", tag), "lung", "unknown")))))
x237 <- k[, .(dataset_id = "GSE237876", source_study_id = "Kim2024_CancerCommun", publication_pmid = "38498378",
              gsm_or_run_id = gsm, expression_column = ifelse(file %in% tar_list, file, NA_character_),
              original_sample_id = sample_id, original_patient_id = patient, canonical_patient_id = patient,
              specimen_id = paste0("GSE237876_", sample_id),
              region_id = fifelse(grepl("^T", tag), tag, tag), technical_replicate_id = "1", visit = NA_character_,
              anatomical_site = site237(tag), specimen_context = "solid_tissue",
              lesion_role = fifelse(grepl("^T", tag), "primary_tumor", fifelse(grepl("^perit", tag), "solid_peritoneal_metastasis",
                            paste0(site237(tag), "_metastasis"))),
              pm_status = fifelse(patient %in% unique(patient[grepl("^perit", tag)]), "PM_patient", "no_PM_specimen"),
              collection_time_or_order = paste0("metastasis timing (Table S2): ", timing),
              therapy_before_sampling = fifelse(chemo == "naïve", "chemotherapy-naive (Table S2)",
                                         fifelse(chemo == "exposed", "chemotherapy-exposed (Table S2)", "unknown")),
              therapy_regimen = "see Kim 2024 Supplementary Fig. S1 (regimens by patient)",
              preservation = "FFPE", assay = "bulk total RNA-seq (TruSeq Stranded Total RNA), PE101",
              genome_build = "hg19 (bcbio transcriptome)", matrix_unit = "Salmon quant.sf per sample (transcript-level NumReads + TPM)",
              histology = paste0("Lauren=", lauren, "; ", histology), stage = tnm, sex = sex, age_if_public = as.character(age),
              purity_published = as.character(purity),
              metadata_source = "GEO SOFT; Kim 2024 Supplementary Tables S1-S2",
              metadata_source_locator = paste0("GSE237876_family.soft.gz:", gsm, "; CAC2-44-514-s002.xlsx Table S2 ", sample_id),
              identity_confidence = "high", site_raw = title, tag, patient, timing, chemo)]
pm_pat237 <- x237[lesion_role == "solid_peritoneal_metastasis", unique(canonical_patient_id)]
x237[canonical_patient_id %in% pm_pat237 & lesion_role %in% c("primary_tumor", "solid_peritoneal_metastasis"),
     primary_pm_pair_id := paste0("GSE237876_", canonical_patient_id)]
x237[!is.na(primary_pm_pair_id), paired_identity_evidence := "GEO patient number = Table S2 PatientID; sister WES phylogeny (Lee 2023 BJC)"]
x237[canonical_patient_id == "GCM05", identity_confidence := "moderate: perit label in GEO/Table S2/WES tree but clinical-course figures label the resected metastasis as liver"]
x237[, eligible_primary := FALSE]
x237[, eligible_validation := !is.na(primary_pm_pair_id)]
x237[, exclusion_reason := fifelse(eligible_validation, "", fifelse(lesion_role == "primary_tumor", "patient without solid PM specimen", "non-peritoneal metastatic site (descriptive / other-site contrast)"))]
x237[, c("tag", "patient", "timing", "chemo") := NULL]

# ---------------- GSE183904 (Kumar 2022; single-cell calibration) --------------------
soft183 <- readLines(gzfile(pm_raw("other_sources", "GSE183904", "GSE183904_family.soft.gz")), warn = FALSE)
tab <- sub("^!Series_summary = ", "", grep("^!Series_summary = NGCII", soft183, value = TRUE))
tab <- rbindlist(lapply(strsplit(tab, "\t"), function(v) data.table(patient = trimws(v[1]), Primary_Tumor = v[2], Primary_Normal = v[3],
                                                                   Peritoneal_Tumor = v[4], Peritoneal_Normal = v[5])))
long <- melt(tab, id.vars = "patient", variable.name = "role", value.name = "file")[file != "-" & !is.na(file)]
long[, sample_n := sub("\\.csv$", "", trimws(file))]
m183 <- meta("GSE183904")
m183[, sample_n := sub(":.*$", "", title)]
m183 <- merge(m183, long, by = "sample_n", all.x = TRUE)
ks1 <- as.data.table(read_excel(pm_raw("other_sources", "GSE183904", "paper", "suppl", "cd-21-0683_supp.xlsx"), sheet = "Supp Table 1", skip = 1))
setnames(ks1, 1:7, c("no", "patient", "age", "sex", "location", "stage", "lauren"))
ks1[, patient := trimws(patient)]
m183 <- merge(m183, ks1[, .(patient, age, sex, stage, lauren)], by = "patient", all.x = TRUE)
r183 <- c(Primary_Tumor = "primary_tumor", Primary_Normal = "normal_adjacent_gastric", Peritoneal_Tumor = "solid_peritoneal_metastasis",
          Peritoneal_Normal = "normal_peritoneum_not_PM")
x183 <- m183[, .(dataset_id = "GSE183904", source_study_id = "Kumar2022_CancerDiscov", publication_pmid = "34642171",
                 gsm_or_run_id = gsm, expression_column = basename(supplementary_file_1), original_sample_id = sample_n,
                 original_patient_id = patient, canonical_patient_id = patient, specimen_id = paste0("GSE183904_", sample_n),
                 region_id = "1", technical_replicate_id = "1", visit = NA_character_,
                 anatomical_site = fifelse(grepl("^Peritoneal", role), "peritoneum", "stomach"), specimen_context = "solid_tissue",
                 lesion_role = r183[as.character(role)], pm_status = fifelse(patient %in% long[role == "Peritoneal_Tumor", patient], "PM_patient", "unknown"),
                 collection_time_or_order = fifelse(grepl("^Peritoneal", role), "diagnostic laparoscopy (Methods; cohort-level)", "resection or endoscopic biopsy (Methods; cohort-level)"),
                 therapy_before_sampling = "unknown", therapy_regimen = "unknown", preservation = "fresh (dissociated for scRNA-seq)",
                 assay = "scRNA-seq 10x 5'", genome_build = "GRCh38 (Cell Ranger 3.0)", matrix_unit = "raw UMI counts (gene x cell CSV)",
                 histology = paste0("Lauren=", lauren), stage = as.character(stage), sex = sex, age_if_public = as.character(age),
                 metadata_source = "GEO SOFT Series_summary mapping table; Kumar 2022 Supp Table 1",
                 metadata_source_locator = paste0("GSE183904_family.soft.gz Series_summary row ", patient, "; ", gsm),
                 identity_confidence = "high", site_raw = title)]
x183[canonical_patient_id %in% long[role == "Peritoneal_Tumor", patient] & canonical_patient_id %in% long[role == "Primary_Tumor", patient] &
       lesion_role %in% c("primary_tumor", "solid_peritoneal_metastasis"), primary_pm_pair_id := paste0("GSE183904_", canonical_patient_id)]
x183[, `:=`(eligible_primary = FALSE, eligible_validation = FALSE)]
x183[, exclusion_reason := fifelse(lesion_role %in% c("primary_tumor", "normal_adjacent_gastric") & pm_status != "PM_patient",
                                   "calibration donor (cell identity/measurement; development side)",
                                   "peritoneal-sampled donor: held out from calibration; descriptive cell-level only")]

# ---------------- GSE163558 (Jiang 2022; calibration across chemistry) --------------
m163 <- meta("GSE163558")
m163[, sample_code := sub("^.* ", "", trimws(title))]
pm163 <- merge(m163, cur("GSE163558_patient_map.csv"), by = "sample_code", all.x = TRUE)
lr163 <- function(code) fifelse(grepl("^PT", code), "primary_tumor", fifelse(grepl("^NT", code), "normal_adjacent_gastric",
                          fifelse(grepl("^P[0-9]", code), "solid_peritoneal_metastasis", fifelse(grepl("^Li", code), "liver_metastasis",
                          fifelse(grepl("^LN", code), "nodal_metastasis", fifelse(grepl("^O", code), "ovary_metastasis", "unknown"))))))
x163 <- pm163[, .(dataset_id = "GSE163558", source_study_id = "Jiang2022_ClinTranslMed", publication_pmid = "35184420",
                  gsm_or_run_id = gsm, expression_column = gsm, original_sample_id = sample_code, original_patient_id = patient,
                  canonical_patient_id = paste0("GSE163558_", patient), specimen_id = paste0("GSE163558_", sample_code),
                  region_id = "1", technical_replicate_id = "1", visit = NA_character_,
                  anatomical_site = c(primary_tumor = "stomach", normal_adjacent_gastric = "stomach", solid_peritoneal_metastasis = "peritoneum",
                                      liver_metastasis = "liver", nodal_metastasis = "lymph_node", ovary_metastasis = "ovary")[lr163(sample_code)],
                  specimen_context = "solid_tissue", lesion_role = lr163(sample_code), pm_status = fifelse(patient == "Patient6", "PM_patient", "no_PM_specimen"),
                  collection_time_or_order = paste0("sampling: ", sampling_method), therapy_before_sampling = "untreated (GEO; paper)",
                  therapy_regimen = "none before sampling", preservation = "fresh (dissociated for scRNA-seq)", assay = "scRNA-seq 10x 3' v3",
                  genome_build = "GRCh38 (Cell Ranger 3.0.2)", matrix_unit = "Cell Ranger filtered UMI counts (MTX)",
                  histology = paste0("Lauren=", lauren), stage = clinical_note, sex = sex, age_if_public = age_decade,
                  metadata_source = "GEO SOFT; Jiang 2022 Supplementary Table S1 (curated map)",
                  metadata_source_locator = paste0(gsm, "; ", metadata_source_locator), identity_confidence = "high", site_raw = title)]
x163[, `:=`(eligible_primary = FALSE, eligible_validation = FALSE)]
x163[, exclusion_reason := fifelse(lesion_role == "solid_peritoneal_metastasis", "single PM without same-patient primary; descriptive cell-level only",
                                   "calibration donor (cell identity/measurement across chemistry)")]

# ---------------- GSE308231 (Wang 2025; solid PM scRNA, donors unknown) -------------
m308 <- meta("GSE308231")
x308 <- m308[, .(dataset_id = "GSE308231", source_study_id = "WangQ2025_npjDigitMed", publication_pmid = "41360923",
                 gsm_or_run_id = gsm, expression_column = sub("^([A-Za-z]+)_?([0-9]+)$", "\\1_\\2", title), original_sample_id = title,
                 original_patient_id = "not reported", canonical_patient_id = paste0("GSE308231_unknown_donor_", title),
                 specimen_id = paste0("GSE308231_", title), region_id = "1", technical_replicate_id = "1", visit = NA_character_,
                 anatomical_site = fifelse(grepl("^Ca", title), "stomach", "peritoneum"), specimen_context = "solid_tissue",
                 lesion_role = fifelse(grepl("^Ca", title), "primary_tumor", "solid_peritoneal_metastasis"), pm_status = "unknown",
                 paired_identity_evidence = "none (no donor key in GEO, BioSample, paper or code)",
                 collection_time_or_order = "fresh surgical specimens (series summary); timing unknown", therapy_before_sampling = "unknown",
                 therapy_regimen = "unknown", preservation = "fresh (dissociated for scRNA-seq)", assay = "scRNA-seq (chemistry conflicting: SeekOne vs cellranger-7.0.1 header)",
                 genome_build = "GRCh38 (36,601 features)", matrix_unit = "filtered UMI counts (MTX integer)",
                 metadata_source = "GEO SOFT; BioSample", metadata_source_locator = paste0("GSE308231_family.soft.gz:", gsm),
                 identity_confidence = "donor identity unknown; libraries may share donors", site_raw = characteristics_raw)]
x308[, `:=`(eligible_primary = FALSE, eligible_validation = FALSE,
            exclusion_reason = "donor identity/pairing unverifiable; library-level descriptive support only")]

# ---------------- GSE239676 + GSE234129 (Cheng 2024 / Wang 2023; ascites context) ---
m239 <- meta("GSE239676")
m239[, code := sub("\\s*from .*$", "", title)]   # note: one GEO title lacks the space ("PBMCfrom Pt-11")
m239[, pt := sub("^.*from (Pt-[0-9]+),.*$", "\\1", title)]
ctx239 <- c(N = "solid_tissue", PRI = "solid_tissue", PC = "ascites", LM = "solid_tissue", OV = "solid_tissue", PB = "blood", PBMC = "blood")
role239 <- c(N = "normal_adjacent_gastric", PRI = "primary_tumor", PC = "malignant_ascites", LM = "liver_metastasis", OV = "ovary_metastasis",
             PB = "peripheral_blood", PBMC = "pbmc")
site239 <- c(N = "stomach", PRI = "stomach", PC = "peritoneal_cavity_fluid", LM = "liver", OV = "ovary", PB = "blood", PBMC = "blood")
x239 <- m239[, .(dataset_id = "GSE239676", source_study_id = "Cheng2024_Gastroenterology", publication_pmid = "39097198",
                 gsm_or_run_id = gsm, expression_column = NA_character_, original_sample_id = title, original_patient_id = pt,
                 canonical_patient_id = paste0("GSE239676_", pt), specimen_id = paste0("GSE239676_", pt, "_", code),
                 region_id = "1", technical_replicate_id = "1", visit = NA_character_, anatomical_site = site239[code],
                 specimen_context = ctx239[code], lesion_role = role239[code], pm_status = "clinical PM status per patient not public",
                 collection_time_or_order = "stage IV at diagnosis; ascites at therapeutic paracentesis (GSE234129 design, subset)",
                 therapy_before_sampling = "treatment-naive (GEO)", therapy_regimen = "none before sampling", preservation = "fresh",
                 assay = "scRNA-seq 10x 5'", genome_build = "GRCh38 (Cell Ranger 3.0)", matrix_unit = "UMI counts, 8,630 genes retained (MTX integer)",
                 metadata_source = "GEO SOFT; GSE239676_meta.tsv", metadata_source_locator = paste0("GSE239676_family.soft.gz:", gsm),
                 identity_confidence = "high", site_raw = source_name_ch1, code)]
bar239 <- fread(pm_raw("other_sources", "GSE239676", "GSE239676_barcodes.tsv.gz"), header = FALSE, col.names = "bc")
meta239 <- fread(pm_raw("other_sources", "GSE239676", "GSE239676_meta.tsv.gz"))
stopifnot(nrow(bar239) == nrow(meta239))
meta239[, bc := bar239$bc]
smp_key <- unique(meta239[, .(Sample, Patient, Tissue)])
# map meta sample codes to GEO titles via patient number and tissue code
tis2code <- c(Ad = "N", P = "PRI", As = "PC", Li = "LM", Ov = "OV", PB = "PB", PBMC = "PBMC")
smp_key[, pt := paste0("Pt-", as.integer(sub("^Pt", "", Patient)))]
smp_key[, code := tis2code[Tissue]]
x239 <- merge(x239, smp_key[, .(pt, code, Sample)], by.x = c("original_patient_id", "code"), by.y = c("pt", "code"), all.x = TRUE)
x239[, expression_column := Sample]
x239[, c("code", "Sample") := NULL]
pt_pri <- x239[lesion_role == "primary_tumor", unique(canonical_patient_id)]
pt_as <- x239[lesion_role == "malignant_ascites", unique(canonical_patient_id)]
x239[canonical_patient_id %in% intersect(pt_pri, pt_as) & lesion_role %in% c("primary_tumor", "malignant_ascites"),
     primary_pm_pair_id := paste0("GSE239676_", original_patient_id, "_primary_ascites")]
x239[!is.na(primary_pm_pair_id), paired_identity_evidence := "same GEO patient label and meta.tsv Patient field"]
x239[, `:=`(eligible_primary = FALSE, eligible_validation = FALSE,
            exclusion_reason = fifelse(!is.na(primary_pm_pair_id), "fluid-context paired contrast (primary vs ascites); not solid PM",
                                       "context/descriptive only"))]

# GSE234129 is a strict subset (same libraries): barcode-level overlap check
m234 <- meta("GSE234129")
meta234 <- fread(pm_raw("other_sources", "GSE234129", "GSE234129_meta.tsv.gz"))
meta234[, core := substr(sub("^.*?([ACGT]{16}).*$", "\\1", cell_barcodes), 1, 16)]
meta239[, core := substr(sub("^.*_([ACGT]{16}).*$", "\\1", bc), 1, 16)]
ov <- rbindlist(lapply(split(meta234, meta234$sample), function(d) {
  sh <- meta239[core %in% d$core, .N, by = Sample][order(-N)]
  data.table(gse234129_sample = d$sample[1], n_cells = nrow(d), best_gse239676_sample = sh$Sample[1],
             shared_fraction = round(sh$N[1] / nrow(d), 3),
             second_best_fraction = if (nrow(sh) > 1) round(sh$N[2] / nrow(d), 3) else 0)
}))
x234 <- m234[, .(dataset_id = "GSE234129", source_study_id = "WangR2023_CancerCell", publication_pmid = "37419119",
                 gsm_or_run_id = gsm, expression_column = sub("^MDA_", "", title), original_sample_id = title,
                 original_patient_id = sub("^MDA_(Pt[0-9]+)-.*$", "\\1", title),
                 canonical_patient_id = paste0("GSE239676_Pt-", sub("^MDA_Pt([0-9]+)-.*$", "\\1", title)),
                 specimen_id = paste0("GSE234129_", title), specimen_context = fifelse(grepl("-As", title), "ascites", fifelse(grepl("-PB", title), "blood", "solid_tissue")),
                 lesion_role = "duplicate_library_of_GSE239676", therapy_before_sampling = "treatment-naive (GEO)",
                 metadata_source = "GEO SOFT; barcode overlap with GSE239676", metadata_source_locator = paste0("GSE234129_family.soft.gz:", gsm),
                 identity_confidence = "same library as GSE239676 (barcode overlap)", site_raw = source_name_ch1,
                 eligible_primary = FALSE, eligible_validation = FALSE,
                 exclusion_reason = "duplicate of GSE239676 libraries; counted once (GSE239676)")]

# ---------------- GSE228598 (Sullivan 2024; peritoneal fluid, no primary) -----------
m228 <- meta("GSE228598")
x228 <- m228[, .(dataset_id = "GSE228598", source_study_id = "Sullivan2024_IJMS", publication_pmid = "38612926", gsm_or_run_id = gsm,
                 expression_column = gsm, original_sample_id = title, original_patient_id = sub(", scRNAseq", "", title),
                 canonical_patient_id = paste0("GSE228598_", gsub(" ", "", sub(", scRNAseq", "", title))), specimen_id = paste0("GSE228598_", gsm),
                 anatomical_site = "peritoneal_cavity_fluid",
                 specimen_context = fifelse(grepl("Washings", char_sample_type), "lavage", "ascites"),
                 lesion_role = fifelse(grepl("Washings", char_sample_type), "peritoneal_washings", "malignant_ascites"),
                 pm_status = "not reported per patient", therapy_before_sampling = "unknown (cohort mixes naive and pretreated)",
                 assay = "scRNA-seq 10x 3'", genome_build = "GRCh38", matrix_unit = "Cell Ranger filtered UMI counts (MTX)",
                 metadata_source = "GEO SOFT", metadata_source_locator = paste0("GSE228598_family.soft.gz:", gsm),
                 identity_confidence = "high (one sample per patient)", site_raw = char_sample_type,
                 eligible_primary = FALSE, eligible_validation = FALSE,
                 exclusion_reason = "fluid context without primary; descriptive only")]

cw <- rbindlist(lapply(list(x314, x237, x183, x163, x308, x239, x234, x228), fill), fill = TRUE)
cw[is.na(eligible_primary), eligible_primary := FALSE]
cw[is.na(eligible_validation), eligible_validation := FALSE]
for (c in COLS) cw[is.na(get(c)) & !c %in% c("eligible_primary", "eligible_validation"), (c) := "unknown"]
o_cw <- pm_out("00_admin", "sample_crosswalk.csv"); write_csv_atomic(cw, o_cw)

# ---------------- overlap register -------------------------------------------------
ovr <- rbindlist(list(
  data.table(source_a = "GSE234129", source_b = "GSE239676", shared = "all 17 libraries (6 patients)",
             evidence = paste0("barcode overlap per sample: min shared ", min(ov$shared_fraction), ", max second-best ", max(ov$second_best_fraction)),
             handling = "count once as GSE239676", independent = FALSE),
  data.table(source_a = "GSE237876", source_b = "PRJNA992126 / Lee 2023 BJC WES", shared = "same 14 patients (duplicate FASTQ submission; sister WES)",
             evidence = "identical SRA spot/base counts; same GCM IDs", handling = "no new patients", independent = FALSE),
  data.table(source_a = "GSE237876", source_b = "msEMT signature (Kim 2024)", shared = "all 66 specimens used to derive msEMT",
             evidence = "Kim 2024 Supplementary methods", handling = "GSE237876 cannot validate msEMT or its subsets", independent = FALSE),
  data.table(source_a = "GSE183904 + GSE163558", source_b = "Charton 2026 deconvolution/pseudobulk reference", shared = "same public scRNA donors",
             evidence = "Charton 2026 Methods, SD9-SD10", handling = "reference reuse must be declared; not independent of Charton cell-type analysis", independent = FALSE),
  data.table(source_a = "GSE183904, GSE163558, GSE308231", source_b = "Wang Q 2025 / Sun M 2026 / Hao 2026 / Go 2026 integrations", shared = "same donors",
             evidence = "papers' data statements", handling = "integrations add no patients", independent = FALSE),
  data.table(source_a = "GSE314812", source_b = "GSE237876 / GSE183904 / GSE163558 / GSE308231 / GSE239676", shared = "none known",
             evidence = "different institutions (SNUBH Korea; Severance Korea; NUH Singapore; ZJU China; Ruijin China; Zhejiang Cancer Hospital)",
             handling = "treated as patient-independent (not genotype-verified)", independent = TRUE),
  data.table(source_a = "GSE308231", source_b = "GSE251950", shared = "none (Ruijin vs Catholic Univ. Korea)", evidence = "GEO contact institutions",
             handling = "independent", independent = TRUE)
))
o_ov <- pm_out("00_admin", "overlap_register.csv"); write_csv_atomic(ovr, o_ov)
o_ov2 <- pm_out("04_results", "audit", "gse234129_gse239676_barcode_overlap.csv"); write_csv_atomic(ov, o_ov2)

# ---------------- counts, eligibility flow and invariants --------------------------
cnt <- cw[dataset_id != "GSE234129", .(
  n_specimens = .N,
  n_patients_identifiable = uniqueN(canonical_patient_id[!grepl("unknown_donor", canonical_patient_id)]),
  n_primary = sum(lesion_role == "primary_tumor"), n_solid_pm = sum(lesion_role == "solid_peritoneal_metastasis"),
  n_ascites = sum(lesion_role == "malignant_ascites"), n_lavage = sum(lesion_role == "peritoneal_washings"),
  n_normal = sum(grepl("^normal", lesion_role)), n_other_met = sum(grepl("metastasis$", lesion_role) & lesion_role != "solid_peritoneal_metastasis"),
  n_patients_primary_and_solid_pm = uniqueN(canonical_patient_id[!is.na(primary_pm_pair_id) & primary_pm_pair_id != "unknown" & grepl("solid|primary", lesion_role) & !grepl("ascites", primary_pm_pair_id)]),
  n_pair_units_primary_solid_pm = uniqueN(primary_pm_pair_id[primary_pm_pair_id != "unknown" & !grepl("ascites", primary_pm_pair_id)]),
  n_patients_primary_and_ascites = uniqueN(canonical_patient_id[grepl("ascites", primary_pm_pair_id)]),
  n_therapy_known = uniqueN(canonical_patient_id[!therapy_before_sampling %in% c("unknown")]),
  n_ffpe = sum(preservation == "FFPE"), n_fresh_frozen = sum(preservation %in% c("fresh-frozen")),
  n_eligible_primary_patients = uniqueN(canonical_patient_id[eligible_primary]),
  n_eligible_validation_patients = uniqueN(canonical_patient_id[eligible_validation])
), by = dataset_id]
o_cnt <- pm_out("04_results", "audit", "source_counts.csv"); write_csv_atomic(cnt, o_cnt)

flow <- rbindlist(list(
  data.table(dataset_id = "GSE314812", step = c("GEO specimens", "independent patients", "patients with any solid PM",
                                                "patient-visit pair units (primary + solid PM)", "patients with >=1 pair unit",
                                                "primary analysis patients (earliest pair per patient)"),
             n = c(nrow(x314), uniqueN(x314$canonical_patient_id), uniqueN(x314[lesion_role == "solid_peritoneal_metastasis", canonical_patient_id]),
                   uniqueN(na.omit(x314$primary_pm_pair_id)), uniqueN(x314[!is.na(primary_pm_pair_id), canonical_patient_id]),
                   uniqueN(x314[eligible_primary == TRUE, canonical_patient_id])),
             unit = c("specimen", "patient", "patient", "pair unit", "patient", "patient")),
  data.table(dataset_id = "GSE237876", step = c("GEO specimens", "patients", "solid PM specimens", "patients with primary + solid PM",
                                                "of which synchronous and chemo-naive PM", "of which identity-conflict-free (excluding GCM05)"),
             n = c(nrow(x237), uniqueN(x237$canonical_patient_id), sum(x237$lesion_role == "solid_peritoneal_metastasis"),
                   uniqueN(x237[eligible_validation == TRUE, canonical_patient_id]),
                   uniqueN(k[grepl("^perit", tag) & timing == "synchronous" & chemo == "naïve", patient]),
                   uniqueN(x237[eligible_validation == TRUE & canonical_patient_id != "GCM05", canonical_patient_id])),
             unit = c("specimen", "patient", "specimen", "patient", "patient", "patient"))
))
o_flow <- pm_out("00_admin", "eligibility_flow.csv"); write_csv_atomic(flow, o_flow)

inv <- rbindlist(list(
  data.table(check = "GSE314812 expression columns matched to GSM", value = sprintf("%d/%d", sum(x314$expression_column != "unknown" & !is.na(x314$expression_column)), nrow(x314)),
             pass = all(!is.na(x314$expression_column)) && length(setdiff(hdr[-(1:2)], x314$original_sample_id)) == 0),
  data.table(check = "GSE314812 duplicate specimen keys", value = as.character(sum(duplicated(x314$specimen_id))), pass = !any(duplicated(x314$specimen_id))),
  data.table(check = "GSE314812 each primary-analysis patient has exactly one primary and one solid PM",
             value = paste(x314[eligible_primary == TRUE, .N, by = .(canonical_patient_id, lesion_role)][, unique(N)], collapse = ","),
             pass = x314[eligible_primary == TRUE, .N, by = .(canonical_patient_id, lesion_role)][, all(N == 1)] &&
                    x314[eligible_primary == TRUE, uniqueN(lesion_role), by = canonical_patient_id][, all(V1 == 2)]),
  data.table(check = "GSE237876 expression files present in RAW.tar", value = sprintf("%d/%d", sum(x237$expression_column %in% tar_list), nrow(x237)),
             pass = all(x237$expression_column %in% tar_list)),
  data.table(check = "GSE237876 all GSM matched to Table S2", value = as.character(sum(!is.na(k$tumor_type))), pass = all(!is.na(k$tumor_type))),
  data.table(check = "GSE183904 all GSM mapped to a patient", value = as.character(sum(!is.na(m183$patient))), pass = all(!is.na(m183$patient))),
  data.table(check = "GSE163558 all GSM mapped to curated patient", value = as.character(sum(!is.na(pm163$patient))), pass = all(!is.na(pm163$patient))),
  data.table(check = "GSE239676 all GSM mapped to meta Sample", value = as.character(sum(!is.na(x239$expression_column))), pass = all(!is.na(x239$expression_column))),
  data.table(check = "GSE234129 libraries are GSE239676 libraries (shared fraction >0.5 and second-best <0.1)",
             value = sprintf("min %.3f / max2 %.3f", min(ov$shared_fraction), max(ov$second_best_fraction)),
             pass = all(ov$shared_fraction > 0.5 & ov$second_best_fraction < 0.1)),
  data.table(check = "Normal peritoneum never labelled PM", value = as.character(cw[grepl("normal_peritoneum", lesion_role) & lesion_role == "solid_peritoneal_metastasis", .N]),
             pass = cw[grepl("normal_peritoneum", lesion_role), all(lesion_role != "solid_peritoneal_metastasis")]),
  data.table(check = "Ascites/lavage never labelled solid tissue", value = as.character(cw[lesion_role %in% c("malignant_ascites", "peritoneal_washings") & specimen_context == "solid_tissue", .N]),
             pass = cw[lesion_role %in% c("malignant_ascites", "peritoneal_washings"), all(specimen_context != "solid_tissue")])
))
o_inv <- pm_out("04_results", "audit", "invariants_sample_audit.csv"); write_csv_atomic(inv, o_inv)
if (!all(inv$pass)) { print(inv[pass == FALSE]); stop("Sample-audit invariants failed") }
print(cnt); print(flow)
step_end(outputs = c(o_cw, o_ov, o_ov2, o_cnt, o_flow, o_inv))
