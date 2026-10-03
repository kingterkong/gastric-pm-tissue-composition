# Stage 13a: verify every cited work against PubMed and build the bibliography.
# Citation keys used in the manuscript are declared here with the claim each supports.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
source(pm_src("03_scripts", "R", "ref_utils.R"))
cfg <- pm_config()
step_begin("13a_references", inputs = pm_src("07_references", "background_clinical_refs.csv"), seed = cfg$project$seed + 131)

CITES <- rbindlist(list(
  data.table(cite_key = "bray2024", pmid = "38572751", supports = "global gastric cancer burden"),
  data.table(cite_key = "thomassen2014", pmid = "23832847", supports = "frequency and prognosis of synchronous peritoneal metastasis"),
  data.table(cite_key = "koemans2021", pmid = "33495964", supports = "rising proportion of synchronous peritoneal metastasis; survival with and without systemic therapy"),
  data.table(cite_key = "rijken2024", pmid = "37916797", supports = "combined synchronous and metachronous peritoneal metastasis incidence"),
  data.table(cite_key = "riihimaki2016", pmid = "27447571", supports = "distribution of metastatic sites in gastric cancer"),
  data.table(cite_key = "guchelaar2025", pmid = "40026831", supports = "share of synchronous metastatic disease involving the peritoneum"),
  data.table(cite_key = "jamel2018", pmid = "28779261", supports = "prognostic value of peritoneal lavage cytology"),
  data.table(cite_key = "fidler2003", pmid = "12778135", supports = "seed and soil hypothesis"),
  data.table(cite_key = "kanda2016", pmid = "27570420", supports = "multistep model of peritoneal dissemination"),
  data.table(cite_key = "sandoval2013", pmid = "24114721", supports = "mesothelial cells as a source of carcinoma-associated fibroblasts"),
  data.table(cite_key = "janjigian2021", pmid = "34102137", supports = "first-line chemotherapy plus nivolumab in advanced gastric cancer"),
  data.table(cite_key = "charton2026", pmid = "41882239", supports = "paired primary and peritoneal metastasis transcriptomes; source of the primary cohort"),
  data.table(cite_key = "kim2024", pmid = "38498378", supports = "metastatic route-specific signals; source of the independent paired cohort; msEMT expressed by stroma"),
  data.table(cite_key = "wangq2025", pmid = "41360923", supports = "macrophage-fibroblast co-infiltration; source of solid peritoneal single-cell libraries"),
  data.table(cite_key = "liy2025", pmid = "39537239", supports = "stromal-myeloid niche in peritoneal metastasis; preadipocyte and endothelial increase in paired bulk"),
  data.table(cite_key = "cheng2024", pmid = "39097198", supports = "metastatic gastric cancer single-cell atlas; ascites specimens"),
  data.table(cite_key = "zhao2024", pmid = "39147169", supports = "spatially resolved peritoneal metastasis microenvironment"),
  data.table(cite_key = "jiang2022", pmid = "35184420", supports = "multi-site gastric cancer single-cell data used for calibration"),
  data.table(cite_key = "kumar2022", pmid = "34642171", supports = "gastric cancer single-cell atlas used for calibration"),
  data.table(cite_key = "hus2024", pmid = "39349474", supports = "ovarian and peritoneal metastasis single-cell study; malignant-cell meta-programmes"),
  data.table(cite_key = "wangr2023", pmid = "37419119", supports = "ascites and primary tumour ecosystems in gastric adenocarcinoma"),
  data.table(cite_key = "ajani2021", pmid = "32345613", supports = "YAP1 in gastric cancer peritoneal carcinomatosis"),
  data.table(cite_key = "canellas2022", pmid = "36352230", supports = "definition of the high-relapse-cell programme"),
  data.table(cite_key = "liberzon2015", pmid = "26771021", supports = "hallmark gene sets"),
  data.table(cite_key = "yoshihara2013", pmid = "24113773", supports = "ESTIMATE stromal signature"),
  data.table(cite_key = "gavish2023", pmid = "37258682", supports = "pan-cancer malignant meta-programmes"),
  data.table(cite_key = "karlsson2021", pmid = "34321199", supports = "Human Protein Atlas single-cell type expression"),
  data.table(cite_key = "tirosh2016", pmid = "27124452", supports = "expression-based inference of copy-number alterations in single cells"),
  data.table(cite_key = "hao2021", pmid = "34062119", supports = "Seurat"),
  data.table(cite_key = "germain2021", pmid = "35814628", supports = "scDblFinder"),
  data.table(cite_key = "patro2017", pmid = "28263959", supports = "Salmon transcript quantification"),
  data.table(cite_key = "soneson2015", pmid = "26925227", supports = "gene-level summarisation of transcript-level estimates"),
  data.table(cite_key = "li2011rsem", pmid = "21816040", supports = "RSEM expected counts"),
  data.table(cite_key = "tylertirosh2021", pmid = "33972543", supports = "most EMT signal in bulk tumours derives from fibroblasts")
))
stopifnot(!anyDuplicated(CITES$cite_key), !anyDuplicated(CITES$pmid))

# every citation key is checked against the fetched title before the bibliography is written
EXPECT <- c(fidler2003 = "seed and soil", li2011rsem = "RSEM", thomassen2014 = "Peritoneal carcinomatosis of gastric origin",
            sandoval2013 = "mesothelial", yoshihara2013 = "stromal and immune cell admixture", liberzon2015 = "hallmark",
            soneson2015 = "transcript-level estimates", tirosh2016 = "metastatic melanoma", kanda2016 = "peritoneal dissemination",
            riihimaki2016 = "Metastatic spread", patro2017 = "Salmon", jamel2018 = "peritoneal lavage cytology",
            ajani2021 = "YAP1", koemans2021 = "Synchronous peritoneal metastases", tylertirosh2021 = "Decoupling epithelial-mesenchymal",
            hao2021 = "multimodal single-cell", janjigian2021 = "nivolumab", karlsson2021 = "single-cell type transcriptomics",
            kumar2022 = "Single-Cell Atlas", jiang2022 = "organ-specific metastasis", germain2021 = "Doublet identification",
            canellas2022 = "EMP1", gavish2023 = "intratumour heterogeneity", wangr2023 = "ecotypes during gastric",
            rijken2024 = "Peritoneal metastases from gastric cancer", kim2024 = "metastatic route-specific",
            bray2024 = "GLOBOCAN", cheng2024 = "Metastatic Gastric Cancer", zhao2024 = "Spatially Resolved",
            hus2024 = "ovarian metastases of gastric", liy2025 = "CAF-macrophage", guchelaar2025 = "prognostic value of peritoneal",
            wangq2025 = "Multi-dimensional omics", charton2026 = "Divergent clonal evolution")

cache <- pm_out("07_references", "pubmed_cache"); dir.create(cache, recursive = TRUE, showWarnings = FALSE)
rec <- eutils_fetch_pubmed(CITES$pmid, cache)
ver <- merge(CITES, rec, by = "pmid", all.x = TRUE)
if (any(!ver$found)) stop("PubMed record not found for: ", paste(ver[found == FALSE, cite_key], collapse = ", "))
ver[, verified_date := "2026-09-26"]
ver[, pubmed_url := paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid, "/")]

mismatch <- ver[!mapply(function(k, t) grepl(EXPECT[[k]], t, fixed = TRUE), cite_key, title)]
if (nrow(mismatch)) { print(mismatch[, .(cite_key, pmid, title)]); stop("fetched title does not match the intended work") }

bad <- ver[retracted == TRUE]
if (nrow(bad)) stop("retracted work cited: ", paste(bad$cite_key, collapse = ", "))

o_ver <- pm_out("07_references", "reference_verification.csv")
write_csv_atomic(ver[, .(citation_key = cite_key, pmid, doi, title, first_author, journal = journal_abbrev, year,
                         volume, issue, pages, n_authors, publication_types, retracted, has_erratum,
                         has_expression_of_concern, comments_corrections, pubmed_url, pmcid,
                         supports_claim = supports, reading_level = "abstract or full text (see audit notes)",
                         verified_date)], o_ver)

bib <- unlist(lapply(seq_len(nrow(ver)), function(i) to_bibtex(ver[i], ver$cite_key[i])))
o_bib <- pm_out("07_references", "references.bib"); write_text_atomic(bib, o_bib)
o_map <- pm_out("07_references", "claim_citation_map.csv")
write_csv_atomic(ver[, .(citation_key = cite_key, claim_supported = supports, pmid, doi, year, first_author)], o_map)

cat(sprintf("%d references verified; %d with errata; %d retracted\n", nrow(ver), sum(ver$has_erratum), sum(ver$retracted)))
print(ver[has_erratum == TRUE, .(cite_key, title = substr(title, 1, 60), comments_corrections)])
step_end(outputs = c(o_ver, o_bib, o_map))
