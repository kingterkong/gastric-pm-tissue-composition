# Stage 14: assemble supplementary tables, run the delivery checks, and package.
# Checks are invariants that would catch a real error, not box-ticking.
source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
suppressPackageStartupMessages(library(openxlsx))
cfg <- pm_config()
R <- function(...) pm_in("04_results", ...)
A <- function(...) pm_in("00_admin", ...)

step_begin("14_verify_delivery", inputs = c(A("sample_crosswalk.csv"), R("manuscript", "key_numbers.json")),
           seed = cfg$project$seed + 14, uses_protocol = TRUE)
outs <- character()

# ---- supplementary tables ------------------------------------------------------------
supp <- list(
  S1_sample_crosswalk = A("sample_crosswalk.csv"),
  S2_overlap_register = A("overlap_register.csv"),
  S3_program_definitions = pm_in("02_data_processed", "frozen_programs", "program_definitions.csv"),
  S4_program_coverage = pm_in("02_data_processed", "frozen_programs", "program_coverage_by_platform.csv"),
  S5_matrix_type_audit = R("audit", "matrix_type_audit.csv"),
  S6_patient_deltas = R("primary", "GSE314812_patient_deltas.csv"),
  S7_hypothesis_registry = A("hypothesis_registry.csv"),
  S8_claims_ledger = A("claims_ledger.csv"),
  S9_number_trace = R("manuscript", "number_trace.csv"),
  S10_reference_verification = pm_in("07_references", "reference_verification.csv")
)
wb <- createWorkbook()
for (nm in names(supp)) {
  if (!file.exists(supp[[nm]])) next
  d <- fread(supp[[nm]])
  addWorksheet(wb, substr(nm, 1, 31)); writeData(wb, substr(nm, 1, 31), d); freezePane(wb, substr(nm, 1, 31), firstRow = TRUE)
}
supp_x <- pm_out("06_manuscript", "journal_submission", "supplementary_tables.xlsx")
assert_writable(supp_x); saveWorkbook(wb, supp_x, overwrite = TRUE); outs <- c(outs, supp_x)

# ---- invariant checks -------------------------------------------------------------
chk <- list()
add <- function(item, pass, evidence) chk[[length(chk) + 1]] <<- data.table(check = item, status = ifelse(pass, "PASS", "FAIL"), evidence = evidence)

cw <- fread(A("sample_crosswalk.csv"), colClasses = "character")
D <- fread(R("primary", "GSE314812_patient_deltas.csv"))
prim <- fread(R("primary", "primary_endpoint.csv"))
sc <- fread(R("primary", "GSE314812_program_scores_by_specimen.csv"))
K <- jsonlite::read_json(R("manuscript", "key_numbers.json"))

add("one delta per patient in the primary analysis", !anyDuplicated(D$canonical_patient_id),
    sprintf("%d rows, %d unique patients", nrow(D), uniqueN(D$canonical_patient_id)))
elig <- sc[eligible_primary == TRUE]
add("every primary-analysis patient has exactly one primary and one solid peritoneal specimen",
    elig[, .N, by = .(canonical_patient_id, lesion_role)][, all(N == 1)] &&
      elig[, uniqueN(lesion_role), by = canonical_patient_id][, all(V1 == 2)],
    sprintf("%d specimens, %d patients", nrow(elig), uniqueN(elig$canonical_patient_id)))
add("primary-analysis specimens are solid tissue only", all(elig$lesion_role %in% c("primary_tumor", "solid_peritoneal_metastasis")),
    paste(unique(elig$lesion_role), collapse = "; "))
add("ascites and lavage never coded as solid tissue",
    cw[lesion_role %in% c("malignant_ascites", "peritoneal_washings"), all(specimen_context != "solid_tissue")],
    sprintf("%d fluid specimens", cw[lesion_role %in% c("malignant_ascites", "peritoneal_washings"), .N]))
add("normal peritoneum never coded as peritoneal metastasis",
    cw[grepl("normal_peritoneum", lesion_role), .N] == 0 || cw[grepl("normal_peritoneum", lesion_role), all(lesion_role != "solid_peritoneal_metastasis")],
    sprintf("%d normal peritoneum specimens", cw[grepl("normal_peritoneum", lesion_role), .N]))
dev_set <- cw[dataset_id == "GSE183904" & grepl("calibration donor", exclusion_reason), unique(canonical_patient_id)]
val_set <- cw[eligible_validation == TRUE, unique(canonical_patient_id)]
add("no undeclared overlap between calibration and validation patient sets", length(intersect(dev_set, val_set)) == 0,
    sprintf("%d calibration donors, %d validation patients, %d shared", length(dev_set), length(val_set), length(intersect(dev_set, val_set))))
pf <- pm_in("02_data_processed", "frozen_programs", "program_definitions.json")
lockf <- A("protocol_lock.json")
lock <- jsonlite::read_json(lockf); defs <- jsonlite::read_json(pf)
hash_ok <- all(vapply(names(lock$program_definition_hashes), function(k)
  identical(lock$program_definition_hashes[[k]][[1]] %||% lock$program_definition_hashes[[k]], defs[[k]]$definition_sha256[[1]] %||% defs[[k]]$definition_sha256), TRUE))
add("frozen programme hashes match the protocol lock", hash_ok, sprintf("%d programmes", length(defs)))
add("effect direction follows the metastasis-minus-primary definition",
    abs((sc[eligible_primary == TRUE & lesion_role == "solid_peritoneal_metastasis"][order(canonical_patient_id), HRC_CORE] -
           sc[eligible_primary == TRUE & lesion_role == "primary_tumor"][order(canonical_patient_id), HRC_CORE]) -
          D[order(canonical_patient_id), HRC_CORE]) |> max() < 1e-12,
    "recomputed from specimen scores")
reg <- fread(A("hypothesis_registry.csv"))
add("every pre-specified hypothesis has a result or a stated reason",
    reg[estimable == FALSE & (is.na(not_executed_reason) | not_executed_reason == ""), .N] == 0,
    sprintf("%d hypotheses, %d not estimable and all with reasons", nrow(reg), reg[estimable == FALSE, .N]))
add("no fabricated P value for a non-estimable hypothesis",
    reg[estimable == FALSE & !is.na(p_value), .N] == 0, "checked hypothesis_registry")

# manuscript language and number checks
man <- paste(readLines(pm_in("06_manuscript", "journal_submission", "manuscript_en.md"), warn = FALSE), collapse = "\n")
banned <- c("TODO", "TBD", "FIXME", "placeholder", "please review", "needs review", "pending author review",
            "待确认", "待核查", "待补", "内部审核")
hits <- banned[vapply(banned, function(b) grepl(b, man, fixed = TRUE), TRUE)]
add("no internal placeholder language in the submission text", length(hits) == 0, if (length(hits)) paste(hits, collapse = "; ") else "none found")
overclaim <- c("independent validation", "cell-intrinsic", "peritoneal-specific", "drives metastasis",
               "therapeutic target", "predicts recurrence", "treatment benefit")
oc <- overclaim[vapply(overclaim, function(b) grepl(b, man, ignore.case = TRUE), TRUE)]
add("no prohibited claim wording in the submission text", length(oc) == 0, if (length(oc)) paste(oc, collapse = "; ") else "none found")
add("no P = 0.000 style reporting", !grepl("P\\s*=\\s*0\\.000([^0-9]|$)", man), "checked manuscript text")
trace <- fread(R("manuscript", "number_trace.csv"))
add("every traced number resolves to a result file",
    all(file.exists(file.path(pm_output_root(), unique(trace$result_file)))),
    sprintf("%d traced numbers across %d files", nrow(trace), uniqueN(trace$result_file)))
add("key manuscript numbers match the result files",
    identical(K$primary_mean[[1]], formatC(prim$mean[1], format = "f", digits = 3)) &&
      identical(K$n_primary_patients[[1]], formatC(prim$n[1], format = "d", big.mark = ",")),
    sprintf("primary mean %s, n %s", K$primary_mean[[1]], K$n_primary_patients[[1]]))
# Chinese internal and delivery documents must quote the same values as the result files.
zh_docs <- c(pm_in("06_manuscript", "internal_review", "\u6587\u7ae0\u5230\u5e95\u53d1\u73b0\u4e86\u4ec0\u4e48_\u4e2d\u6587\u5185\u90e8\u8bf4\u660e.md"),
             pm_in("06_manuscript", "internal_review", "\u5185\u90e8\u8d28\u91cf\u8bc4\u4f30\u4e0e\u53ef\u80fd\u8d28\u7591.md"),
             pm_in("08_delivery", "START_HERE_\u4e2d\u6587\u4ea4\u4ed8\u8bf4\u660e.md"),
             pm_in("README.md"), pm_in("RESUME.md"))
zh_docs <- zh_docs[file.exists(zh_docs)]
zh_keys <- c("primary_mean", "primary_ci", "adj_fib_mean", "cons_median", "implied_epi", "ext_mean")
norm_txt <- function(x) gsub("\u2212", "-", x)
zh_bad <- character()
for (f in zh_docs) {
  t <- norm_txt(paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  for (k in zh_keys) {
    v <- K[[k]][[1]]
    # intervals are written with a Chinese connector, so compare the endpoints, not the whole string
    parts <- regmatches(v, gregexpr("-?[0-9]+\\.[0-9]+", v))[[1]]
    for (pnum in parts) {
      stem <- sub("^-", "", pnum)
      if (grepl(stem, t, fixed = TRUE) && !grepl(pnum, t, fixed = TRUE))
        zh_bad <- c(zh_bad, paste0(basename(f), ":", k, ":", pnum))
    }
  }
}
add("Chinese documents quote the same values as the result files", length(zh_bad) == 0,
    if (length(zh_bad)) paste(zh_bad, collapse = "; ") else sprintf("%d documents, %d keys checked", length(zh_docs), length(zh_keys)))

# coverage numbers quoted in the text must be the scoring-universe values, not annotation coverage
covf <- fread(R("primary", "GSE314812_program_coverage.csv"))
cov_ok <- identical(K$meso_measured[[1]], as.character(covf[program_id == "COMPART_MESOTHELIUM", n_used])) &&
  identical(K$hrc_measured[[1]], as.character(covf[program_id == "HRC_CORE", n_used]))
fig3cap <- any(grepl(sprintf("%d of %d genes measured", covf[program_id == "COMPART_MESOTHELIUM", n_used],
                             covf[program_id == "COMPART_MESOTHELIUM", n_defined]),
                     readLines(pm_in("06_manuscript", "journal_submission", "manuscript_en.md"), warn = FALSE), fixed = TRUE))
add("programme coverage quoted in the text equals the scoring-universe coverage", cov_ok && fig3cap,
    sprintf("mesothelium %s, primary programme %s; figure legend consistent: %s",
            K$meso_measured[[1]], K$hrc_measured[[1]], fig3cap))

refs <- fread(pm_in("07_references", "reference_verification.csv"))
add("no retracted work cited", sum(refs$retracted) == 0, sprintf("%d references verified", nrow(refs)))
add("all cited works carry a PMID and DOI", all(nzchar(refs$pmid)) && all(!is.na(refs$doi) & nzchar(refs$doi)),
    sprintf("%d references", nrow(refs)))
# old project untouched
snap <- fread(A("old_project_snapshot.csv"), colClasses = c(path = "character", size = "numeric", mtime = "character"))
old <- pm_old_project()
f <- list.files(old, recursive = TRUE, all.files = TRUE, full.names = TRUE, no.. = TRUE)
f <- f[!grepl("/renv/(library|staging|sandbox)/", f)]
f <- f[!grepl("(^|/)\\.DS_Store$", f)]   # macOS Finder metadata; touched by the OS, carries no research content
fi <- file.info(f)
now <- data.table(path = sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", old), "/"), "", f),
                  size = as.numeric(fi$size), mtime = format(fi$mtime, "%Y-%m-%d %H:%M:%OS3"))
m <- merge(snap, now, by = "path", all = TRUE, suffixes = c(".base", ".now"))
changed <- m[is.na(size.base) | is.na(size.now) | size.base != size.now | mtime.base != mtime.now]
add("previous project left unmodified", nrow(changed) == 0,
    sprintf("%d files compared (macOS .DS_Store excluded), %d changed", nrow(m), nrow(changed)))
# figures and documents exist and are non-trivial
figs <- list.files(pm_in("05_figures", "main"), pattern = "\\.pdf$", full.names = TRUE)
add("six main figures rendered as vector PDF", length(figs) == 6 && all(file.info(figs)$size > 3000),
    sprintf("%d PDFs, smallest %d bytes", length(figs), min(file.info(figs)$size)))
docs <- list.files(pm_in("06_manuscript", "journal_submission"), pattern = "\\.(docx|pdf)$", full.names = TRUE)
add("submission documents rendered", length(docs) >= 8 && all(file.info(docs)$size > 5000),
    sprintf("%d rendered files", length(docs)))

checks <- rbindlist(chk)
o_chk <- pm_out("00_admin", "final_checklist.csv"); write_csv_atomic(checks, o_chk); outs <- c(outs, o_chk)
print(checks[, .(check = substr(check, 1, 62), status)])
if (any(checks$status == "FAIL")) { print(checks[status == "FAIL"]); stop("delivery checks failed") }
step_end(outputs = outs)
