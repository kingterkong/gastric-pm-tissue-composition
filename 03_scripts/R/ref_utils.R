# Reference verification helpers: fetch PubMed records by PMID via E-utilities,
# parse bibliographic fields and publication-type flags (retraction/correction),
# and render BibTeX. Results are cached so re-runs do not re-query the network.

suppressPackageStartupMessages(library(xml2))

eutils_fetch_pubmed <- function(pmids, cache_dir) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  pmids <- unique(as.character(pmids))
  todo <- pmids[!file.exists(file.path(cache_dir, paste0(pmids, ".xml")))]
  for (chunk in split(todo, ceiling(seq_along(todo) / 150))) {
    if (!length(chunk)) next
    url <- sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&tool=gastric_pm_refs&id=%s",
                   paste(chunk, collapse = ","))
    tmp <- tempfile(fileext = ".xml")
    for (a in 1:3) {
      st <- system2("curl", c("-sS", "-f", "--connect-timeout", "30", "--max-time", "120", "-o", shQuote(tmp), shQuote(url)))
      if (identical(st, 0L)) break
      Sys.sleep(2 * a)
    }
    if (!identical(st, 0L)) stop("E-utilities efetch failed")
    doc <- read_xml(tmp)
    for (art in xml_find_all(doc, "//PubmedArticle")) {
      id <- xml_text(xml_find_first(art, ".//MedlineCitation/PMID"))
      write_xml(art, file.path(cache_dir, paste0(id, ".xml")))
    }
    Sys.sleep(0.4)
  }
  rbindlist(lapply(pmids, function(id) {
    f <- file.path(cache_dir, paste0(id, ".xml"))
    if (!file.exists(f)) return(data.table(pmid = id, found = FALSE))
    a <- read_xml(f)
    g <- function(xp) { n <- xml_find_first(a, xp); if (inherits(n, "xml_missing")) NA_character_ else xml_text(n) }
    au <- xml_find_all(a, ".//AuthorList/Author")
    auth <- vapply(au, function(n) {
      ln <- xml_text(xml_find_first(n, "LastName")); ini <- xml_text(xml_find_first(n, "Initials"))
      cn <- xml_text(xml_find_first(n, "CollectiveName"))
      if (!is.na(ln)) paste0(ln, ", ", ifelse(is.na(ini), "", ini)) else cn
    }, "")
    ptypes <- xml_text(xml_find_all(a, ".//PublicationTypeList/PublicationType"))
    cc <- xml_find_all(a, ".//CommentsCorrectionsList/CommentsCorrections")
    cc_types <- xml_attr(cc, "RefType")
    year <- g(".//Article/Journal/JournalIssue/PubDate/Year")
    if (is.na(year)) year <- substr(g(".//Article/Journal/JournalIssue/PubDate/MedlineDate"), 1, 4)
    if (is.na(year)) year <- g(".//ArticleDate/Year")
    data.table(
      pmid = id, found = TRUE,
      title = g(".//ArticleTitle"), journal = g(".//Journal/Title"), journal_abbrev = g(".//Journal/ISOAbbreviation"),
      year = year, volume = g(".//JournalIssue/Volume"), issue = g(".//JournalIssue/Issue"),
      pages = g(".//Pagination/MedlinePgn"), elocation = g(".//ELocationID[@EIdType='pii']"),
      doi = g(".//ArticleIdList/ArticleId[@IdType='doi']") %||% g(".//ELocationID[@EIdType='doi']"),
      pmcid = g(".//ArticleIdList/ArticleId[@IdType='pmc']"),
      first_author = if (length(auth)) auth[1] else NA_character_, authors = paste(auth, collapse = " and "),
      n_authors = length(auth),
      publication_types = paste(ptypes, collapse = "; "),
      retracted = any(grepl("Retracted Publication", ptypes)) || any(cc_types %in% c("RetractionIn", "RetractionOf")),
      has_erratum = any(cc_types %in% c("ErratumIn")),
      has_expression_of_concern = any(cc_types %in% c("ExpressionOfConcernIn")),
      comments_corrections = paste(unique(cc_types), collapse = "; ")
    )
  }), fill = TRUE)
}

bibtex_escape <- function(x) {
  x <- gsub("([&%$#_])", "\\\\\\1", x)
  x
}

to_bibtex <- function(r, key) {
  f <- function(n, v) if (!is.null(v) && length(v) && !is.na(v) && nzchar(v)) sprintf("  %s = {%s},", n, bibtex_escape(v)) else NULL
  c(sprintf("@article{%s,", key),
    f("author", r$authors), f("title", sub("\\.$", "", r$title)), f("journal", r$journal),
    f("year", r$year), f("volume", r$volume), f("number", r$issue),
    f("pages", if (!is.na(r$pages)) r$pages else r$elocation), f("doi", r$doi), f("pmid", r$pmid),
    "}", "")
}
