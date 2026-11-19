# Literature search, step 4a: OA full text and candidate passages ####
#
# ms/LITSEARCH_PLAN.md §4.1-4.2. Outside R/RUNALL.R.
# Run from the repository root: Rscript --vanilla ms/litsearch/04_fulltext.R
#
# Shortlist = screened records with relevance >= 2 (screen.csv), plus library papers
# tagged relevance 3 (their passages come from JV's highlights, not from the PDF).
# Full text is OA only (JV, answer 4): pdf_url / oa_url from OpenAlex. JV's Zotero
# PDFs are never read.
#
# Writes candidates.csv            every screened record with its screening fields
#        cache/pdf/, cache/text/   downloads and page-split text (gitignored)
#        cache/passages.csv        candidate passages per paper (id, page, text)

source("ms/litsearch/oa.R")

claims <- read_csv(file.path(LS, "claims.csv"), show_col_types = FALSE)
queries <- read_csv(file.path(LS, "queries.csv"), show_col_types = FALSE)
lib <- read_csv(file.path(LS, "library.csv"), show_col_types = FALSE)
tags <- read_csv(file.path(LS, "library_tags.csv"), show_col_types = FALSE)
pool <- readRDS(file.path(CACHE, "pool.rds"))
ranked <- read_csv(file.path(LS, "ranked.csv"), show_col_types = FALSE)
screen <- read_csv(file.path(LS, "screen.csv"), show_col_types = FALSE)
screened_ids <- read_csv(file.path(CACHE, "screen_ids.csv"), show_col_types = FALSE)
ANNOT <- Sys.getenv("ZOTERO_ANNOT_DIR", "ms/ZoteroAnnotations")

## Citation keys for new records, in the Better BibTeX style of the library ####
ascii <- function(x) iconv(x, to = "ASCII//TRANSLIT") |> str_remove_all("[^A-Za-z]")
stop_w <- c("a", "an", "the", "of", "and", "or", "in", "on", "for", "to", "with", "from",
            "by", "at", "as", "is", "are", "do", "does", "can", "how", "what", "when", "why",
            "will", "be", "into", "under", "its", "their", "using", "use", "toward", "towards")
make_key <- function(author, title, year) {
  sur <- str_to_lower(ascii(str_extract(author, "[^ ]+$")))
  words <- str_split(str_to_lower(str_remove_all(title, "[^[:alnum:] -]")), "[ -]+")[[1]]
  words <- words[!words %in% stop_w & nzchar(words)][1:3]
  paste0(sur, paste0(str_to_title(ascii(words[!is.na(words)])), collapse = ""), year)
}

## candidates.csv: every screened record ####
candidates <- screened_ids |>
  left_join(screen, by = "openalex_id") |>
  left_join(
    pool |> select(openalex_id, doi, title, year, first_author, venue, cited_by_count, is_oa, oa_url, pdf_url, abstract),
    by = "openalex_id"
  ) |>
  left_join(
    ranked |> group_by(openalex_id) |> summarise(best_score = max(score), routes = paste(unique(routes), collapse = "|")),
    by = "openalex_id"
  ) |>
  mutate(
    relevance = coalesce(relevance, 0L),
    reason = coalesce(reason, "off-topic for the claim register"),
    abstract_missing = is.na(abstract),
    key = pmap_chr(list(first_author, title, year), make_key)
  ) |>
  rename(top_ranked_for = claims.x, claims = claims.y)
stopifnot(!anyDuplicated(candidates$key[candidates$relevance >= 2]))
candidates |>
  select(openalex_id, key, doi, title, year, first_author, venue, cited_by_count, is_oa,
         top_ranked_for, relevance, claims, role, reason, abstract_missing, best_score, routes) |>
  arrange(desc(relevance), claims, desc(best_score)) |>
  write_csv(file.path(LS, "candidates.csv"), na = "")

## Fetch OA full text for the new shortlist ####
short <- candidates |> filter(relevance >= 2)
dir.create(file.path(CACHE, "pdf"), showWarnings = FALSE)
dir.create(file.path(CACHE, "text"), showWarnings = FALSE)

get_raw <- function(u) {
  tryCatch({
    resp <- request(u) |>
      req_user_agent("Mozilla/5.0 (literature search; R httr2)") |>
      req_timeout(60) |>
      req_error(is_error = \(r) FALSE) |>
      req_perform()
    if (resp_status(resp) == 200) resp_body_raw(resp) else NULL
  }, error = function(e) NULL)
}
is_pdf_raw <- function(r) !is.null(r) && length(r) > 4 && identical(r[1:4], charToRaw("%PDF"))

fetch_pdf <- function(id, urls) {
  out <- file.path(CACHE, "pdf", paste0(id, ".pdf"))
  if (file.exists(out)) return(out)
  for (u in unique(na.omit(urls))) {
    r <- get_raw(u)
    if (!is_pdf_raw(r) && !is.null(r)) {
      # landing page: follow the scholarly citation_pdf_url meta tag if present
      html <- tryCatch(rawToChar(r), error = function(e) "")
      m <- str_match(html, '<meta[^>]+name="citation_pdf_url"[^>]+content="([^"]+)"')[, 2]
      if (is.na(m)) m <- str_match(html, '<meta[^>]+content="([^"]+)"[^>]+name="citation_pdf_url"')[, 2]
      r <- if (!is.na(m)) get_raw(m) else NULL
    }
    if (is_pdf_raw(r)) {
      writeBin(r, out)
      return(out)
    }
  }
  NA_character_
}

# Europe PMC open-access full text (XML), for records without a PDF
fetch_epmc <- function(id, doi) {
  out <- file.path(CACHE, "text", paste0(id, ".epmc.txt"))
  if (file.exists(out)) return(out)
  if (is.na(doi)) return(NA_character_)
  s <- tryCatch(
    request("https://www.ebi.ac.uk/europepmc/webservices/rest/search") |>
      req_url_query(query = sprintf('DOI:"%s"', doi), format = "json", resultType = "lite") |>
      req_perform() |>
      resp_body_json(),
    error = function(e) NULL
  )
  hit <- s$resultList$result
  if (!length(hit) || is.null(hit[[1]]$pmcid) || !identical(hit[[1]]$isOpenAccess, "Y")) return(NA_character_)
  r <- get_raw(sprintf("https://www.ebi.ac.uk/europepmc/webservices/rest/%s/fullTextXML", hit[[1]]$pmcid))
  if (is.null(r)) return(NA_character_)
  x <- rawToChar(r)
  x <- str_remove(x, "(?s)<ref-list.*$")
  paras <- str_match_all(x, "(?s)<p[ >].*?</p>")[[1]][, 1] |>
    str_remove_all("<[^>]+>") |>
    str_squish()
  write_lines(paras, out)
  out
}

pdf_pages <- function(pdf) {
  txt <- file.path(CACHE, "text", sub("\\.pdf$", ".txt", basename(pdf)))
  if (!file.exists(txt)) system2("pdftotext", c("-layout", shQuote(pdf), shQuote(txt)), stdout = FALSE, stderr = FALSE)
  if (!file.exists(txt)) return(character())
  str_split(read_file(txt), "\f")[[1]]
}

short$pdf <- map2_chr(short$openalex_id, map2(short$pdf_url, short$oa_url, c), fetch_pdf)
short$epmc <- ifelse(is.na(short$pdf), map2_chr(short$openalex_id, short$doi, fetch_epmc), NA_character_)
short$source <- case_when(
  !is.na(short$pdf) ~ "OA full text (PDF)",
  !is.na(short$epmc) ~ "OA full text (Europe PMC XML, no page numbers)",
  TRUE ~ "abstract only"
)

## Candidate passages ####
# Sentences scored by the claim vocabulary of the paper's claims; references cut off
kw <- queries |>
  mutate(terms = str_extract_all(query, '"[^"]+"|\\b[A-Za-z][A-Za-z0-9-]{3,}\\b')) |>
  tidyr::unnest(terms) |>
  mutate(terms = str_to_lower(str_remove_all(terms, '"'))) |>
  filter(!terms %in% c("and", "or", "not", "model", "data")) |>
  distinct(claim_id, terms)
extra <- c("warm", "temperature", "climat", "extrapolat", "persist", "suitab", "project",
           "loss", "bog", "peat", "sphagnum", "transfer", "uncertain", "novel", "prevalence")

sentences <- function(txt) {
  txt |>
    str_replace_all("-\\s*\\n\\s*", "") |>
    str_replace_all("\\s+", " ") |>
    str_split("(?<=[.!?])\\s+(?=[A-Z])") |>
    unlist()
}
score_sent <- function(s, terms) {
  s_l <- str_to_lower(s)
  map_int(s_l, \(x) sum(str_detect(x, fixed(terms)))) + 0.5 * map_int(s_l, \(x) sum(str_detect(x, fixed(extra))))
}

cand_passages <- function(id, claims_str, pdf, epmc, abstract, n = 6) {
  terms <- kw$terms[kw$claim_id %in% str_split(claims_str, ";")[[1]]]
  if (is.na(pdf) && !is.na(epmc)) {
    d <- tibble(page = NA_integer_, text = sentences(paste(read_lines(epmc), collapse = " ")))
  } else if (!is.na(pdf)) {
    pages <- pdf_pages(pdf)
    ref_page <- which(str_detect(pages, "(?m)^\\s*(References|REFERENCES|Literature [Cc]ited|Bibliography)\\s*$"))
    if (length(ref_page)) pages <- pages[seq_len(min(ref_page))]
    d <- imap(pages, \(p, i) tibble(page = i, text = sentences(p))) |> list_rbind()
    # scanned or unreadable PDF: fall back to the abstract
    if (is.null(d) || sum(nchar(d$text)) < 2000) {
      if (is.na(abstract)) return(NULL)
      d <- tibble(page = NA_integer_, text = sentences(abstract))
    }
  } else if (!is.na(abstract)) {
    d <- tibble(page = NA_integer_, text = sentences(abstract))
  } else {
    return(NULL)
  }
  d |>
    filter(nchar(text) >= 80, nchar(text) <= 600, !str_detect(text, "https?://|doi\\.org|©|\\bet al\\.,? \\(?\\d{4}")) |>
    mutate(score = score_sent(text, terms)) |>
    filter(score > 0) |>
    slice_max(score, n = n, with_ties = FALSE) |>
    mutate(openalex_id = id)
}
passages_new <- pmap(
  list(short$openalex_id, short$claims, short$pdf, short$epmc, short$abstract),
  cand_passages
) |>
  list_rbind() |>
  left_join(short |> select(openalex_id, key, source), by = "openalex_id") |>
  mutate(origin = source)

## Library papers: candidate passages from JV's highlights ####
lib_short <- tags |>
  filter(relevance == 3, role != "data") |>
  group_by(key) |>
  summarise(claims = paste(claim_id, collapse = ";"), .groups = "drop")
hl_passages <- map2(lib_short$key, lib_short$claims, \(k, cl) {
  f <- file.path(ANNOT, paste0(k, ".md"))
  if (!file.exists(f)) return(NULL)
  h <- str_subset(read_lines(f), "^>")
  if (!length(h)) return(NULL)
  tibble(
    raw = h,
    page = as.integer(str_match(h, "\\*\\(p\\. ?(\\d+)")[, 2]),
    text = h |> str_remove("^> ?") |> str_remove("\\*\\*\\[[^]]+\\]\\*\\* ?") |>
      str_remove("\\s*\\*\\(p\\..*$") |> str_squish()
  ) |>
    filter(nchar(text) >= 40) |>
    mutate(score = score_sent(text, kw$terms[kw$claim_id %in% str_split(cl, ";")[[1]]])) |>
    slice_max(score, n = 5, with_ties = FALSE) |>
    mutate(key = k, origin = "JV highlight", source = "JV highlight", openalex_id = NA_character_)
}) |>
  list_rbind()

passages <- bind_rows(passages_new, hl_passages |> select(-raw)) |>
  group_by(key) |>
  mutate(pid = paste0(key, "#", row_number())) |>
  ungroup()
write_csv(passages, file.path(CACHE, "passages.csv"))
write_csv(short |> select(openalex_id, key, claims, source), file.path(CACHE, "shortlist.csv"))

cat(
  "Shortlist:", nrow(short), "new records;", sum(!is.na(short$pdf)), "OA PDF;",
  sum(!is.na(short$epmc)), "Europe PMC XML;",
  nrow(lib_short), "library papers via highlights;", nrow(passages), "candidate passages\n"
)

# sessionInfo ####
sessioninfo::session_info()
