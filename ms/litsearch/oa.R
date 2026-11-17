# Minimal OpenAlex client for the literature search ####
# Sourced by ms/litsearch/03_search.R and later steps. Every request is logged to
# ms/litsearch/search_log.csv with its hit count and cost.

library(httr2)
library(dplyr)
library(purrr)
library(stringr)
library(readr)
library(tibble)

LS <- "ms/litsearch"
CACHE <- file.path(LS, "cache")
dir.create(CACHE, showWarnings = FALSE)

## API key: Rscript --vanilla reads no .Renviron; load JV's user file explicitly ####
readRenviron(Sys.getenv(
  "LITSEARCH_RENVIRON",
  file.path(Sys.getenv("USERPROFILE"), "OneDrive - Høgskulen på Vestlandet", "Documents", ".Renviron")
))
stopifnot("openalexR.apikey not set" = nzchar(Sys.getenv("openalexR.apikey")))

OA_SELECT <- paste(
  "id", "doi", "display_name", "publication_year", "type", "language", "cited_by_count",
  "open_access", "best_oa_location", "primary_location", "abstract_inverted_index",
  "referenced_works", "related_works", "is_retracted", "authorships",
  sep = ","
)

log_query <- function(step, query, n_hits, cost = NA_real_) {
  f <- file.path(LS, "search_log.csv")
  tibble(
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    step = step, query = query, n_hits = n_hits
  ) |>
    write_csv(f, append = file.exists(f))
}

# One GET against /works; returns the parsed JSON (results + meta)
oa_get <- function(params) {
  req <- request("https://api.openalex.org/works") |>
    req_url_query(!!!params, api_key = Sys.getenv("openalexR.apikey")) |>
    req_retry(max_tries = 4, backoff = \(i) 2^i) |>
    req_throttle(rate = 8)
  resp <- req_perform(req)
  resp_body_json(resp, simplifyVector = FALSE)
}

abstract_from_index <- function(ii) {
  if (is.null(ii) || length(ii) == 0) return(NA_character_)
  pos <- unlist(ii)
  words <- rep(names(ii), lengths(ii))
  paste(words[order(pos)], collapse = " ")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

flatten_work <- function(w) {
  auth <- map_chr(w$authorships %||% list(), \(a) a$author$display_name %||% NA_character_)
  tibble(
    openalex_id = str_remove(w$id, "^https://openalex.org/"),
    doi = str_to_lower(str_remove(w$doi %||% NA_character_, "^https://doi.org/")),
    title = w$display_name %||% NA_character_,
    year = w$publication_year %||% NA_integer_,
    type = w$type %||% NA_character_,
    language = w$language %||% NA_character_,
    first_author = auth[1] %||% NA_character_,
    n_authors = length(auth),
    venue = w$primary_location$source$display_name %||% NA_character_,
    cited_by_count = w$cited_by_count %||% NA_integer_,
    is_oa = w$open_access$is_oa %||% NA,
    oa_url = w$open_access$oa_url %||% NA_character_,
    pdf_url = w$best_oa_location$pdf_url %||% NA_character_,
    is_retracted = w$is_retracted %||% NA,
    abstract = abstract_from_index(w$abstract_inverted_index),
    referenced_works = list(str_remove(unlist(w$referenced_works), "^https://openalex.org/")),
    related_works = list(str_remove(unlist(w$related_works), "^https://openalex.org/"))
  )
}

# Fetch up to `n_pages` pages of 200 for a filter string, logging each page
oa_works <- function(filter, step, label = filter, sort = NULL, n_pages = 1) {
  out <- list()
  for (p in seq_len(n_pages)) {
    params <- list(filter = filter, `per-page` = 200, page = p, select = OA_SELECT)
    if (!is.null(sort)) params$sort <- sort
    js <- tryCatch(oa_get(params), error = function(e) {
      message("OpenAlex error on ", step, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(js)) break
    res <- js$results
    log_query(step, sprintf("%s [page %d]", label, p), length(res))
    if (length(res) == 0) break
    out[[p]] <- map(res, flatten_work) |> list_rbind()
    if (p * 200 >= (js$meta$count %||% 0)) break
  }
  list_rbind(out)
}

# Fetch works by OpenAlex ID in batches of 50
oa_by_ids <- function(ids, step) {
  ids <- unique(ids[!is.na(ids)])
  split(ids, ceiling(seq_along(ids) / 50)) |>
    map(\(b) oa_works(
      paste0("openalex_id:", paste(b, collapse = "|")), step,
      label = sprintf("%d ids", length(b))
    )) |>
    list_rbind()
}
