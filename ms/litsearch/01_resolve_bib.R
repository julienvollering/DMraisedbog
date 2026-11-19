# Resolve the manuscript library (ms/Repeat-raisedbogDM.bib) to OpenAlex ####
#
# Literature search, step 2 (ms/LITSEARCH_PLAN.md §2.1). Outside R/RUNALL.R.
# Run from the repository root: Rscript --vanilla ms/litsearch/01_resolve_bib.R
#
# Writes
#   ms/litsearch/library.csv     one row per .bib entry: OpenAlex match, whether
#                                ms.qmd cites it, number of JV highlights
#   ms/litsearch/search_log.csv  appended: one row per OpenAlex call batch
#   ms/litsearch/cache/          raw OpenAlex records (gitignored)
#
# Reads JV's highlights from ms/ZoteroAnnotations/ (untracked; exists only in the
# main checkout). Set ZOTERO_ANNOT_DIR when running from a worktree.

library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tibble)

## API key ####
# Rscript --vanilla reads no .Renviron, and command-line R's home is not where
# Positron keeps it, so load the user file explicitly. The key is never printed.
renviron <- Sys.getenv(
  "LITSEARCH_RENVIRON",
  file.path(
    Sys.getenv("USERPROFILE"),
    "OneDrive - Høgskulen på Vestlandet",
    "Documents",
    ".Renviron"
  )
)
readRenviron(renviron)
stopifnot("openalexR.apikey not set" = nzchar(Sys.getenv("openalexR.apikey")))

LS <- "ms/litsearch"
CACHE <- file.path(LS, "cache")
dir.create(CACHE, showWarnings = FALSE)
ANNOT <- Sys.getenv("ZOTERO_ANNOT_DIR", "ms/ZoteroAnnotations")
stopifnot("annotation folder not found" = dir.exists(ANNOT))

log_query <- function(step, query, n_hits) {
  tibble(
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    step = step,
    query = query,
    n_hits = n_hits
  ) |>
    write_csv(
      file.path(LS, "search_log.csv"),
      append = file.exists(file.path(LS, "search_log.csv"))
    )
}

## Parse the .bib ####
# Better BibLaTeX export: one field per line, "  field = {value},".
bib_txt <- read_file("ms/Repeat-raisedbogDM.bib")
entries <- str_split(bib_txt, "\n(?=@)")[[1]] |> str_subset("^@")

field <- function(e, f) {
  m <- str_match(e, regex(paste0("\\n  ", f, " = \\{(.*?)\\},?\\n"), dotall = TRUE))
  m[, 2]
}
clean <- function(x) {
  x |>
    str_remove_all("[{}]") |>
    str_replace_all("\\\\mkbibemph", "") |>
    str_squish()
}

bib <- tibble(
  key = str_match(entries, "^@\\w+\\{([^,]+),")[, 2],
  type = str_match(entries, "^@(\\w+)\\{")[, 2],
  title = map_chr(entries, field, "title") |> clean(),
  year = str_extract(map_chr(entries, field, "date"), "^\\d{4}") |> as.integer(),
  doi = map_chr(entries, field, "doi") |> str_to_lower(),
  has_abstract = !is.na(map_chr(entries, field, "abstract"))
)
stopifnot(!anyNA(bib$key), !anyDuplicated(bib$key))

## Cited in the manuscript? ####
ms_keys <- read_file("ms/ms.qmd") |>
  str_extract_all("@[A-Za-z][A-Za-z0-9_:-]*[A-Za-z0-9]") |>
  unlist() |>
  str_remove("^@") |>
  unique()
bib <- bib |> mutate(cited_in_ms = key %in% ms_keys)

## JV's highlights ####
n_highlights <- function(k) {
  f <- file.path(ANNOT, paste0(k, ".md"))
  if (!file.exists(f)) return(NA_integer_)
  sum(str_detect(read_lines(f), "^>"))
}
bib <- bib |> mutate(n_highlights = map_int(key, n_highlights))

## OpenAlex: by DOI ####
sel <- c(
  "id", "doi", "display_name", "publication_year", "type", "cited_by_count",
  "open_access", "referenced_works", "related_works", "abstract_inverted_index",
  "primary_location", "is_retracted"
)
dois <- bib$doi[!is.na(bib$doi)]
by_doi <- openalexR::oa_fetch(
  entity = "works",
  doi = dois,
  options = list(select = sel),
  verbose = FALSE
)
log_query("bib_resolve_doi", sprintf("doi filter, %d DOIs", length(dois)), nrow(by_doi))

by_doi <- by_doi |>
  mutate(doi_key = str_to_lower(str_remove(doi, "^https://doi.org/"))) |>
  distinct(doi_key, .keep_all = TRUE)

## OpenAlex: by title, for entries without a DOI or with an unmatched DOI ####
title_sim <- function(a, b) {
  norm <- function(x) str_squish(str_remove_all(str_to_lower(x), "[^a-z0-9 ]"))
  1 - adist(norm(a), norm(b))[1, 1] / max(nchar(norm(a)), nchar(norm(b)))
}
need_title <- bib |> filter(is.na(doi) | !doi %in% by_doi$doi_key)
by_title <- map(seq_len(nrow(need_title)), \(i) {
  ti <- need_title$title[i]
  q <- str_sub(str_remove_all(ti, "[^[:alnum:] ]"), 1, 200)
  w <- tryCatch(
    openalexR::oa_fetch(
      entity = "works",
      title.search = q,
      options = list(select = sel, per_page = 5),
      pages = 1,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  log_query("bib_resolve_title", q, if (is.null(w)) 0L else nrow(w))
  if (is.null(w) || nrow(w) == 0) return(NULL)
  w$sim <- map_dbl(w$display_name, \(t) title_sim(ti, t))
  w <- w |> filter(sim >= 0.9) |> slice_max(sim, n = 1, with_ties = FALSE)
  if (nrow(w) == 0) return(NULL)
  w |> mutate(key = need_title$key[i])
}) |>
  compact() |>
  bind_rows()

## Assemble ####
works <- bind_rows(
  bib |>
    select(key, doi) |>
    inner_join(by_doi, by = c("doi" = "doi_key")) |>
    mutate(match = "doi"),
  by_title |> mutate(match = "title")
) |>
  distinct(key, .keep_all = TRUE)
saveRDS(works, file.path(CACHE, "library_openalex.rds"))

library_tbl <- bib |>
  left_join(
    works |>
      transmute(
        key,
        openalex_id = str_remove(id, "^https://openalex.org/"),
        match,
        oa_title = display_name,
        cited_by_count,
        is_oa,
        oa_url,
        n_refs = map_int(referenced_works, \(x) sum(!is.na(x)))
      ),
    by = "key"
  ) |>
  mutate(match = coalesce(match, "unresolved")) |>
  arrange(key)

write_csv(library_tbl, file.path(LS, "library.csv"), na = "")

cat(
  "Entries:", nrow(library_tbl),
  "| matched by DOI:", sum(library_tbl$match == "doi"),
  "| by title:", sum(library_tbl$match == "title"),
  "| unresolved:", sum(library_tbl$match == "unresolved"),
  "| cited in ms.qmd:", sum(library_tbl$cited_in_ms),
  "| with highlights:", sum(library_tbl$n_highlights > 0, na.rm = TRUE), "\n"
)

# sessionInfo ####
sessioninfo::session_info()
