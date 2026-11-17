# Literature search, step 4b: verify quotes, write evidence cards and claim matrix ####
#
# ms/LITSEARCH_PLAN.md §4.3-4.4. Outside R/RUNALL.R.
# Run from the repository root (Git Bash, so pdftotext output from 04 exists):
#   Rscript --vanilla ms/litsearch/05_evidence.R
#
# Reads  evidence_selection.csv (hand-picked quotes), caveats.csv, candidates.csv,
#        library.csv, library_tags.csv, claims.csv, cache/ (04_fulltext.R)
# Writes evidence/<key>.md, claim_matrix.md, new_refs_dois.txt,
#        evidence_verified.csv (each quote with where it was found, or dropped)

library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tidyr)

LS <- "ms/litsearch"
CACHE <- file.path(LS, "cache")
ANNOT <- Sys.getenv("ZOTERO_ANNOT_DIR", "ms/ZoteroAnnotations")

claims <- read_csv(file.path(LS, "claims.csv"), show_col_types = FALSE)
sel <- read_csv(file.path(LS, "evidence_selection.csv"), show_col_types = FALSE)
cav <- read_csv(file.path(LS, "caveats.csv"), show_col_types = FALSE)
cand <- read_csv(file.path(LS, "candidates.csv"), show_col_types = FALSE)
lib <- read_csv(file.path(LS, "library.csv"), show_col_types = FALSE)
tags <- read_csv(file.path(LS, "library_tags.csv"), show_col_types = FALSE)
short <- read_csv(file.path(CACHE, "shortlist.csv"), show_col_types = FALSE)
pool <- readRDS(file.path(CACHE, "pool.rds"))

## Normalisation for matching: whitespace, hyphenation, dashes, quotes, line numbers ####
norm <- function(x) {
  x |>
    str_replace_all("-\\s*\\n\\s*", "") |>
    str_replace_all("[‐‑‒–—−]", "-") |>
    str_replace_all("[‘’‛′]", "'") |>
    str_replace_all("[“”„]", '"') |>
    str_replace_all("\\s+", " ") |>
    str_to_lower() |>
    str_trim()
}

## Sources per key: list of (location, text) ####
sources_for <- function(key) {
  out <- list()
  s <- short[short$key == key, ]
  if (nrow(s)) {
    id <- s$openalex_id
    pdf_txt <- file.path(CACHE, "text", paste0(id, ".txt"))
    if (file.exists(pdf_txt)) {
      pages <- str_split(read_file(pdf_txt), "\f")[[1]]
      out <- c(out, imap(pages, \(p, i) list(loc = paste0("p. ", i, " (OA PDF page)"), text = p)))
    }
    epmc <- file.path(CACHE, "text", paste0(id, ".epmc.txt"))
    if (file.exists(epmc)) out <- c(out, list(list(loc = "Europe PMC full text (no page numbers)", text = read_file(epmc))))
    ab <- pool$abstract[pool$openalex_id == id]
    if (length(ab) && !is.na(ab)) out <- c(out, list(list(loc = "abstract", text = ab)))
  }
  f <- file.path(ANNOT, paste0(key, ".md"))
  if (file.exists(f)) {
    h <- str_subset(read_lines(f), "^>")
    out <- c(out, map(h, \(x) list(
      loc = paste0("p. ", coalesce(str_match(x, "\\*\\(p\\. ?([0-9]+)")[, 2], "?"), " (JV highlight)"),
      text = x
    )))
  }
  out
}

locate <- function(key, quote) {
  q <- norm(quote)
  # PDF layout can insert line numbers or column breaks; also try without digits-only tokens
  strip_num <- function(x) str_squish(str_remove_all(x, "\\b\\d{1,3}\\b"))
  for (s in sources_for(key)) {
    t <- norm(s$text)
    if (str_detect(t, fixed(q))) return(s$loc)
    if (str_detect(strip_num(t), fixed(strip_num(q)))) return(paste(s$loc, "[line numbers removed]"))
  }
  NA_character_
}

ver <- sel |>
  mutate(found_at = map2_chr(key, quote, locate), verified = !is.na(found_at))
write_csv(ver, file.path(LS, "evidence_verified.csv"), na = "")
cat("Quotes verified:", sum(ver$verified), "of", nrow(ver), "\n")
if (any(!ver$verified)) print(ver |> filter(!verified) |> select(key, claim_id), n = Inf)
ok <- ver |> filter(verified)

## Reference metadata (new records from OpenAlex, library from library.csv) ####
meta_new <- cand |>
  filter(relevance >= 2) |>
  transmute(key, title, year, first_author, venue, doi, in_bib = FALSE, relevance, why = reason, claims_all = claims)
meta_lib <- lib |>
  transmute(key, title, year, first_author = NA_character_, venue = NA_character_, doi, in_bib = TRUE,
            relevance = NA_integer_, why = NA_character_, claims_all = NA_character_) |>
  left_join(
    tags |> group_by(key) |> summarise(why2 = first(note), claims2 = paste(claim_id, collapse = ";"), rel2 = max(relevance)),
    by = "key"
  ) |>
  mutate(why = why2, claims_all = claims2, relevance = rel2) |>
  select(-why2, -claims2, -rel2)
meta <- bind_rows(meta_new, meta_lib) |> distinct(key, .keep_all = TRUE)

src_label <- function(key) {
  s <- short$source[short$key == key]
  if (length(s)) return(s)
  if (key %in% lib$key) return("JV's Zotero highlights")
  "unknown"
}

## Evidence cards ####
dir.create(file.path(LS, "evidence"), showWarnings = FALSE)
unlink(list.files(file.path(LS, "evidence"), full.names = TRUE))
card_keys <- union(ok$key, meta_new$key)
for (k in card_keys) {
  m <- meta[meta$key == k, ][1, ]
  q <- ok |> filter(key == k)
  cv <- cav$caveat[cav$key == k]
  src <- src_label(k)
  lines <- c(
    sprintf("# %s", k),
    "",
    sprintf("**%s** (%s)%s", coalesce(m$title, "?"), coalesce(as.character(m$year), "?"),
            ifelse(is.na(m$venue), "", paste0(". *", m$venue, "*"))),
    "",
    paste0(
      ifelse(is.na(m$first_author), "", paste0(m$first_author, " et al. · ")),
      ifelse(is.na(m$doi), "no DOI", paste0("DOI: ", m$doi)), " · ",
      ifelse(m$in_bib, "**already in the .bib**", "new (not in the .bib)"),
      ifelse(is.na(m$relevance), "", paste0(" · relevance ", m$relevance))
    ),
    "",
    sprintf("**Why it matters:** %s (claims: %s).", coalesce(m$why, "see passages"), coalesce(m$claims_all, paste(unique(q$claim_id), collapse = ";"))),
    "",
    sprintf("**Text read:** %s.", src),
    ""
  )
  if (nrow(q)) {
    lines <- c(lines, "## Passages", "")
    for (i in seq_len(nrow(q))) {
      lines <- c(
        lines,
        sprintf("> %s", q$quote[i]),
        "",
        sprintf("— %s · claim **%s** · %s · *%s*", q$found_at[i], q$claim_id[i], q$role[i], q$use[i]),
        ""
      )
    }
  } else {
    lines <- c(lines, "*No passage selected: abstract missing or not informative; screened on title and venue.*", "")
  }
  lines <- c(lines, sprintf("**Caveat:** %s", if (length(cv)) cv else "none noted beyond the study system above."))
  write_lines(lines, file.path(LS, "evidence", paste0(k, ".md")))
}

## Claim matrix ####
first_quote <- ok |> group_by(claim_id, key) |> slice(1) |> ungroup()
mx <- c(
  "# Claim → evidence matrix",
  "",
  sprintf("Generated by `ms/litsearch/05_evidence.R` on %s. One row per reference and claim, with its", Sys.Date()),
  "best verified passage. *New* = not yet in the `.bib`. Library references tagged to a claim but",
  "without a selected passage are listed after each table. Cards: `ms/litsearch/evidence/<key>.md`.",
  ""
)
for (i in seq_len(nrow(claims))) {
  cid <- claims$claim_id[i]
  d <- first_quote |> filter(claim_id == cid) |> left_join(meta |> select(key, in_bib, year), by = "key")
  lib_only <- tags |> filter(claim_id == cid, relevance >= 2, !key %in% d$key) |> pull(key)
  mx <- c(mx, sprintf("## %s — %s", cid, claims$claim[i]), "")
  if (cid == "M6") {
    m6 <- read_csv(file.path(LS, "m6_citations.csv"), show_col_types = FALSE)
    mx <- c(mx, "| Item | Citation | DOI | Status |", "|---|---|---|---|",
            sprintf("| %s | %s | %s | %s |", m6$item, m6$citation, coalesce(m6$doi, ""), m6$status), "")
    next
  }
  if (nrow(d)) {
    d <- d |> arrange(match(role, c("support", "comparison", "method", "contrast", "data")), desc(year))
    mx <- c(mx, "| Reference | Role | Passage | Where |", "|---|---|---|---|")
    mx <- c(mx, sprintf(
      "| `%s`%s | %s | %s | %s |",
      d$key, ifelse(d$in_bib, "", " *new*"), d$role,
      str_replace_all(str_trunc(d$quote, 220), "\\|", "/"), d$found_at
    ), "")
  } else {
    mx <- c(mx, "*No passage selected.*", "")
  }
  if (length(lib_only)) mx <- c(mx, paste0("Also in the library (relevance >= 2, no passage selected): ", paste0("`", lib_only, "`", collapse = ", ")), "")
  n_new <- sum(!d$in_bib, na.rm = TRUE)
  gap <- if (nrow(d) == 0 && !length(lib_only)) "**Gap:** nothing found in the library or the search." else
    if (nrow(d) > 0 && !any(d$role %in% c("support", "comparison", "method"))) "**Gap:** only contrasting evidence found." else NULL
  if (!is.null(gap)) mx <- c(mx, gap, "")
}
write_lines(mx, file.path(LS, "claim_matrix.md"))

## DOIs for Zotero 'Add Item by Identifier' ####
new_used <- meta_new |> filter(key %in% ok$key)
new_rest <- meta_new |> filter(!key %in% ok$key)
write_lines(
  c(
    "# New references with a verified passage in the claim matrix (paste into Zotero: Add Item by Identifier)",
    na.omit(new_used$doi),
    "",
    "# Data and software citations for the [TODO cite] markers (m6_citations.csv)",
    read_csv(file.path(LS, "m6_citations.csv"), show_col_types = FALSE) |>
      filter(!is.na(doi), !str_detect(status, "closed")) |>
      pull(doi),
    "",
    "# Further candidates (relevance >= 2, no passage selected); see candidates.csv",
    paste0("# ", na.omit(new_rest$doi))
  ),
  file.path(LS, "new_refs_dois.txt")
)
cat("Cards:", length(card_keys), "| new references with passages:", nrow(new_used),
    "| without DOI:", sum(is.na(new_used$doi)), "\n")

# sessionInfo ####
sessioninfo::session_info()
