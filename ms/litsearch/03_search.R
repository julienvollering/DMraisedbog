# Literature search, step 3: keyword searches + citation chaining ####
#
# ms/LITSEARCH_PLAN.md §2.2-2.3 and §3 (automatic part). Outside R/RUNALL.R.
# Run from the repository root: Rscript --vanilla ms/litsearch/03_search.R
#
# Reads  queries.csv, claims.csv, library.csv, library_tags.csv,
#        cache/library_openalex.rds (01_resolve_bib.R)
# Writes cache/pool.rds     every unique record found (OpenAlex fields)
#        cache/routes.csv   how each record was reached, per claim
#        ranked.csv         per claim x record: automatic score and rank
#        search_log.csv     appended, one row per API page

source("ms/litsearch/oa.R")

claims <- read_csv(file.path(LS, "claims.csv"), show_col_types = FALSE)
queries <- read_csv(file.path(LS, "queries.csv"), show_col_types = FALSE)
lib <- read_csv(file.path(LS, "library.csv"), show_col_types = FALSE)
tags <- read_csv(file.path(LS, "library_tags.csv"), show_col_types = FALSE)
lib_oa <- readRDS(file.path(CACHE, "library_openalex.rds")) |>
  mutate(openalex_id = str_remove(id, "^https://openalex.org/"))

## 1. Keyword searches (title + abstract), two pages of 200 per query ####
cached <- function(name, expr) {
  f <- file.path(CACHE, paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  x <- expr
  saveRDS(x, f)
  x
}

q_hits <- cached("q_hits", map(seq_len(nrow(queries)), \(i) {
  q <- queries[i, ]
  f <- sprintf(
    "title_and_abstract.search:%s,from_publication_date:%d-01-01,type:article|review|book-chapter,is_retracted:false",
    q$query, q$from_year
  )
  w <- oa_works(f, step = paste0("search_", q$query_id), label = q$query, n_pages = 2)
  if (is.null(w) || nrow(w) == 0) return(NULL)
  w |> mutate(claim_id = q$claim_id, route = "search", via = q$query_id, rank = row_number())
}) |>
  list_rbind())

## 2. Anchors: core (relevance 3) library papers per claim, data citations excluded ####
anchors <- tags |>
  filter(relevance == 3, role != "data") |>
  inner_join(lib |> filter(!is.na(openalex_id)) |> select(key, openalex_id), by = "key")

## 3. Forward chaining: most-cited citers, and most-cited recent citers ####
anchor_ids <- unique(anchors$openalex_id)
fwd <- cached("forward", map(anchor_ids, \(a) {
  all <- oa_works(paste0("cites:", a), "forward_all", label = a, sort = "cited_by_count:desc")
  recent <- oa_works(
    paste0("cites:", a, ",from_publication_date:2021-01-01"), "forward_recent",
    label = a, sort = "cited_by_count:desc"
  )
  both <- bind_rows(all, recent)
  if (nrow(both) == 0) return(NULL)
  both |>
    distinct(openalex_id, .keep_all = TRUE) |>
    mutate(anchor = a)
}) |>
  list_rbind())
fwd_routes <- fwd |>
  select(-any_of("claim_id")) |>
  inner_join(anchors |> distinct(openalex_id, claim_id), by = c("anchor" = "openalex_id"),
             relationship = "many-to-many") |>
  mutate(route = "forward", via = anchor, rank = NA_integer_)

## 4. Backward chaining and related works: keep records reached from 2+ anchors ####
link_tbl <- function(col, route) {
  lib_oa |>
    filter(openalex_id %in% anchor_ids) |>
    transmute(anchor = openalex_id, target = map({{ col }}, \(x) str_remove(x, "^https://openalex.org/"))) |>
    tidyr::unnest(target) |>
    filter(!is.na(target)) |>
    mutate(route = route)
}
links <- bind_rows(link_tbl(referenced_works, "backward"), link_tbl(related_works, "related"))
multi <- links |>
  group_by(target) |>
  filter(n_distinct(anchor) >= 2) |>
  ungroup()
bwd_works <- cached("backward", oa_by_ids(setdiff(unique(multi$target), lib_oa$openalex_id), "backward_related_fetch"))
bwd_routes <- multi |>
  inner_join(anchors |> distinct(openalex_id, claim_id), by = c("anchor" = "openalex_id"),
             relationship = "many-to-many") |>
  transmute(openalex_id = target, claim_id, route, via = anchor, rank = NA_integer_)

## Pool ####
pool <- bind_rows(q_hits, fwd, bwd_works) |>
  select(-any_of(c("claim_id", "route", "via", "rank", "anchor"))) |>
  distinct(openalex_id, .keep_all = TRUE)
saveRDS(pool, file.path(CACHE, "pool.rds"))

routes <- bind_rows(
  q_hits |> select(openalex_id, claim_id, route, via, rank),
  fwd_routes |> select(openalex_id, claim_id, route, via, rank),
  bwd_routes
) |>
  distinct()
write_csv(routes, file.path(CACHE, "routes.csv"))

## Automatic ranking per claim ####
# Keyword vocabulary per claim: every quoted phrase and bare term in its queries
kw <- queries |>
  mutate(terms = str_extract_all(query, '"[^"]+"|\\b[A-Za-z][A-Za-z0-9-]{3,}\\b')) |>
  tidyr::unnest(terms) |>
  mutate(terms = str_to_lower(str_remove_all(terms, '"'))) |>
  filter(!terms %in% c("and", "or", "not", "model", "climate", "change", "data")) |>
  distinct(claim_id, terms)

in_bib_ids <- na.omit(lib$openalex_id)
in_bib_dois <- na.omit(str_to_lower(lib$doi))

pool_f <- pool |>
  mutate(
    in_bib = openalex_id %in% in_bib_ids | (!is.na(doi) & doi %in% in_bib_dois),
    text = str_to_lower(paste(title, coalesce(abstract, ""))),
    age = pmax(1, 2027 - coalesce(year, 2000L))
  ) |>
  filter(
    !coalesce(is_retracted, FALSE),
    coalesce(language, "en") == "en",
    type %in% c("article", "review", "book-chapter", "book", "preprint"),
    !is.na(title)
  )

route_sum <- routes |>
  group_by(openalex_id, claim_id) |>
  summarise(
    q_score = max(1 - coalesce(rank, 400L) / 400, 0),
    n_anchor = n_distinct(via[route != "search"]),
    routes = paste(sort(unique(route)), collapse = "+"),
    .groups = "drop"
  )

kw_hits <- pool_f |>
  select(openalex_id, text) |>
  tidyr::crossing(claim_id = unique(kw$claim_id)) |>
  left_join(kw |> group_by(claim_id) |> summarise(terms = list(terms)), by = "claim_id") |>
  mutate(n_kw = map2_int(text, terms, \(t, k) sum(str_detect(t, fixed(k))))) |>
  select(openalex_id, claim_id, n_kw)

ranked <- route_sum |>
  inner_join(pool_f |> select(openalex_id, in_bib, cited_by_count, age, abstract), by = "openalex_id") |>
  left_join(kw_hits, by = c("openalex_id", "claim_id")) |>
  mutate(
    has_abstract = !is.na(abstract),
    score = 2 * q_score + pmin(n_anchor, 4) + 0.5 * pmin(coalesce(n_kw, 0L), 4) +
      0.3 * log1p(coalesce(cited_by_count, 0L) / age)
  ) |>
  select(-abstract) |>
  group_by(claim_id) |>
  arrange(desc(score), .by_group = TRUE) |>
  mutate(rank_in_claim = row_number()) |>
  ungroup()
write_csv(ranked, file.path(LS, "ranked.csv"))

cat(
  "Pool:", nrow(pool), "records;", nrow(pool_f), "after pre-filter;",
  sum(pool_f$in_bib), "already in the .bib\n"
)
ranked |> count(claim_id, name = "n_records") |> print(n = Inf)

# sessionInfo ####
sessioninfo::session_info()
