# How generous is the bog map, and what counts as "no longer suitable"? ####

# PURPOSE: Reads the production projection against the Lyngstad survey footprint, where bog prevalence is known, to give a real precision for the argmax bog map, a prior-shifted (calibrated) P(bog), and three competing criteria for loss under each scenario.

# The production forest is trained at equal class priors (make.wt / make.size), so argmax
# answers "which class is most likely at equal priors" and over-predicts bog against real
# prevalence; balanced P(bog) is a likelihood ratio, not an occupancy chance (notebook
# 2026-10-07 (b), plan_calibration_check.md). The OOB precision of 0.19 is a property of
# the designed training prevalence, not of the map. Inside the Lyngstad survey footprint
# the mapping is complete, so prevalence is known there and three things can be measured
# from rasters that already exist (no refit, no new prediction):
#
#  A. Commission. Of the cells the present-day map calls bog inside the footprint, how
#     many are mapped bog? Split by what the others are (mapped other-peat vs non-peat).
#  B. Calibrated P(bog). Prior shift of the three balanced probabilities to the
#     footprint's class shares, and a reliability diagram against mapped bog.
#  C. Loss criteria. Class flip (the current pipeline), falling below a percentile of
#     balanced P(bog) at present-day mapped bogs, and calibrated odds relative to present
#     bogs -- per scenario, with where they disagree.
#
# IN-SAMPLE CAVEAT. Every Norwegian bog cell, and a subsample of the footprint's absence
# cells, is a training row of the production forest, so the projected probability there is
# in-sample and flatters recall. At those cells the OOB probability from the production fit
# is substituted for the current-climate value ("oob" below); the plain projection is kept
# alongside ("projection") so the size of the effect is visible. Future values are always
# the projection: no training row carries future climate.
#
# REFERENCE FIT. The production forest draws 3 x n_bog rows per tree without replacement
# at equal class weight, so each tree takes ~98% of the bog rows and a bog row is OOB in
# only ~1.5% of its 1,000 trees. OOB P(bog) at the mapped bogs then rests on a handful of
# trees, and the 10th percentile used as the strict loss cut is coarse (0.556 = 5/9 on
# 2026-10-07). The OOB values here come instead from a sidecar fit: the same call as
# pl2_predict.R (same frame, same make.wt / make.size / swor recipe) with many more trees,
# so it is the same estimator with the Monte-Carlo noise taken out. It is used only for
# OOB values at training cells; the projection is still the production forest's.
# plan_calibration_check.md, "Note 2026-11-19".

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(sf)
library(terra)
library(randomForestSRC)

source("R/functions.R")
source("R/config.R")

# Band edges for artype_60 -- must match pl0_labelNorwayBlock.R, which defines what
# "otherpeat" and "nonpeat" mean in the Norway block. Cells in between were excluded from
# training ("mixed"); cells without artype_60 are "unlabelled".
OTHERPEAT_MIN <- 0.5
NONPEAT_MAX <- 0.1

# Reference percentile of balanced P(bog) at present-day mapped bogs below which a cell
# counts as lost (criterion 2), and the share of reference odds a cell must keep under the
# calibrated criterion (criterion 3).
REF_PERCENTILE <- 0.10
ODDS_RETAINED <- 0.5

# Thresholds on balanced P(bog) for the context sweep in step A.
SWEEP <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)

# Bins of balanced P(bog) for the reliability diagram.
REL_BREAKS <- seq(0, 1, by = 0.1)

# Cuts on balanced P(bog) for the persistence sweep in step C: how the share of present
# bog that stays suitable depends on where the strict cut is put.
PERSIST_CUTS <- seq(0.30, 0.70, by = 0.05)

# Sidecar reference fit (see REFERENCE FIT above). At ~1.5% OOB per bog row, 13,000 trees
# leave each bog row OOB in ~200. The diagnostic fit counts the actual OOB rate from the
# in-bag record, which the production fit did not keep.
REF_NTREE <- 13000
REF_SEED <- -4242
DIAG_NTREE <- 200

record_settings(
  "R/pl2_calibrationCheck.R",
  otherpeat_min_artype60 = OTHERPEAT_MIN,
  nonpeat_max_artype60 = NONPEAT_MAX,
  ref_percentile = REF_PERCENTILE,
  odds_retained = ODDS_RETAINED,
  threshold_sweep = SWEEP,
  reliability_breaks = REL_BREAKS,
  persistence_cuts = PERSIST_CUTS,
  prior_source = "Lyngstad footprint: bog share; non-bog split otherpeat:nonpeat by labelled counts",
  current_probabilities = "OOB of the sidecar reference fit at Norwegian training cells, production projection elsewhere",
  reference_ntree = REF_NTREE,
  reference_seed = REF_SEED,
  diagnostic_ntree = DIAG_NTREE
)

## Load predictions ####

SCEN <- tibble(
  scenario = c("current", FUTURE_SCENARIO, MILD_SCENARIO),
  prefix = c(
    "current_", "future_", paste0("future", scenario_suffix(MILD_SCENARIO), "_")
  )
)

pred <- c(
  rast("output/pl2/rf_local_pred_brf.tif"),
  rast(paste0("output/pl2/rf_local_pred_brf", scenario_suffix(MILD_SCENARIO), ".tif"))
)
stopifnot(all(unlist(map(SCEN$prefix, ~ paste0(.x, RESPONSE_LEVELS))) %in% names(pred)))
cell_km2 <- prod(res(pred)) / 1e6

# Full-domain probabilities, one matrix per scenario, rows = cells with a prediction.
P_all <- values(pred)
ok <- which(complete.cases(P_all))
P_all <- P_all[ok, , drop = FALSE]
scen_mat <- function(M, prefix) {
  out <- M[, paste0(prefix, RESPONSE_LEVELS), drop = FALSE]
  colnames(out) <- RESPONSE_LEVELS
  out
}
P <- set_names(map(SCEN$prefix, ~ scen_mat(P_all, .x)), SCEN$scenario)
rm(P_all)
row_of_cell <- rep(NA_integer_, ncell(pred))
row_of_cell[ok] <- seq_along(ok)

argmax <- function(M) {
  factor(RESPONSE_LEVELS[max.col(M, ties.method = "first")], levels = RESPONSE_LEVELS)
}

## OOB substitution at the Norwegian training cells ####

mf <- read_csv("output/pl2/modeling_frame_regional.csv", show_col_types = FALSE) |>
  filter(scenario == "current")
model <- readRDS("output/pl2/model_production_brf.rds")
oob_prod <- model$predicted.oob[, RESPONSE_LEVELS]
prod_xvars <- model$xvar.names
# pl2_predict.R fits on the current rows of the frame in file order, so OOB rows line up.
stopifnot(
  nrow(oob_prod) == nrow(mf),
  identical(as.character(model$yvar), mf$response)
)
rm(model)
gc()

# The training frame exactly as pl2_predict.R builds it.
train_df <- mf |>
  select(-scenario, -x, -y) |>
  mutate(response = factor(response, levels = RESPONSE_LEVELS)) |>
  as.data.frame()
predictor_names <- setdiff(names(train_df), c("response", "dataset"))
train_df <- train_df[, c("response", predictor_names)]
stopifnot(setequal(predictor_names, prod_xvars))

fit_reference <- function(ntree, seed, membership = FALSE) {
  rfsrc(
    formula = response ~ .,
    data = train_df,
    ntree = ntree,
    case.wt = randomForestSRC:::make.wt(train_df$response),
    sampsize = randomForestSRC:::make.size(train_df$response),
    importance = FALSE,
    forest = FALSE,
    membership = membership,
    seed = seed
  )
}

# Diagnostic: share of trees in which a row is OOB, by class, from the in-bag record.
diag_fit <- fit_reference(DIAG_NTREE, REF_SEED - 1, membership = TRUE)
oob_rate <- rowMeans(diag_fit$inbag == 0)
rm(diag_fit)
gc()

t_ref <- Sys.time()
ref_fit <- fit_reference(REF_NTREE, REF_SEED)
cat(
  "Reference fit,", REF_NTREE, "trees:",
  round(as.numeric(difftime(Sys.time(), t_ref, units = "mins")), 1), "min\n"
)
oob <- ref_fit$predicted.oob[, RESPONSE_LEVELS]
rm(ref_fit)
gc()
stopifnot(nrow(oob) == nrow(mf), !anyNA(oob))

oob_trees <- tibble(response = train_df$response, oob_rate = oob_rate) |>
  group_by(response) |>
  summarise(
    n_rows = n(),
    oob_rate_mean = mean(oob_rate),
    oob_rate_min = min(oob_rate),
    oob_trees_production = oob_rate_mean * 1000,
    oob_trees_reference = oob_rate_mean * REF_NTREE,
    .groups = "drop"
  )
cat("\nShare of trees in which a training row is OOB (diagnostic fit):\n")
print(as.data.frame(oob_trees), row.names = FALSE, digits = 3)

no_rows <- which(mf$dataset == "NO")
train_cell <- cellFromXY(pred, as.matrix(mf[no_rows, c("x", "y")]))
stopifnot(
  !anyNA(train_cell),
  !anyDuplicated(train_cell),
  max(abs(xyFromCell(pred, train_cell) - as.matrix(mf[no_rows, c("x", "y")]))) < 1e-6
)
train_row <- row_of_cell[train_cell]
stopifnot(!anyNA(train_row))

P_cur_oob <- P[["current"]]
P_cur_oob[train_row, ] <- oob[no_rows, ]
cat("OOB probabilities substituted at", length(train_row), "Norwegian training cells\n")

## Footprint population and mapped class ####

# The footprint is the presence + absence population of pl0_rasterizePresenceAbsence.R:
# presence = cells whose centre falls in a Lyngstad raised-bog polygon, absence = the rest
# of the NiB survey footprint (cells whose centre falls in it). A cell partly covered by a
# bog polygon whose centre is outside counts as absence -- see the comment there.
presence <- read_csv("output/presence_coords_regional.csv", show_col_types = FALSE)
absence <- read_csv("output/absence_coords_regional.csv", show_col_types = FALSE)
fp_xy <- bind_rows(
  presence |> mutate(mapped_bog = TRUE),
  absence |> mutate(mapped_bog = FALSE)
)
fp_cell <- cellFromXY(pred, as.matrix(fp_xy[c("x", "y")]))
stopifnot(!anyNA(fp_cell), !anyDuplicated(fp_cell))

preds_reg <- rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")
artype <- terra::extract(preds_reg[["artype_60"]], as.matrix(fp_xy[c("x", "y")]))[[1]]

fp <- fp_xy |>
  mutate(
    cell = fp_cell,
    row = row_of_cell[fp_cell],
    artype_60 = artype,
    mapped = factor(
      case_when(
        mapped_bog ~ "bog",
        is.na(artype_60) ~ "unlabelled",
        artype_60 >= OTHERPEAT_MIN ~ "otherpeat",
        artype_60 < NONPEAT_MAX ~ "nonpeat",
        TRUE ~ "mixed"
      ),
      levels = c(RESPONSE_LEVELS, "mixed", "unlabelled")
    )
  )
cat(
  "Footprint cells:", nrow(fp), "| without a prediction:", sum(is.na(fp$row)), "\n"
)
fp <- fp |> filter(!is.na(row))

# The land-use / water screen the training rows passed. The map itself is applied to every
# cell, so the unscreened footprint is the main population; the screened one shows how much
# of the commission sits on cells the training data never sampled.
fp$screened <- landuse_screen_keep(fp, label = "footprint")

cat("\nFootprint cells by mapped class:\n")
print(as.data.frame(count(fp, mapped)), row.names = FALSE)

## A. Commission inside the footprint ####

pbog_fp <- list(
  oob = P_cur_oob[fp$row, "bog"],
  projection = P[["current"]][fp$row, "bog"]
)
cls_fp <- list(
  oob = argmax(P_cur_oob[fp$row, , drop = FALSE]),
  projection = argmax(P[["current"]][fp$row, , drop = FALSE])
)

score <- function(called, mapped) {
  tp <- sum(called & mapped == "bog")
  n_called <- sum(called)
  n_mapped <- sum(mapped == "bog")
  tibble(
    n_cells = length(mapped),
    n_mapped_bog = n_mapped,
    prevalence = n_mapped / length(mapped),
    n_pred_bog = n_called,
    tp = tp,
    recall = tp / n_mapped,
    precision = tp / n_called,
    commission = 1 - tp / n_called,
    pred_to_mapped_ratio = n_called / n_mapped,
    on_otherpeat = sum(called & mapped == "otherpeat"),
    on_nonpeat = sum(called & mapped == "nonpeat"),
    on_mixed = sum(called & mapped == "mixed"),
    on_unlabelled = sum(called & mapped == "unlabelled")
  )
}

commission <- expand_grid(
  population = c("footprint", "footprint, screened"),
  probabilities = c("oob", "projection")
) |>
  pmap_dfr(function(population, probabilities) {
    keep <- if (population == "footprint") rep(TRUE, nrow(fp)) else fp$screened
    pb <- pbog_fp[[probabilities]][keep]
    mapped <- fp$mapped[keep]
    n_mapped <- sum(mapped == "bog")
    # Area-matched cut: the balanced P(bog) above which exactly as many cells are called
    # bog as are mapped bog.
    cut_am <- sort(pb, decreasing = TRUE)[n_mapped]
    bind_rows(
      score(cls_fp[[probabilities]][keep] == "bog", mapped) |>
        mutate(rule = "argmax", threshold = NA_real_),
      score(pb >= cut_am, mapped) |>
        mutate(rule = "area-matched cut", threshold = cut_am),
      map_dfr(SWEEP, ~ score(pb >= .x, mapped) |>
        mutate(rule = "P(bog) cut", threshold = .x))
    ) |>
      mutate(population = population, probabilities = probabilities, .before = 1)
  }) |>
  relocate(rule, threshold, .after = probabilities)

cat("\nA. Present-day bog calls inside the footprint:\n")
commission |>
  filter(rule != "P(bog) cut") |>
  select(population, probabilities, rule, threshold, n_mapped_bog, n_pred_bog,
         recall, precision, pred_to_mapped_ratio, on_otherpeat, on_nonpeat) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)
write_csv(commission, "output/pl2/calibration_footprint.csv", append = FALSE)

## B. Calibrated P(bog) by prior shift ####

# Priors. The bog prior is the mapped-bog share of the whole footprint. The remainder is
# split between otherpeat and nonpeat in the ratio of the cells artype_60 labels as such,
# because the model has no class for the mixed and unlabelled cells; they are real non-bog
# cells, so leaving them out would inflate the bog prior.
n_fp <- nrow(fp)
n_lab <- table(fp$mapped)
prior_bog <- n_lab[["bog"]] / n_fp
prior <- c(
  nonpeat = (1 - prior_bog) * n_lab[["nonpeat"]] / (n_lab[["nonpeat"]] + n_lab[["otherpeat"]]),
  otherpeat = (1 - prior_bog) * n_lab[["otherpeat"]] / (n_lab[["nonpeat"]] + n_lab[["otherpeat"]]),
  bog = prior_bog
)[RESPONSE_LEVELS]
cat("\nB. Footprint priors:\n")
print(round(prior, 5))

# p_cal(c) proportional to p_bal(c) x prior(c) / (1/3); the 1/3 cancels on renormalising.
calibrate <- function(M) {
  W <- sweep(M, 2, prior, "*")
  W / rowSums(W)
}
pcal_fp <- calibrate(P_cur_oob[fp$row, , drop = FALSE])[, "bog"]

qs <- c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1)
dist_summary <- function(x, label, scale) {
  tibble(
    set = label,
    scale = scale,
    n = length(x),
    mean = mean(x),
    !!!set_names(as.list(quantile(x, qs, names = FALSE)), paste0("q", qs * 100)),
    share_above_0.5 = mean(x > 0.5)
  )
}
at_mapped <- fp$mapped == "bog"
at_called <- cls_fp$oob == "bog"
prior_shift <- bind_rows(
  dist_summary(pbog_fp$oob[at_mapped], "mapped bog", "balanced"),
  dist_summary(pcal_fp[at_mapped], "mapped bog", "calibrated"),
  dist_summary(pbog_fp$oob[at_called], "argmax bog", "balanced"),
  dist_summary(pcal_fp[at_called], "argmax bog", "calibrated"),
  dist_summary(pcal_fp, "footprint", "calibrated")
) |>
  mutate(
    prior_nonpeat = prior[["nonpeat"]],
    prior_otherpeat = prior[["otherpeat"]],
    prior_bog = prior[["bog"]]
  )
cat("\nP(bog) at present-day mapped and argmax-bog cells in the footprint (OOB):\n")
print(as.data.frame(select(prior_shift, set:share_above_0.5)), row.names = FALSE, digits = 3)
write_csv(prior_shift, "output/pl2/calibration_prior_shift.csv", append = FALSE)

# Reliability diagram. Binned on balanced P(bog); within each bin, mean calibrated P(bog)
# against the observed mapped-bog share. The prior is taken from the same footprint, so
# the overall level agrees by construction; what the diagram tests is whether the
# prior-shift FORM holds across the range of balanced P.
wilson <- function(k, n, side, z = 1.96) {
  p <- k / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  centre + side * half
}
reliability <- tibble(
  pbal = pbog_fp$oob,
  pcal = pcal_fp,
  bog = at_mapped
) |>
  mutate(bin = cut(pbal, REL_BREAKS, include.lowest = TRUE)) |>
  group_by(bin) |>
  summarise(
    n = n(),
    n_bog = sum(bog),
    mean_pbal = mean(pbal),
    mean_pcal = mean(pcal),
    observed = n_bog / n,
    observed_lo = wilson(n_bog, n, -1),
    observed_hi = wilson(n_bog, n, 1),
    .groups = "drop"
  )
cat("\nReliability by balanced P(bog) bin:\n")
print(as.data.frame(reliability), row.names = FALSE, digits = 3)
write_csv(reliability, "output/pl2/calibration_reliability_bins.csv", append = FALSE)

floor_p <- 1e-6
p_rel <- reliability |>
  filter(n_bog > 0) |>
  ggplot(aes(mean_pcal, observed)) +
  geom_abline(linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = pmax(observed_lo, floor_p), ymax = observed_hi), width = 0) +
  geom_point(aes(size = n)) +
  geom_text(aes(label = bin), vjust = -1, size = 2.8) +
  scale_x_log10() +
  scale_y_log10() +
  scale_size_area(max_size = 5) +
  labs(
    x = "Mean calibrated P(bog) in bin",
    y = "Observed mapped-bog share (95% Wilson)",
    size = "Cells",
    title = "Prior-shifted P(bog) against mapped bog, Lyngstad footprint",
    subtitle = "Bins of balanced P(bog), OOB at training cells; bins with no mapped bog omitted"
  ) +
  theme_bw()
ggsave("output/pl2/calibration_reliability.png", p_rel, width = 7, height = 6, dpi = 150)

## C. Loss criteria at present-day bog ####

# Reference: balanced P(bog), OOB, at present-day mapped bog cells in the footprint.
ref_pbal <- pbog_fp$oob[at_mapped]
ref_pcal <- pcal_fp[at_mapped]
ref <- c(
  pbal_p = unname(quantile(ref_pbal, REF_PERCENTILE)),
  pbal_min = min(ref_pbal),
  pcal_median = median(ref_pcal)
)
odds <- function(p) p / (1 - p)
cat("\nC. Reference at mapped bog (OOB):\n")
print(signif(ref, 4))

# The same reference from the production forest's own OOB (~15 trees per bog row), to
# show how much of the 2026-10-07 value was Monte-Carlo noise. Built like P_cur_oob: OOB
# at training cells, projection at the mapped-bog cells that are not training rows.
P_cur_oob_prod <- P[["current"]]
P_cur_oob_prod[train_row, ] <- oob_prod[no_rows, ]
ref_pbal_prod <- P_cur_oob_prod[fp$row[at_mapped], "bog"]
rm(P_cur_oob_prod)
cat(
  "Mapped bog cells in the footprint that are training rows:",
  sum(fp$cell[at_mapped] %in% train_cell), "of", sum(at_mapped), "\n"
)
ref_q <- c(0, 0.05, REF_PERCENTILE, 0.25, 0.5)
bog_oob_rate <- oob_trees$oob_rate_mean[oob_trees$response == "bog"]
ref_row <- function(source, ntree, pbog) {
  tibble(
    source = source,
    ntree = ntree,
    n_mapped_bog = length(pbog),
    oob_trees_per_bog_row = bog_oob_rate * ntree,
    !!!set_names(as.list(quantile(pbog, ref_q, names = FALSE)), paste0("q", ref_q * 100)),
    share_at_1 = mean(pbog == 1)
  )
}
reference_tbl <- bind_rows(
  ref_row("reference fit", REF_NTREE, ref_pbal),
  ref_row("production fit", 1000, ref_pbal_prod)
)
cat("\nOOB balanced P(bog) at mapped bogs, reference vs production fit:\n")
print(as.data.frame(reference_tbl), row.names = FALSE, digits = 3)
write_csv(oob_trees, "output/pl2/calibration_oob_trees.csv", append = FALSE)
write_csv(reference_tbl, "output/pl2/calibration_reference.csv", append = FALSE)

# Criteria, each a function of a probability matrix (rows = cells).
criteria <- list(
  "class flip (argmax bog)" = function(M) argmax(M) == "bog",
  "balanced P(bog) >= reference P10" = function(M) M[, "bog"] >= ref[["pbal_p"]],
  "balanced P(bog) >= reference minimum" = function(M) M[, "bog"] >= ref[["pbal_min"]],
  "calibrated odds >= half reference median" = function(M) {
    odds(calibrate(M)[, "bog"]) >= ODDS_RETAINED * odds(ref[["pcal_median"]])
  }
)

# Two populations of present-day bog:
#  - mapped bog in the footprint: current values are OOB (the reference cells themselves);
#  - argmax bog over Norway today: defined on the projection, as in
#    pl2_compareScenarios.R, so the class-flip rows reproduce its persistence table.
pop_rows <- list(
  "mapped bog (footprint)" = fp$row[at_mapped],
  "argmax bog today (Norway)" = which(argmax(P[["current"]]) == "bog")
)
values_for <- function(pop, s, rows) {
  if (s == "current" && pop == "mapped bog (footprint)") {
    P_cur_oob[rows, , drop = FALSE]
  } else {
    P[[s]][rows, , drop = FALSE]
  }
}

flags <- list()
criteria_tbl <- expand_grid(population = names(pop_rows), scenario = SCEN$scenario) |>
  pmap_dfr(function(population, scenario) {
    rows <- pop_rows[[population]]
    M <- values_for(population, scenario, rows)
    f <- map(criteria, ~ .x(M))
    flags[[paste(population, scenario, sep = "|")]] <<- f
    tibble(
      population = population,
      scenario = scenario,
      criterion = names(criteria),
      threshold = c(NA, ref[["pbal_p"]], ref[["pbal_min"]],
                    ODDS_RETAINED * odds(ref[["pcal_median"]])),
      n_population = length(rows),
      cells_suitable = map_int(f, sum),
      km2_suitable = cells_suitable * cell_km2,
      pct_suitable = 100 * cells_suitable / length(rows)
    )
  })

# Expected bog area over all of Norway: sum of P(bog) x cell area, balanced and calibrated.
# The balanced sum rises under both scenarios (notebook 2026-10-07) and is not a prevalence;
# the calibrated sum carries the footprint prior to the whole domain.
expected_tbl <- map_dfr(SCEN$scenario, function(s) {
  M <- if (s == "current") P_cur_oob else P[[s]]
  tibble(
    population = "Norway",
    scenario = s,
    criterion = c("expected area, balanced P(bog)", "expected area, calibrated P(bog)"),
    n_population = nrow(M),
    km2_suitable = c(sum(M[, "bog"]), sum(calibrate(M)[, "bog"])) * cell_km2
  )
})

criteria_out <- bind_rows(criteria_tbl, expected_tbl)
cat("\nPresent-day bog that stays suitable, by criterion:\n")
criteria_out |>
  select(population, scenario, criterion, cells_suitable, km2_suitable, pct_suitable) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 4)
write_csv(criteria_out, "output/pl2/scenario_comparison_criteria.csv", append = FALSE)

# Consistency with pl2_compareScenarios.R: the class-flip rows for the argmax population
# must reproduce its persistence table.
persist_old <- read_csv(
  "output/pl2/scenario_comparison_persistence_cells.csv", show_col_types = FALSE
)
for (s in c(MILD_SCENARIO, FUTURE_SCENARIO)) {
  old <- persist_old |>
    filter(persist %in% c("bog under both", paste("bog under", s, "only"))) |>
    pull(cells) |>
    sum()
  new <- criteria_tbl |>
    filter(population == "argmax bog today (Norway)", scenario == s,
           criterion == "class flip (argmax bog)") |>
    pull(cells_suitable)
  cat("Class-flip persistence,", s, ": here", new, "| pl2_compareScenarios.R", old, "\n")
}

# Where the criteria disagree: present-day bog cells by which of the class-flip, P10 and
# calibrated criteria still call them suitable under each future scenario.
short <- c(
  "class flip (argmax bog)" = "flip",
  "balanced P(bog) >= reference P10" = "p10",
  "calibrated odds >= half reference median" = "cal"
)
agreement <- expand_grid(
  population = names(pop_rows), scenario = c(MILD_SCENARIO, FUTURE_SCENARIO)
) |>
  pmap_dfr(function(population, scenario) {
    f <- flags[[paste(population, scenario, sep = "|")]][names(short)]
    as_tibble(set_names(f, short)) |>
      count(flip, p10, cal, name = "cells") |>
      mutate(population = population, scenario = scenario, km2 = cells * cell_km2,
             .before = 1)
  })
cat("\nAgreement between criteria (TRUE = still suitable):\n")
print(as.data.frame(agreement), row.names = FALSE, digits = 4)
write_csv(agreement, "output/pl2/scenario_comparison_criteria_agreement.csv", append = FALSE)

# Persistence against the cut. Under SSP1-2.6 the share of present bog kept falls steeply
# between cuts of 0.4 and 0.5 (notebook 2026-11-19), so the strict bound is quoted with
# the sweep, not on the reference P10 alone.
sweep_tbl <- expand_grid(population = names(pop_rows), scenario = SCEN$scenario) |>
  pmap_dfr(function(population, scenario) {
    rows <- pop_rows[[population]]
    pb <- values_for(population, scenario, rows)[, "bog"]
    tibble(
      population = population,
      scenario = scenario,
      cut = c(PERSIST_CUTS, ref[["pbal_p"]]),
      is_reference = c(rep(FALSE, length(PERSIST_CUTS)), TRUE),
      n_population = length(rows),
      cells_suitable = map_int(cut, ~ sum(pb >= .x)),
      km2_suitable = cells_suitable * cell_km2,
      pct_suitable = 100 * cells_suitable / length(rows),
      max_pbog = max(pb)
    )
  })
cat("\nShare of present-day bog with balanced P(bog) >= cut:\n")
sweep_tbl |>
  select(population, scenario, cut, pct_suitable) |>
  mutate(cut = round(cut, 3)) |>
  pivot_wider(names_from = cut, values_from = pct_suitable) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)
write_csv(sweep_tbl, "output/pl2/scenario_comparison_criteria_sweep.csv", append = FALSE)

# sessionInfo ####

sessioninfo::session_info()
