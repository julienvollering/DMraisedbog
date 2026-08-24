# Weighted-PCA exploration of training vs prediction domains ####

# PURPOSE: Diagnostic: weighted-PCA view of training versus prediction domains, to see the extrapolation geometry directly.

# Visualizes a feature-importance-weighted PCA of the predictors, stratified by training
# vs prediction domain, to size the modelling challenge (how far the future prediction
# domains sit from the training data). Companion to pl2_exploreFeatureSpaceDistances.R,
# which measures the same thing one predictor at a time.
#
# Moved out of pl3 when that stage was collapsed into pl2. The EU-anchor layer below is
# vestigial: EU presences are pooled into TRAINING now rather than held out as high-DI
# anchors, so its file.exists() guard finds nothing and the layer is simply not drawn.
#
# The weighted feature space is built to match the CAST dissimilarity index used
# elsewhere: each predictor is standardized on the TRAINING data (current
# scenario) and multiplied by its BRF variable-importance weight. Euclidean
# distance in this space is therefore the weighted distance that DI is based on,
# so the PCA is a 2-D projection of the same space DI/LPD operate in.

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

## Configuration ####

# Subsample size per domain for plotting (avoids overplotting large domains)
plot_sample_size <- 5000

# Number of top loadings (by |contribution| on PC1/PC2) to annotate
n_loadings <- 10

seed <- 42
set.seed(seed)

dir.create("output/pl2", showWarnings = FALSE, recursive = TRUE)

## Read shared upstream data (modeling frame and weights) ####

mf <- read_csv("output/pl2/modeling_frame_regional.csv", show_col_types = FALSE)

weights_features <- read_csv(
  "output/pl2/weights_feature_data_partitioning.csv",
  show_col_types = FALSE
)

# Use the same weighting method as the DI / W_sample work
weights <- weights_features |>
  filter(method == "Balanced Random Forest") |>
  select(feature, median) |>
  arrange(desc(median)) |>
  tibble::deframe()

feat <- names(weights)

# Guard: every weighted feature must be present in the modeling frame
missing <- setdiff(feat, names(mf))
if (length(missing) > 0) {
  stop(
    "Weighted features missing from modeling frame: ",
    paste(missing, collapse = ", ")
  )
}

## Define domains ####

# Training = current scenario (the data the model is fit on: presence + absence)
train <- mf |>
  filter(scenario == "current") |>
  select(response, all_of(feat)) |>
  drop_na()

# Prediction domain 1: full mapped spatial domain under future climate
future_all <- mf |>
  filter(scenario == "future") |>
  select(all_of(feat)) |>
  drop_na()

# Prediction domain 2: future climate at current presence locations
# (the persistence research question's domain)
future_rows <- mf |>
  filter(scenario == "future") |>
  select(x, y, all_of(feat))
future_presence <- mf |>
  filter(scenario == "current", response == "bog") |>
  select(x, y) |>
  left_join(future_rows, by = c("x", "y")) |>
  select(all_of(feat)) |>
  drop_na()

# EU presences (current climate), reconstructed onto the 44 predictors in
# pl3_extractEUpredictors.R. These are leakage-free high-DI ANCHORS, not a
# prediction domain and not in the training/DI reference set; shown here only to
# see where the EU support falls relative to the Norwegian training envelope.
# Optional: only included if that extraction step has already been run.
eu_path <- "output/pl2/eu_presence_predictors.csv"
have_eu <- file.exists(eu_path)
if (have_eu) {
  eu_presence <- read_csv(eu_path, show_col_types = FALSE) |>
    select(all_of(feat)) |>
    drop_na()
}

cat(
  "Domain sizes -- training:",
  nrow(train),
  " future (all):",
  nrow(future_all),
  " future @ presence:",
  nrow(future_presence),
  if (have_eu) paste(" EU presences:", nrow(eu_presence)) else "",
  "\n"
)

## Build the weighted feature-space transform (fit on training only) ####

train_mat <- as.matrix(train[, feat])
ctr <- colMeans(train_mat)
scl <- apply(train_mat, 2, sd)
scl[scl == 0] <- 1 # guard against constant features

w <- weights[feat] # align weight order to feature columns

# z = (x - train_mean) / train_sd, then scale each column by its weight
weighted_transform <- function(m) {
  z <- sweep(m, 2, ctr, "-")
  z <- sweep(z, 2, scl, "/")
  sweep(z, 2, w, "*")
}

## Fit weighted PCA on training, project all domains ####

X_train <- weighted_transform(train_mat)
pca <- prcomp(X_train, center = FALSE, scale. = FALSE)

var_expl <- pca$sdev^2 / sum(pca$sdev^2)
cat(
  "Variance explained -- PC1:",
  round(100 * var_expl[1], 1),
  "%",
  " PC2:",
  round(100 * var_expl[2], 1),
  "%\n"
)

project_pc <- function(df) {
  scores <- weighted_transform(as.matrix(df[, feat])) %*% pca$rotation[, 1:2]
  as.data.frame(scores) |> setNames(c("PC1", "PC2"))
}

scores_train <- project_pc(train) |> mutate(domain = "training")
scores_future_all <- project_pc(future_all) |> mutate(domain = "future (all)")
scores_future_pres <- project_pc(future_presence) |>
  mutate(domain = "future @ presence")
if (have_eu) {
  scores_eu <- project_pc(eu_presence) |> mutate(domain = "EU presences")
}

## Quantify extrapolation magnitude (fraction beyond training PC range) ####

pc1_rng <- range(scores_train$PC1)
pc2_rng <- range(scores_train$PC2)

outside_box <- function(scores) {
  out <- scores$PC1 < pc1_rng[1] |
    scores$PC1 > pc1_rng[2] |
    scores$PC2 < pc2_rng[1] |
    scores$PC2 > pc2_rng[2]
  mean(out)
}

extrap_summary <- tibble(
  domain = c("future (all)", "future @ presence"),
  n = c(nrow(scores_future_all), nrow(scores_future_pres)),
  frac_outside_train_PC12_box = c(
    outside_box(scores_future_all),
    outside_box(scores_future_pres)
  )
)

# EU presences are anchors, not a prediction domain, but reporting how far they
# fall outside the training envelope shows whether they reach the high-DI band.
if (have_eu) {
  extrap_summary <- bind_rows(
    extrap_summary,
    tibble(
      domain = "EU presences",
      n = nrow(scores_eu),
      frac_outside_train_PC12_box = outside_box(scores_eu)
    )
  )
}

cat("\nExtrapolation summary (fraction outside training PC1/PC2 range):\n")
print(as.data.frame(extrap_summary), digits = 3)

write_csv(extrap_summary, "output/pl2/weighted_pca_extrapolation_summary.csv")

## Loadings table ####

loadings <- as.data.frame(pca$rotation[, 1:2]) |>
  tibble::rownames_to_column("feature") |>
  mutate(magnitude = sqrt(PC1^2 + PC2^2)) |>
  arrange(desc(magnitude))

write_csv(loadings, "output/pl2/weighted_pca_loadings.csv")

cat("\nTop loadings on PC1/PC2:\n")
print(as.data.frame(head(loadings, n_loadings)), digits = 3)

## Plot: domains in weighted PC space ####

# Subsample each domain for legible plotting
sample_domain <- function(scores, n) {
  if (nrow(scores) > n) scores[sample.int(nrow(scores), n), ] else scores
}

domain_scores <- list(
  sample_domain(scores_train, plot_sample_size),
  sample_domain(scores_future_all, plot_sample_size),
  sample_domain(scores_future_pres, plot_sample_size)
)
# EU presences are a small anchor set (~100s); plot all of them, unsampled
if (have_eu) {
  domain_scores <- c(domain_scores, list(scores_eu))
}

plot_df <- bind_rows(domain_scores) |>
  mutate(
    domain = factor(
      domain,
      levels = c(
        "training",
        "future (all)",
        "future @ presence",
        "EU presences"
      )
    )
  )

# Training convex hull, to show the envelope the prediction domains must fall in
hull <- scores_train[chull(scores_train$PC1, scores_train$PC2), ]

xlab <- paste0("PC1 (", round(100 * var_expl[1], 1), "%)")
ylab <- paste0("PC2 (", round(100 * var_expl[2], 1), "%)")

g <- ggplot(plot_df, aes(PC1, PC2, color = domain)) +
  geom_polygon(
    data = hull,
    aes(PC1, PC2),
    inherit.aes = FALSE,
    fill = NA,
    color = "grey30",
    linetype = "dashed"
  ) +
  # Dense domains as faint clouds; EU anchors (few points) emphasized on top
  geom_point(
    data = function(d) dplyr::filter(d, domain != "EU presences"),
    alpha = 0.25,
    size = 0.5
  ) +
  geom_point(
    data = function(d) dplyr::filter(d, domain == "EU presences"),
    alpha = 0.9,
    size = 1.6
  ) +
  scale_color_manual(
    values = c(
      "training" = "grey50",
      "future (all)" = "#1b9e77",
      "future @ presence" = "#d95f02",
      "EU presences" = "#7570b3"
    )
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(
    title = "Weighted-PCA: training, future domains, and EU anchors",
    subtitle = "Dashed hull = training envelope; points beyond it are extrapolation",
    x = xlab,
    y = ylab,
    color = "Domain"
  ) +
  theme_minimal()

g

ggsave(
  "output/pl2/weighted_pca_domains.png",
  g,
  width = 8,
  height = 6,
  dpi = 130
)
cat("\nPlot written to output/pl2/weighted_pca_domains.png\n")

# sessionInfo ####

sessioninfo::session_info()
