library(readr)
library(dplyr)
library(terra)
library(cluster)
library(tidyr)
library(ggplot2)
library(sf)
library(purrr)
library(twosamples)

## Functions ####

# Land-use and water screen ####

# ONE rule, stated once, applied to BOTH blocks and to ALL THREE classes -- the same
# principle draw_stratified_absences() exists for, extended from the absence draw to
# absence *eligibility*. Built by pl0_buildLandUseScreen.R from ESA WorldCover 10 m.
#
# WHY IT EXISTS. The first full projection (notebook 2026-08-19) put all 564 Lyngstad
# raised bogs into `nonpeat`, and the cells anchoring that answer were coastal urban
# fringe: 48.5% water and 23.8% built-up or cultivated on average, against 5.1% and 7.3%
# for ordinary Norwegian non-peat. `nonpeat` was conflating three different things --
# land whose climate cannot support peat, land that is half sea, and land that was drained
# and built on. Only the first is evidence about climate, and the error ran in the unsafe
# direction: it charged anthropogenic conversion and open ocean to warming.
#
# WHY BOTH BLOCKS. Screening Norway alone would clean the cold block and leave the warm
# one dirty, and the EU block is what supplies the warm flank -- so the asymmetry would
# land exactly where the projection reads its answer. This is not hypothetical: EU
# non-peat carries MORE human cover than Norway's (16.2% vs 3.1% mean), while Norway
# carries more water (5.6% vs 2.4%). The blocks fail differently, and a screen built from
# only one block's symptom would have missed the other's.
#
# WHY ALL THREE CLASSES. Bog presences carry more human cover than random non-peat
# absences (0.092 vs 0.073 by AR50), because bogs sit in flat lowlands near settlement.
# Applying the screen to absences alone would therefore apply it to the LESS affected
# population and quietly reshape the presence set. At the 50% threshold it removes 4.8% of
# EU bog and 1.8% of Norwegian bog: small, but not nothing, and not silent.
LANDUSE_SCREEN_PATH <- "output/pl0/landuse_screen_250m.tif"

# Chosen 2026-08-19, one threshold for both quantities for simplicity. The sweep behind
# it is recorded in notebook.md (removal at 25/33/50/67%) and is NOT rerun as a
# sensitivity arm -- the numbers stand as reported.
LANDUSE_MAX_HUMAN_PCT <- 50
LANDUSE_MAX_NONLAND_PCT <- 50

# Returns a logical vector, TRUE = keep. `human` is cropland + built-up; `nonland` is
# mapped permanent water plus WorldCover no-data, which is open ocean (the product's
# footprint stops offshore). The two are summed here but kept as separate bands on disk,
# because they are different measurements even though both mean "not usable land".
#
# Cells with no screen value are DROPPED, not kept: absence of WorldCover means absence of
# land, and rule 9's principle applies -- uncovered cells are excluded, never defaulted.
landuse_screen_keep <- function(
  xy,
  label = "",
  screen_path = LANDUSE_SCREEN_PATH,
  max_human = LANDUSE_MAX_HUMAN_PCT,
  max_nonland = LANDUSE_MAX_NONLAND_PCT,
  verbose = TRUE
) {
  scr <- terra::rast(screen_path)
  names(scr) <- c("cropland", "builtup", "water", "nodata")

  v <- terra::extract(scr, as.matrix(xy[, c("x", "y")]))
  human <- v$cropland + v$builtup
  nonland <- v$water + v$nodata

  keep <- !is.na(human) & !is.na(nonland) &
    human <= max_human & nonland <= max_nonland

  if (verbose) {
    cat(sprintf(
      "landuse screen%s: keep %d of %d (%.1f%%) | dropped: human>%d%% %d, nonland>%d%% %d, no data %d\n",
      if (nzchar(label)) paste0(" [", label, "]") else "",
      sum(keep), length(keep), 100 * mean(keep),
      max_human, sum(!is.na(human) & human > max_human),
      max_nonland, sum(!is.na(nonland) & nonland > max_nonland),
      sum(is.na(human) | is.na(nonland))
    ))
  }

  keep
}

# Identify predictors that actually move under the scenario ####

# Returns the predictors whose values differ between the current and future frames at the
# same cell. Everything else is static by construction -- terrain, chiefly -- and a static
# axis contributes exactly zero to any current-to-future distance no matter how much
# weight it carries.
#
# This matters because the DI ruler is VI-weighted and the pooled 3-class fit puts ~44% of
# its weight mass on slope + elevation. That weight is real (slope separates other-peat
# from the rest almost perfectly) but it is inert for a projection: including it dilutes
# the metric and compresses future novelty toward zero, which biases the reliability
# reading in the direction that makes the projection look safer than it is. The predictors
# stay in the MODEL; they are excluded only from the distance metric.
#
# Derived from the data rather than hard-coded, so it stays correct if the predictor set
# changes. Classification is on the FRACTION of matched cells that change, not the maximum
# shift: exactly one cell (x = 4462125, y = 4161875) belongs to both the EU block and the
# Norwegian projection domain -- a single-cell overlap in the domain clip -- and its two
# terrain sources disagree by 0.08 degrees of slope, which is enough to make a max-based
# test call slope "dynamic" off one row in 84,707. The fractions separate cleanly in
# practice (about 1e-5 for static predictors against 1.0 for dynamic ones), so anything
# landing in between is reported rather than silently binned.
# Two distinct restrictions come out of this, and they are used for different things:
#
#  - `dynamic`: the predictor differs at all between scenarios. This is the right filter
#    for the DISTANCE METRIC, because a static axis contributes exactly zero to a
#    current-to-future distance while still inflating the metric's denominator.
#  - `material`: the predictor's mean projected shift is at least `material_sd` standard
#    deviations of its own training spread. This is the right filter for the PARTITION
#    AXIS, because "changes between scenarios" is not the same as "moves far". bio02
#    changes in essentially every cell yet shifts only 0.26 SD, so folds cut on it would
#    separate the data along a gradient the projection barely traverses -- the same
#    failure as cutting on static terrain, just less extreme.
identify_dynamic_predictors <- function(
  mf,
  predictor_names,
  tol = 0,
  id_cols = c("x", "y"),
  dynamic_frac = 0.5,
  material_sd = 0.5
) {
  current <- mf |> dplyr::filter(scenario == "current")
  future <- mf |> dplyr::filter(scenario == "future")

  matched <- dplyr::inner_join(
    current[, c(id_cols, predictor_names)],
    future[, c(id_cols, predictor_names)],
    by = id_cols,
    suffix = c("_cur", "_fut")
  )

  if (nrow(matched) == 0) {
    stop(
      "No cells matched between the current and future frames on ",
      paste(id_cols, collapse = "/")
    )
  }

  changed_frac <- vapply(
    predictor_names,
    function(v) {
      d <- abs(matched[[paste0(v, "_fut")]] - matched[[paste0(v, "_cur")]])
      mean(d > tol, na.rm = TRUE)
    },
    numeric(1)
  )

  # Several genuinely dynamic climate predictors are discrete enough that a few percent
  # of cells land on the same value in both scenarios (ngd0, scd, swe), so the band only
  # flags predictors that are neither clearly static nor clearly dynamic.
  # Mean projected shift, expressed in standard deviations of the predictor's own spread
  # across the training frame so the scale is comparable between, say, gdd5 (degree-days)
  # and bio02 (degrees C).
  shift_sd <- vapply(
    predictor_names,
    function(v) {
      shift <- mean(
        matched[[paste0(v, "_fut")]] - matched[[paste0(v, "_cur")]],
        na.rm = TRUE
      )
      spread <- stats::sd(current[[v]], na.rm = TRUE)
      if (is.finite(spread) && spread > 0) shift / spread else 0
    },
    numeric(1)
  )

  ambiguous <- names(changed_frac)[changed_frac > 0.01 & changed_frac < 0.9]
  if (length(ambiguous) > 0) {
    warning(
      "Predictors change in some but not most cells, so static/dynamic is not clear-cut: ",
      paste(ambiguous, collapse = ", ")
    )
  }

  dynamic <- names(changed_frac)[changed_frac >= dynamic_frac]

  list(
    dynamic = dynamic,
    static = names(changed_frac)[changed_frac < dynamic_frac],
    material = dynamic[abs(shift_sd[dynamic]) >= material_sd],
    changed_frac = changed_frac,
    shift_sd = shift_sd,
    n_matched = nrow(matched)
  )
}

# Feature scaling for the weighted distance metric ####

# The centre/spread used to put features on a common footing before the VI weights are
# applied. Deliberately a property of the FEATURE SPACE, not of a fold: derived once from
# the pooled training frame and handed to every call of calculate_weighted_di().
#
# Deriving it per fold instead is what CAST's trainDI does, and it is defensible when
# folds are broad -- but the presence-sorted partitions here are extremely uneven in
# width (partitions 2-4 span 0.5, 0.28 and 0.20 degrees of bio10 against partition 1's
# 5.2), and inside the narrow ones some predictors are near-constant. `gsp` has an
# essentially zero within-fold SD in partition 4, which put partition-1 rows 1.5 million
# fold-SDs away on that one axis and drove DI to 1.9e5. Freezing the scaling removes a
# numerical artefact and makes DI comparable between folds, which the pairwise design
# needs anyway. The normalisation constant stays per-fold: THAT is what makes DI relative
# to a particular training set.
compute_feature_scaling <- function(data, feature_names) {
  m <- as.matrix(data[, feature_names, drop = FALSE])
  center <- colMeans(m, na.rm = TRUE)
  scale <- apply(m, 2, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(center = center, scale = scale)
}

# Calculate weighted Euclidean Dissimilarity Index from test to train ####
# Implements custom DI calculation using weighted Euclidean distance
# Similar to CAST::trainDI but optimized for pairwise train-test comparisons
#
# DI is per row of `test_data` -- one value per cell, never one per fold or per cut.
# Aggregates over it (a fold's mean, a cut's distribution) are summaries computed
# downstream; this function's output is always the full per-row vector.
#
# Two arguments exist for the frozen-ruler requirement (plan_EUintegration.md section
# 1.4): the pl3 sweep must measure every cut against ONE ruler, so `scaling` and
# `train_avg_dist` can be supplied from the production fit instead of being re-derived
# per training set. Left NULL (the pl2 default) they are derived from `train_data`, which
# is the ordinary definition of DI relative to that training set.
calculate_weighted_di <- function(
  train_data,
  test_data,
  weights = NULL,
  scaling = NULL,
  train_avg_dist = NULL,
  norm_sample = 2000,
  seed = 42,
  verbose = TRUE
) {
  # Inputs:
  #   train_data: data.frame/tibble with training features (no response column)
  #   test_data: data.frame/tibble with test features (no response column)
  #   weights: optional named vector of variable importance weights
  #   scaling: optional list(center =, scale =) to freeze the feature scaling
  #   train_avg_dist: optional pre-computed normalisation constant to freeze the ruler
  #   norm_sample: training rows sampled to estimate the normalisation constant
  #   seed: seed for that subsample
  #   verbose: logical, whether to report progress
  #
  # Returns:
  #   List with DI vector for test data, train_avg_dist, and the scaling used

  # Ensure both datasets have same columns in same order
  common_cols <- intersect(colnames(train_data), colnames(test_data))
  train_data <- as.matrix(train_data[, common_cols, drop = FALSE])
  test_data <- as.matrix(test_data[, common_cols, drop = FALSE])

  # Scale features (training data center/scale for both, unless frozen scaling is given)
  if (is.null(scaling)) {
    train_scaled <- scale(train_data)
    center_vec <- attr(train_scaled, "scaled:center")
    scale_vec <- attr(train_scaled, "scaled:scale")
  } else {
    center_vec <- scaling$center[common_cols]
    scale_vec <- scaling$scale[common_cols]
    train_scaled <- scale(train_data, center = center_vec, scale = scale_vec)
  }

  # A constant feature has scale 0 and would propagate NaN through every distance. A
  # NEAR-constant one is worse than that, because it silently survives: dividing by a
  # near-zero spread turns an ordinary difference into an enormous one, and a single such
  # feature can dominate the whole metric. Warn rather than repair, since the fix belongs
  # upstream (freeze the scaling with compute_feature_scaling()).
  degenerate <- !is.finite(scale_vec) | scale_vec == 0
  scale_vec[degenerate] <- 1
  if (any(degenerate)) {
    warning(
      "Zero or non-finite spread for: ",
      paste(names(scale_vec)[degenerate], collapse = ", "),
      " -- set to 1"
    )
  }
  train_scaled <- scale(train_data, center = center_vec, scale = scale_vec)
  test_scaled <- scale(test_data, center = center_vec, scale = scale_vec)

  # Apply variable importance weights if provided
  if (!is.null(weights)) {
    # Match weight order to column order
    weight_vec <- weights[common_cols]
    # Set negative weights to 0
    weight_vec[weight_vec < 0] <- 0
    # Apply weights by multiplying each column by weight (following trainDI approach)
    train_scaled <- sweep(train_scaled, 2, weight_vec, "*")
    test_scaled <- sweep(test_scaled, 2, weight_vec, "*")
  }

  # Average training-to-training distance (the normalisation constant) ####
  # Estimated from a random subsample of training rows against ALL training rows. The
  # constant is the grand mean of the pairwise distances, so a subsample of rows gives an
  # unbiased estimate of it; computing all n^2 pairs costs hours at pooled-frame partition
  # sizes and buys a decimal place that no downstream comparison can see.
  if (is.null(train_avg_dist)) {
    n_train <- nrow(train_scaled)
    m <- min(norm_sample, n_train)

    if (!is.null(seed)) {
      set.seed(seed)
    }
    anchor_idx <- if (m < n_train) sort(sample.int(n_train, m)) else seq_len(n_train)

    if (verbose) {
      message(
        "Computing DI normalisation constant from ", m, " of ", n_train,
        " training rows..."
      )
    }

    # Distances via ||a-b||^2 = ||a||^2 + ||b||^2 - 2a.b, so the inner loop is one BLAS
    # matrix product rather than m R-level sweeps over the whole training matrix.
    # Anchors are processed in blocks to bound the intermediate at block x n_train.
    train_sq <- rowSums(train_scaled^2)
    dist_total <- 0
    dist_count <- 0
    block <- 250

    for (start in seq(1, m, by = block)) {
      idx <- anchor_idx[start:min(start + block - 1, m)]
      anchors <- train_scaled[idx, , drop = FALSE]

      d2 <- outer(rowSums(anchors^2), train_sq, "+") -
        2 * tcrossprod(anchors, train_scaled)
      d2[d2 < 0] <- 0 # Guard against negative values from rounding
      d <- sqrt(d2)

      # Exclude each anchor's distance to itself
      d[cbind(seq_along(idx), idx)] <- NA

      dist_total <- dist_total + sum(d, na.rm = TRUE)
      dist_count <- dist_count + sum(!is.na(d))
    }

    # Every anchor contributes the same number of distances (n_train - 1), so the grand
    # mean over all pairs equals the mean of the per-anchor means the loop used to form.
    train_avg_dist <- dist_total / dist_count
  }

  # Minimum test-to-train distance for each test point ####
  # Exact nearest neighbour via a kd-tree; the previous row-by-row scan was the same
  # quantity computed in O(n_test x n_train) R-level operations.
  if (verbose) {
    message("Computing DI of ", nrow(test_scaled), " test rows...")
  }

  nn <- FNN::get.knnx(
    data = train_scaled,
    query = test_scaled,
    k = 1,
    algorithm = "kd_tree"
  )
  test_min_dist <- nn$nn.dist[, 1]

  # Dissimilarity Index = normalized minimum distance
  DI <- test_min_dist / train_avg_dist

  # Return results
  return(list(
    DI = DI,
    train_avg_dist = train_avg_dist,
    scaling = list(center = center_vec, scale = scale_vec)
  ))
}

# Quantile cut-based partitioning with probabilistic interleaving/segregation
create_cut_based_partitions <- function(
  data,
  feature_name,
  target_col,
  n_cuts,
  n_partitions = 5,
  segregation_prob = 0,
  seed = 42
) {
  cat(
    "Creating",
    n_partitions,
    "partitions using",
    n_cuts,
    "quantile cuts on feature '",
    feature_name,
    "' (segregation_prob =",
    segregation_prob,
    ")...\n"
  )

  # Extract feature values
  feature_values <- data[[feature_name]]

  # Create quantile-based cuts
  quantile_breaks <- quantile(
    feature_values,
    probs = seq(0, 1, length.out = n_cuts + 1),
    na.rm = TRUE
  )

  # Handle case where quantiles are not unique (can happen with discrete data)
  quantile_breaks <- unique(quantile_breaks)
  actual_cuts <- length(quantile_breaks) - 1

  if (actual_cuts < n_cuts) {
    cat(
      "Warning: Only",
      actual_cuts,
      "unique quantile breaks found (requested",
      n_cuts,
      ")\n"
    )
  }

  intervals <- cut(
    feature_values,
    breaks = quantile_breaks,
    labels = FALSE,
    include.lowest = TRUE
  )

  # Probabilistic interleaving/segregation assignment
  # segregation_prob = 0: Pure interleaving (all intervals distributed round-robin)
  # segregation_prob = 1: Maximum segregation (intervals form contiguous blocks)

  # Guarantee all n_partitions are non-empty by assigning first n_partitions intervals
  # one-per-partition, then apply probabilistic logic to remaining intervals

  set.seed(seed)
  unique_intervals_ordered <- sort(unique(intervals))
  n_intervals <- length(unique_intervals_ordered)
  interval_to_partition <- integer(n_intervals)

  # Step 1: Assign first n_partitions intervals (one per partition) to guarantee coverage
  if (n_intervals >= n_partitions) {
    interval_to_partition[1:n_partitions] <- 1:n_partitions

    # Step 2: Apply probabilistic logic to remaining intervals
    if (n_intervals > n_partitions) {
      for (i in (n_partitions + 1):n_intervals) {
        if (runif(1) < segregation_prob) {
          # SEGREGATE: Join contiguous block (inherit previous partition)
          interval_to_partition[i] <- interval_to_partition[i - 1]
        } else {
          # INTERLEAVE: Round-robin assignment
          interval_to_partition[i] <- ((i - 1) %% n_partitions) + 1
        }
      }
    }
  } else {
    # Edge case: fewer intervals than partitions
    # Assign round-robin (some partitions will be empty)
    cat(
      "Warning: Only",
      n_intervals,
      "intervals but",
      n_partitions,
      "partitions requested.\n"
    )
    for (i in seq_along(unique_intervals_ordered)) {
      interval_to_partition[i] <- ((i - 1) %% n_partitions) + 1
    }
  }

  # Map intervals to their assigned partitions
  partition_assignments <- interval_to_partition[intervals]

  # Calculate partition statistics
  partition_stats <- data |>
    mutate(partition = factor(partition_assignments)) |>
    group_by(partition) |>
    summarise(
      n_total = n(),
      n_presence = sum(.data[[target_col]] == 1),
      prevalence = n_presence / n_total,
      .groups = "drop"
    )

  cat("\nPartition summary:\n")
  print(partition_stats)

  cat("\nPartition assignment:\n")
  for (p in 1:n_partitions) {
    cat(sprintf(
      "Partition %d: %d obs (%.1f%%), %.3f prevalence\n",
      p,
      partition_stats$n_total[p],
      partition_stats$n_total[p] / sum(partition_stats$n_total) * 100,
      partition_stats$prevalence[p]
    ))
  }

  return(list(
    partitions = partition_assignments,
    partition_stats = partition_stats
  ))
}

# Spacer-based partitioning ####
create_spacer_based_partitions <- function(
  data,
  feature_name,
  target_col,
  n_cuts,
  spacer_proportion = 0,
  n_partitions = 5,
  seed = 42
) {
  cat(
    "Creating",
    n_partitions,
    "partitions using",
    n_cuts,
    "quantile cuts on feature '",
    feature_name,
    "' with spacer_proportion =",
    spacer_proportion,
    "...\n"
  )

  # Extract feature values
  feature_values <- data[[feature_name]]

  # Create quantile-based cuts (same as cut-based approach)
  quantile_breaks <- quantile(
    feature_values,
    probs = seq(0, 1, length.out = n_cuts + 1),
    na.rm = TRUE
  )

  # Handle case where quantiles are not unique (can happen with discrete data)
  quantile_breaks <- unique(quantile_breaks)
  actual_cuts <- length(quantile_breaks) - 1

  if (actual_cuts < n_cuts) {
    cat(
      "Warning: Only",
      actual_cuts,
      "unique quantile breaks found (requested",
      n_cuts,
      ")\n"
    )
  }

  intervals <- cut(
    feature_values,
    breaks = quantile_breaks,
    labels = FALSE,
    include.lowest = TRUE
  )

  # Validate inputs
  if (spacer_proportion < 0 || spacer_proportion > 1) {
    stop("spacer_proportion must be in [0, 1]")
  }

  # Spacer-based assignment: Create contiguous blocks with spacer gaps
  set.seed(seed)
  unique_intervals_ordered <- sort(unique(intervals))
  n_intervals <- length(unique_intervals_ordered)

  # Check if we have enough intervals
  if (n_intervals < n_partitions) {
    stop(sprintf(
      "Insufficient intervals (%d) for %d partitions",
      n_intervals, n_partitions
    ))
  }

  interval_to_partition <- integer(n_intervals)

  # Calculate spacer and partition allocations
  n_spacers <- n_partitions - 1
  total_intervals <- n_intervals

  spacer_intervals_total <- floor(total_intervals * spacer_proportion)
  spacer_intervals_per_gap <- floor(spacer_intervals_total / n_spacers)
  remaining_spacer_intervals <- spacer_intervals_total - (spacer_intervals_per_gap * n_spacers)

  partition_intervals_total <- total_intervals - spacer_intervals_total
  partition_intervals_per_block <- floor(partition_intervals_total / n_partitions)
  remaining_partition_intervals <- partition_intervals_total - (partition_intervals_per_block * n_partitions)

  # Build partition and spacer size vectors
  partition_sizes <- rep(partition_intervals_per_block, n_partitions)
  if (remaining_partition_intervals > 0) {
    partition_sizes[1:remaining_partition_intervals] <-
      partition_sizes[1:remaining_partition_intervals] + 1
  }

  spacer_sizes <- rep(spacer_intervals_per_gap, n_spacers)
  if (remaining_spacer_intervals > 0) {
    spacer_sizes[1:remaining_spacer_intervals] <-
      spacer_sizes[1:remaining_spacer_intervals] + 1
  }

  # Build interval-to-partition mapping with spacers
  current_interval <- 1

  for (p in 1:n_partitions) {
    # Assign partition block
    end_interval <- current_interval + partition_sizes[p] - 1
    if (end_interval > n_intervals) end_interval <- n_intervals

    interval_to_partition[current_interval:end_interval] <- p
    current_interval <- end_interval + 1

    # Assign spacer gap (if not last partition)
    if (p < n_partitions && current_interval <= n_intervals) {
      end_interval <- current_interval + spacer_sizes[p] - 1
      if (end_interval > n_intervals) end_interval <- n_intervals

      interval_to_partition[current_interval:end_interval] <- NA
      current_interval <- end_interval + 1
    }
  }

  # Map intervals to their assigned partitions
  partition_assignments <- interval_to_partition[intervals]

  # Calculate partition statistics (excluding spacer rows)
  non_spacer_mask <- !is.na(partition_assignments)

  partition_stats <- data[non_spacer_mask, ] |>
    mutate(partition = factor(partition_assignments[non_spacer_mask])) |>
    group_by(partition) |>
    summarise(
      n_total = n(),
      n_presence = sum(.data[[target_col]] == 1),
      prevalence = n_presence / n_total,
      .groups = "drop"
    )

  # Calculate spacer statistics
  n_spacer_rows <- sum(!non_spacer_mask)
  spacer_prevalence <- if (n_spacer_rows > 0) {
    sum(data[[target_col]][!non_spacer_mask] == 1) / n_spacer_rows
  } else {
    0
  }

  # Print diagnostics
  cat("\nSpacer-based partitioning summary:\n")
  cat("  Total intervals:", n_intervals, "\n")
  cat("  Spacer proportion:", round(spacer_proportion, 3), "\n")
  cat("  Spacer intervals:", spacer_intervals_total, "(",
      round(spacer_intervals_total / n_intervals * 100, 1), "%)\n")
  cat("  Partition intervals:", partition_intervals_total, "(",
      round(partition_intervals_total / n_intervals * 100, 1), "%)\n")
  cat("  Spacer rows (data excluded):", n_spacer_rows, "(",
      round(n_spacer_rows / nrow(data) * 100, 1), "%)\n")
  if (n_spacer_rows > 0) {
    cat("  Spacer prevalence:", round(spacer_prevalence, 3), "\n")
  }
  cat("\n")

  cat("Partition summary (excluding spacers):\n")
  print(partition_stats)

  cat("\nPartition assignment:\n")
  for (p in 1:n_partitions) {
    cat(sprintf(
      "Partition %d: %d obs (%.1f%%), %.3f prevalence\n",
      p,
      partition_stats$n_total[p],
      partition_stats$n_total[p] / sum(partition_stats$n_total) * 100,
      partition_stats$prevalence[p]
    ))
  }

  # Validation: Warn if any partition is too small
  min_partition_size <- min(partition_stats$n_total)
  if (min_partition_size < 100) {
    warning(sprintf(
      "Smallest partition has only %d observations (< 100)",
      min_partition_size
    ))
  }

  return(list(
    partitions = partition_assignments,
    partition_stats = partition_stats
  ))
}

# Evaluate spacer-based partitioning for a given n_cuts and spacer_proportion
evaluate_spacer_partitioning <- function(
  n_cuts,
  spacer_proportion = 0,
  data,
  precomputed_distances,
  feature_name,
  featuredist_sample_size = 1e4,
  seed = 42
) {
  cat(
    "Evaluating n_cuts =",
    n_cuts,
    ", spacer_proportion =",
    spacer_proportion,
    "...\n"
  )

  # Create partitions using spacer-based approach
  partition_result <- create_spacer_based_partitions(
    data = data,
    feature_name = feature_name,
    target_col = "response",
    n_cuts = n_cuts,
    n_partitions = 5,
    spacer_proportion = spacer_proportion,
    seed = seed
  )

  # Filter out spacer rows (partition == NA)
  cv_folds <- partition_result$partitions
  non_spacer_mask <- !is.na(cv_folds)

  cat(
    "Non-spacer data:", sum(non_spacer_mask), "rows",
    "(", round(sum(non_spacer_mask) / length(non_spacer_mask) * 100, 1), "%)\n"
  )

  # Calculate CV distances using only non-spacer data
  cv_dists <- calculate_cv_distances(
    training_data = data[non_spacer_mask, ],
    variable_name = feature_name,
    cv_folds = cv_folds[non_spacer_mask],
    sample_size = featuredist_sample_size,
    standardize = FALSE,
    seed = seed
  )

  # Debug: Print CV distance statistics
  cat(
    "CV distances - n:",
    length(cv_dists),
    ", mean:",
    round(mean(cv_dists, na.rm = TRUE), 4),
    ", range:",
    paste(round(range(cv_dists, na.rm = TRUE), 4), collapse = " - "),
    "\n"
  )

  # Combine all distances and create result structure (unscaled)
  all_dists <- c(
    precomputed_distances$sample_to_sample,
    precomputed_distances$prediction_to_sample,
    cv_dists
  )

  distance_types <- factor(
    c(
      rep("sample-to-sample", length(precomputed_distances$sample_to_sample)),
      rep(
        "prediction-to-sample",
        length(precomputed_distances$prediction_to_sample)
      ),
      rep("CV-distances", length(cv_dists))
    ),
    levels = c("sample-to-sample", "prediction-to-sample", "CV-distances")
  )

  # Create feature_dists result structure for compatibility (no scaling)
  feature_dists <- tibble::tibble(
    dist = all_dists,
    what = distance_types,
    dist_type = "feature"
  )

  # Set class and attributes similar to featuredist output
  class(feature_dists) <- c("geodist", class(feature_dists))
  attr(feature_dists, "type") <- "feature"

  # Calculate W_CV statistic
  W_CV <- twosamples::wass_stat(
    feature_dists[feature_dists$what == "CV-distances", "dist"][[1]],
    feature_dists[feature_dists$what == "prediction-to-sample", "dist"][[1]]
  )

  # Print results
  cat(
    "n_cuts =",
    n_cuts,
    ", spacer_proportion =",
    spacer_proportion,
    "-> W_CV =",
    round(W_CV, 4),
    "\n\n"
  )

  return(list(
    n_cuts = n_cuts,
    spacer_proportion = spacer_proportion,
    partition_result = partition_result,
    W_CV = W_CV,
    feature_name = feature_name,
    feature_dists = feature_dists
  ))
}

# Evaluate cut-based partitioning for a given n_cuts and segregation_prob
evaluate_cut_partitioning <- function(
  n_cuts,
  segregation_prob = 0,
  data,
  precomputed_distances,
  feature_name,
  featuredist_sample_size = 1e4,
  seed = 42
) {
  cat(
    "Evaluating n_cuts =",
    n_cuts,
    ", segregation_prob =",
    segregation_prob,
    "...\n"
  )

  # Create partitions using cut-based approach
  partition_result <- create_cut_based_partitions(
    data = data,
    feature_name = feature_name,
    target_col = "response",
    n_cuts = n_cuts,
    n_partitions = 5, # Keep 5 partitions (1 test + 4 CV)
    segregation_prob = segregation_prob,
    seed = seed
  )

  # Use partition assignments directly as CV folds
  cv_folds <- partition_result$partitions

  # Calculate CV distances across all partitions
  cv_dists <- calculate_cv_distances(
    training_data = data,
    variable_name = feature_name,
    cv_folds = cv_folds,
    sample_size = featuredist_sample_size,
    standardize = FALSE,
    seed = seed
  )

  # Debug: Print CV distance statistics
  cat(
    "CV distances - n:",
    length(cv_dists),
    ", mean:",
    round(mean(cv_dists, na.rm = TRUE), 4),
    ", range:",
    round(range(cv_dists, na.rm = TRUE), 4),
    "\n"
  )

  # Combine all distances and create result structure (unscaled)
  all_dists <- c(
    precomputed_distances$sample_to_sample,
    precomputed_distances$prediction_to_sample,
    cv_dists
  )

  distance_types <- factor(
    c(
      rep("sample-to-sample", length(precomputed_distances$sample_to_sample)),
      rep(
        "prediction-to-sample",
        length(precomputed_distances$prediction_to_sample)
      ),
      rep("CV-distances", length(cv_dists))
    ),
    levels = c("sample-to-sample", "prediction-to-sample", "CV-distances")
  )

  # Create feature_dists result structure for compatibility (no scaling)
  feature_dists <- tibble::tibble(
    dist = all_dists,
    what = distance_types,
    dist_type = "feature"
  )

  # Set class and attributes similar to featuredist output
  class(feature_dists) <- c("geodist", class(feature_dists))
  attr(feature_dists, "type") <- "feature"

  # Calculate W_CV statistic
  W_CV <- twosamples::wass_stat(
    feature_dists[feature_dists$what == "CV-distances", "dist"][[1]],
    feature_dists[feature_dists$what == "prediction-to-sample", "dist"][[1]]
  )

  # Print results
  cat(
    "n_cuts =",
    n_cuts,
    ", segregation_prob =",
    segregation_prob,
    "-> W_CV =",
    round(W_CV, 4),
    "\n\n"
  )

  return(list(
    n_cuts = n_cuts,
    segregation_prob = segregation_prob,
    partition_result = partition_result,
    W_CV = W_CV,
    feature_name = feature_name,
    feature_dists = feature_dists
  ))
}

featuredist <- function(
  training_data,
  prediction_data,
  variable_name,
  cvfold_column = NULL,
  sample_size = 1e4,
  seed = NULL,
  standardize = FALSE,
  scale_distances = FALSE
) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Sample training and prediction data
  n_train <- min(nrow(training_data), sample_size)
  n_pred <- min(nrow(prediction_data), sample_size)

  train_sampled <- training_data[sample.int(nrow(training_data), n_train), ]
  pred_sampled <- prediction_data[sample.int(nrow(prediction_data), n_pred), ]

  # Extract variable values
  train_var <- train_sampled[[variable_name]]
  pred_var <- pred_sampled[[variable_name]]

  # Remove NA values
  train_clean <- train_var[!is.na(train_var)]
  pred_clean <- pred_var[!is.na(pred_var)]

  # Standardize if requested
  if (standardize) {
    # Min-max standardization using combined training + prediction data
    combined_values <- c(train_clean, pred_clean)
    min_val <- min(combined_values)
    max_val <- max(combined_values)

    if (max_val > min_val) {
      train_clean <- (train_clean - min_val) / (max_val - min_val)
      pred_clean <- (pred_clean - min_val) / (max_val - min_val)
    } else {
      # All values identical - set to 0.5 (middle of [0,1] range)
      train_clean <- rep(0.5, length(train_clean))
      pred_clean <- rep(0.5, length(pred_clean))
    }
  }

  # Convert to matrices for FNN
  train_matrix <- matrix(train_clean, ncol = 1)
  pred_matrix <- matrix(pred_clean, ncol = 1)

  # Calculate sample-to-sample distances
  s2s_dists <- numeric(length(train_clean))
  for (i in seq_along(train_clean)) {
    dists <- FNN::knnx.dist(
      train_matrix[i, , drop = FALSE],
      train_matrix,
      k = 1
    )
    dists[i] <- NA # Exclude self-distance
    s2s_dists[i] <- min(dists, na.rm = TRUE)
  }

  # Calculate prediction-to-sample distances
  s2p_dists <- numeric(length(pred_clean))
  for (i in seq_along(pred_clean)) {
    dists <- FNN::knnx.dist(pred_matrix[i, , drop = FALSE], train_matrix, k = 1)
    s2p_dists[i] <- min(dists, na.rm = TRUE)
  }

  # Combine results
  result <- tibble::tibble(
    dist = c(s2s_dists, s2p_dists),
    what = factor(
      c(
        rep("sample-to-sample", length(s2s_dists)),
        rep("prediction-to-sample", length(s2p_dists))
      ),
      levels = c("sample-to-sample", "prediction-to-sample", "CV-distances")
    ),
    dist_type = "feature"
  )

  # Add CV distances if fold column is provided
  if (!is.null(cvfold_column)) {
    cv_folds <- train_sampled[[cvfold_column]]
    cv_folds <- cv_folds[!is.na(train_var)] # Match cleaned training data

    cv_dists <- numeric(0)
    unique_folds <- unique(cv_folds)

    for (fold in unique_folds) {
      test_idx <- which(cv_folds == fold)
      train_idx <- which(cv_folds != fold)

      if (length(test_idx) > 0 && length(train_idx) > 0) {
        test_matrix <- matrix(train_clean[test_idx], ncol = 1)
        fold_train_matrix <- matrix(train_clean[train_idx], ncol = 1)

        for (i in seq_along(test_idx)) {
          dists <- FNN::knnx.dist(
            test_matrix[i, , drop = FALSE],
            fold_train_matrix,
            k = 1
          )
          cv_dists <- c(cv_dists, min(dists, na.rm = TRUE))
        }
      }
    }

    # Add CV distances to result
    cv_result <- tibble::tibble(
      dist = cv_dists,
      what = factor(rep("CV-distances", length(cv_dists))),
      dist_type = "feature"
    )

    result <- dplyr::bind_rows(result, cv_result)
  }

  # Scale distances if requested
  if (scale_distances) {
    max_dist <- max(result$dist, na.rm = TRUE)
    if (max_dist > 0) {
      result$dist <- result$dist / max_dist
    }
  }

  # Set class and attributes similar to geodist
  class(result) <- c("geodist", class(result))
  attr(result, "type") <- "feature"

  # Calculate W statistics
  W_sample <- twosamples::wass_stat(
    result[result$what == "sample-to-sample", "dist"][[1]],
    result[result$what == "prediction-to-sample", "dist"][[1]]
  )
  attr(result, "W_sample") <- W_sample

  if (!is.null(cvfold_column)) {
    W_CV <- twosamples::wass_stat(
      result[result$what == "CV-distances", "dist"][[1]],
      result[result$what == "prediction-to-sample", "dist"][[1]]
    )
    attr(result, "W_CV") <- W_CV
  }

  return(result)
}

# Calculate prediction-to-sample distances (independent of k)
calculate_prediction_distances <- function(
  training_data,
  prediction_data,
  variable_name,
  sample_size = 1e4,
  standardize = FALSE,
  seed = NULL
) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Sample training and prediction data
  n_train <- min(nrow(training_data), sample_size)
  n_pred <- min(nrow(prediction_data), sample_size)

  train_sampled <- training_data[sample.int(nrow(training_data), n_train), ]
  pred_sampled <- prediction_data[sample.int(nrow(prediction_data), n_pred), ]

  # Extract variable values
  train_var <- train_sampled[[variable_name]]
  pred_var <- pred_sampled[[variable_name]]

  # Remove NA values
  train_clean <- train_var[!is.na(train_var)]
  pred_clean <- pred_var[!is.na(pred_var)]

  # Standardize if requested
  if (standardize) {
    # Min-max standardization using combined training + prediction data
    combined_values <- c(train_clean, pred_clean)
    min_val <- min(combined_values)
    max_val <- max(combined_values)

    if (max_val > min_val) {
      train_clean <- (train_clean - min_val) / (max_val - min_val)
      pred_clean <- (pred_clean - min_val) / (max_val - min_val)
    } else {
      # All values identical - set to 0.5 (middle of [0,1] range)
      train_clean <- rep(0.5, length(train_clean))
      pred_clean <- rep(0.5, length(pred_clean))
    }
  }

  # Convert to matrices for FNN
  train_matrix <- matrix(train_clean, ncol = 1)
  pred_matrix <- matrix(pred_clean, ncol = 1)

  # Calculate sample-to-sample distances
  s2s_dists <- numeric(length(train_clean))
  for (i in seq_along(train_clean)) {
    dists <- FNN::knnx.dist(
      train_matrix[i, , drop = FALSE],
      train_matrix,
      k = 1
    )
    dists[i] <- NA # Exclude self-distance
    s2s_dists[i] <- min(dists, na.rm = TRUE)
  }

  # Calculate prediction-to-sample distances
  s2p_dists <- numeric(length(pred_clean))
  for (i in seq_along(pred_clean)) {
    dists <- FNN::knnx.dist(pred_matrix[i, , drop = FALSE], train_matrix, k = 1)
    s2p_dists[i] <- min(dists, na.rm = TRUE)
  }

  return(list(
    sample_to_sample = s2s_dists,
    prediction_to_sample = s2p_dists,
    train_matrix = train_matrix,
    train_clean = train_clean
  ))
}

# Calculate CV distances only (depends on k via fold assignments)
calculate_cv_distances <- function(
  training_data,
  variable_name,
  cv_folds,
  sample_size = 1e4,
  standardize = FALSE,
  seed = NULL
) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Sample training data (consistent with precomputation approach)
  n_train <- min(nrow(training_data), sample_size)
  sample_idx <- sample.int(nrow(training_data), n_train)
  train_sampled <- training_data[sample_idx, ]
  cv_folds_sampled <- cv_folds[sample_idx]

  # Extract variable values
  train_var <- train_sampled[[variable_name]]

  # Remove NA values
  train_clean <- train_var[!is.na(train_var)]
  cv_folds_clean <- cv_folds_sampled[!is.na(train_var)]

  # Standardize if requested (consistent with precomputation)
  if (standardize) {
    # Min-max standardization using full training data
    min_val <- min(train_clean)
    max_val <- max(train_clean)

    if (max_val > min_val) {
      train_clean <- (train_clean - min_val) / (max_val - min_val)
    } else {
      # All values identical - set to 0.5 (middle of [0,1] range)
      train_clean <- rep(0.5, length(train_clean))
    }
  }

  cv_dists <- numeric(0)
  unique_folds <- unique(cv_folds_clean)

  # Remove NA folds if any
  unique_folds <- unique_folds[!is.na(unique_folds)]

  for (fold in unique_folds) {
    test_idx <- which(cv_folds_clean == fold)
    train_idx <- which(cv_folds_clean != fold)

    if (length(test_idx) > 0 && length(train_idx) > 0) {
      test_matrix <- matrix(train_clean[test_idx], ncol = 1)
      fold_train_matrix <- matrix(train_clean[train_idx], ncol = 1)

      for (i in seq_along(test_idx)) {
        dists <- FNN::knnx.dist(
          test_matrix[i, , drop = FALSE],
          fold_train_matrix,
          k = 1
        )
        cv_dists <- c(cv_dists, min(dists, na.rm = TRUE))
      }
    }
  }

  return(cv_dists)
}


# Function to calculate G-mean from predictions and true values
# Robust to factor level ordering by explicitly identifying positive/negative classes
calculate_gmean <- function(predicted_class, true_class) {
  # Convert to character to avoid factor level issues
  predicted_class <- as.character(predicted_class)
  true_class <- as.character(true_class)

  # Create confusion matrix
  cm <- table(Actual = true_class, Predicted = predicted_class)

  if (nrow(cm) == 2 && ncol(cm) == 2) {
    # Explicitly identify which row/column corresponds to positive class ("1")
    pos_row <- which(rownames(cm) == "1")
    pos_col <- which(colnames(cm) == "1")
    neg_row <- which(rownames(cm) == "0")
    neg_col <- which(colnames(cm) == "0")

    # Handle case where positive or negative class is missing
    if (
      length(pos_row) == 0 ||
        length(pos_col) == 0 ||
        length(neg_row) == 0 ||
        length(neg_col) == 0
    ) {
      tpr <- 0
      tnr <- 0
      gmean <- 0
    } else {
      # Extract confusion matrix values using explicit indexing
      tp <- cm[pos_row, pos_col] # True positives
      fn <- cm[pos_row, neg_col] # False negatives
      fp <- cm[neg_row, pos_col] # False positives
      tn <- cm[neg_row, neg_col] # True negatives

      # Calculate TPR and TNR
      tpr <- tp / (tp + fn) # Sensitivity
      tnr <- tn / (tn + fp) # Specificity

      # Calculate G-mean
      gmean <- sqrt(tpr * tnr)
    }
  } else {
    # Handle edge cases (only one class present)
    tpr <- 0
    tnr <- 0
    gmean <- 0
  }

  return(list(gmean = gmean, tpr = tpr, tnr = tnr))
}

# Multiclass evaluation metrics ####

# The 3-class response (nonpeat / otherpeat / bog) has no "positive class", so the
# binary TPR/TNR/G-mean triple does not apply. Reported instead:
#
#  - the full J x J confusion matrix, which is the point of the 3-class reframing: it
#    forces the hard cold contrast (bog vs otherpeat) to be shown separately instead of
#    dissolved into a pooled negative class (notebook 2026-06-18);
#  - per-class recall, precision and one-vs-rest specificity and AUC;
#  - macro G-mean, the geometric mean of the per-class recalls. This is the standard
#    multiclass extension and keeps the property G-mean was chosen for: it goes to zero
#    if ANY class is abandoned, so a model that ignores bog cannot score well by getting
#    the two abundant classes right.
#
# Levels are passed explicitly so a test fold that happens to contain no bog cells still
# produces a full-size matrix rather than a silently smaller one.
calculate_multiclass_metrics <- function(
  predicted_class,
  true_class,
  predicted_prob = NULL,
  levels = c("nonpeat", "otherpeat", "bog")
) {
  predicted_class <- factor(as.character(predicted_class), levels = levels)
  true_class <- factor(as.character(true_class), levels = levels)

  cm <- table(Actual = true_class, Predicted = predicted_class)

  n <- sum(cm)
  by_class <- purrr::map_dfr(levels, function(lvl) {
    tp <- cm[lvl, lvl]
    fn <- sum(cm[lvl, ]) - tp
    fp <- sum(cm[, lvl]) - tp
    tn <- n - tp - fn - fp

    # A class absent from the test fold has no recall to report; NA propagates into the
    # macro summaries rather than being silently scored as 0 or 1.
    recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_

    auc <- NA_real_
    if (!is.null(predicted_prob) && lvl %in% colnames(predicted_prob)) {
      auc <- calculate_ovr_auc(predicted_prob[, lvl], true_class == lvl)
    }

    tibble::tibble(
      class = lvl,
      n_true = tp + fn,
      n_predicted = tp + fp,
      recall = recall,
      precision = precision,
      specificity = specificity,
      auc = auc
    )
  })

  recalls <- by_class$recall

  summary <- tibble::tibble(
    n = n,
    accuracy = sum(diag(cm)) / n,
    balanced_accuracy = mean(recalls, na.rm = TRUE),
    # exp(mean(log(.))) rather than prod(.)^(1/J): the product underflows nothing here,
    # but this form keeps NA handling identical to the mean above.
    Gmean_macro = if (any(recalls == 0, na.rm = TRUE)) {
      0
    } else {
      exp(mean(log(recalls), na.rm = TRUE))
    },
    macro_auc = mean(by_class$auc, na.rm = TRUE)
  )

  list(confusion = cm, by_class = by_class, summary = summary)
}

# One-vs-rest AUC via the Mann-Whitney rank identity. Ties get average ranks, which is
# the standard trapezoidal treatment; no extra package needed.
calculate_ovr_auc <- function(prob, is_positive) {
  is_positive <- as.logical(is_positive)
  keep <- !is.na(prob) & !is.na(is_positive)
  prob <- prob[keep]
  is_positive <- is_positive[keep]

  # Doubles, not integers. sum() on a logical returns an integer, and both n_pos * n_neg
  # and n_pos * (n_pos + 1) overflow R's 32-bit integer past about 46,000 positives --
  # silently, as NA with a warning. That is well inside the range here: binning by bog-DI
  # quantiles leaves the low-novelty bins holding six figures of rows.
  n_pos <- as.numeric(sum(is_positive))
  n_neg <- as.numeric(sum(!is_positive))
  if (n_pos == 0 || n_neg == 0) {
    return(NA_real_)
  }

  r <- rank(prob)
  (sum(r[is_positive]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

# Function to calculate comprehensive metrics at a specific threshold
# Robust to factor level ordering and edge cases
get_metrics_at_threshold <- function(predicted_prob, true_class, threshold) {
  # Convert probabilities to class predictions using threshold
  predicted_class <- ifelse(predicted_prob >= threshold, "1", "0")

  # Convert to character to avoid factor level issues
  predicted_class <- as.character(predicted_class)
  true_class <- as.character(true_class)

  # Create confusion matrix
  cm <- table(Actual = true_class, Predicted = predicted_class)

  if (nrow(cm) == 2 && ncol(cm) == 2) {
    # Explicitly identify which row/column corresponds to positive class ("1")
    pos_row <- which(rownames(cm) == "1")
    pos_col <- which(colnames(cm) == "1")
    neg_row <- which(rownames(cm) == "0")
    neg_col <- which(colnames(cm) == "0")

    # Handle case where positive or negative class is missing
    if (
      length(pos_row) == 0 ||
        length(pos_col) == 0 ||
        length(neg_row) == 0 ||
        length(neg_col) == 0
    ) {
      tpr <- 0
      fpr <- 0
      tnr <- 0
      fnr <- 0
      gmean <- 0
    } else {
      # Extract confusion matrix values using explicit indexing
      tp <- cm[pos_row, pos_col] # True positives
      fn <- cm[pos_row, neg_col] # False negatives
      fp <- cm[neg_row, pos_col] # False positives
      tn <- cm[neg_row, neg_col] # True negatives

      # Calculate all metrics
      tpr <- tp / (tp + fn) # Sensitivity = TPR = 1 - FNR
      tnr <- tn / (tn + fp) # Specificity = TNR = 1 - FPR
      fpr <- fp / (fp + tn) # False positive rate = 1 - TNR
      fnr <- fn / (fn + tp) # False negative rate = 1 - TPR

      # Calculate G-mean
      gmean <- sqrt(tpr * tnr)
    }
  } else {
    # Handle edge cases (only one class present in predictions or true labels)
    # Check which class is missing to provide appropriate metrics
    if (nrow(cm) == 1 || ncol(cm) == 1) {
      # Some calculations may still be possible
      actual_classes <- rownames(cm)
      predicted_classes <- colnames(cm)

      if ("1" %in% actual_classes && "1" %in% predicted_classes) {
        # Only positive class present
        tp <- cm["1", "1"]
        tpr <- 1.0 # All positives correctly predicted
        fnr <- 0.0
        fpr <- NA # No negatives to calculate
        tnr <- NA
      } else if ("0" %in% actual_classes && "0" %in% predicted_classes) {
        # Only negative class present
        tn <- cm["0", "0"]
        tnr <- 1.0 # All negatives correctly predicted
        fpr <- 0.0
        tpr <- NA # No positives to calculate
        fnr <- NA
      } else {
        # Complete mismatch
        tpr <- 0
        fpr <- 1
        tnr <- 0
        fnr <- 1
      }
      gmean <- 0 # Cannot calculate meaningful G-mean
    } else {
      # Completely empty or malformed
      tpr <- 0
      fpr <- 0
      tnr <- 0
      fnr <- 0
      gmean <- 0
    }
  }

  return(data.frame(
    TPR = tpr,
    FPR = fpr,
    TNR = tnr,
    FNR = fnr,
    Gmean = gmean
  ))
}

# Function to calculate classification thresholds using ROCR
# Returns three thresholds: default (G-mean optimal or implicit), sens95, spec95
calculate_classification_thresholds <- function(
  predicted_prob,
  true_labels,
  method = c("gmean", "implicit"),
  implicit_class = NULL,
  target_sensitivity = 0.95,
  target_specificity = 0.95
) {
  method <- match.arg(method)

  # Convert to numeric for ROCR
  true_numeric <- as.numeric(as.character(true_labels))

  # Create ROCR prediction object
  pred_obj <- ROCR::prediction(predicted_prob, true_numeric)

  # Extract performance metrics
  # Note: ROCR orders thresholds from high to low (Inf -> -Inf)
  # TPR increases as threshold decreases (low -> high index)
  # TNR decreases as threshold decreases (low -> high index)
  tpr_values <- ROCR::performance(pred_obj, "tpr")@y.values[[1]]
  tnr_values <- ROCR::performance(pred_obj, "tnr")@y.values[[1]]
  all_thresholds <- ROCR::performance(pred_obj, "tpr")@x.values[[1]]

  # 1. Calculate default threshold based on method
  if (method == "gmean") {
    # G-mean optimal threshold (for BRF)
    gmean_values <- sqrt(tpr_values * tnr_values)
    optimal_idx <- which.max(gmean_values)
    threshold_default <- all_thresholds[optimal_idx]
  } else if (method == "implicit") {
    # Implicit threshold from class predictions (for RFQ)
    if (is.null(implicit_class)) {
      stop("implicit_class must be provided when method = 'implicit'")
    }
    predicted_positives <- implicit_class == "1"
    if (sum(predicted_positives) > 0) {
      threshold_default <- min(predicted_prob[predicted_positives])
    } else {
      threshold_default <- 0.5 # Fallback
    }
  }

  # 2. High sensitivity threshold: TPR >= target_sensitivity
  # Get FIRST threshold where TPR >= target (highest threshold = max specificity)
  sens_candidates <- which(tpr_values >= target_sensitivity)
  if (length(sens_candidates) > 0) {
    threshold_sens <- all_thresholds[sens_candidates[1]]
  } else {
    # If target sensitivity not achievable, use highest sensitivity
    threshold_sens <- all_thresholds[which.max(tpr_values)]
    cat(
      "Warning:",
      target_sensitivity * 100,
      "% sensitivity not achievable, using max TPR =",
      round(max(tpr_values), 3),
      "\n"
    )
  }

  # 3. High specificity threshold: TNR >= target_specificity
  # Get LAST threshold where TNR >= target (lowest threshold = max sensitivity)
  spec_candidates <- which(tnr_values >= target_specificity)
  if (length(spec_candidates) > 0) {
    threshold_spec <- all_thresholds[spec_candidates[length(spec_candidates)]]
  } else {
    # If target specificity not achievable, use highest specificity
    threshold_spec <- all_thresholds[which.max(tnr_values)]
    cat(
      "Warning:",
      target_specificity * 100,
      "% specificity not achievable, using max TNR =",
      round(max(tnr_values), 3),
      "\n"
    )
  }

  return(list(
    threshold_default = threshold_default,
    threshold_sens = threshold_sens,
    threshold_spec = threshold_spec
  ))
}

## Model training functions ####

# Function to train a BRF model with threshold optimization
train_brf_model <- function(
  cv_fold,
  mtry,
  nodesize,
  job_id,
  train_data,
  test_partition = 1
) {
  cat(
    "Job",
    job_id,
    "- BRF Fold",
    cv_fold,
    "mtry:",
    mtry,
    "nodesize:",
    nodesize,
    "...\n"
  )

  # Split data by partition (exclude test_partition from all tuning)
  # train_fold: all partitions except cv_fold and test_partition
  # val_fold: only cv_fold partition
  train_fold <- train_data |>
    filter(partition != test_partition, partition != cv_fold) |>
    select(-partition) |>
    as.data.frame()

  val_fold <- train_data |>
    filter(partition == cv_fold) |>
    select(-partition) |>
    as.data.frame()

  # Convert response to factor with explicit levels (0, 1) for consistency
  train_fold$response <- factor(
    as.character(train_fold$response),
    levels = c("0", "1")
  )
  val_fold$response <- factor(
    as.character(val_fold$response),
    levels = c("0", "1")
  )

  # Ensure ar50 is properly handled as factor if it exists
  if ("ar50" %in% names(train_fold)) {
    train_fold$ar50 <- as.factor(as.character(train_fold$ar50))
    val_fold$ar50 <- as.factor(as.character(val_fold$ar50))
  }

  # Calculate class sizes for balanced sampling
  n_presence <- sum(train_fold$response == "1")
  n_absence <- sum(train_fold$response == "0")
  balanced_size <- min(n_presence, n_absence)

  # Train BRF model with balanced sampling
  model <- randomForest::randomForest(
    formula = response ~ .,
    data = train_fold,
    ntree = 1000,
    mtry = mtry,
    nodesize = nodesize,
    sampsize = c(balanced_size, balanced_size), # Balanced sampling
    replace = TRUE,
    importance = FALSE,
    do.trace = FALSE
  )

  # Get predictions on training fold for threshold optimization
  train_pred_prob <- predict(model, newdata = train_fold, type = "prob")[, "1"]
  train_true <- as.numeric(as.character(train_fold$response))

  # Optimize threshold using ROCR on training folds
  pred_obj <- ROCR::prediction(train_pred_prob, train_true)

  # Calculate G-mean for all possible thresholds
  tpr_values <- ROCR::performance(pred_obj, "tpr")@y.values[[1]]
  tnr_values <- ROCR::performance(pred_obj, "tnr")@y.values[[1]]
  gmean_values <- sqrt(tpr_values * tnr_values)

  # Find optimal threshold
  optimal_idx <- which.max(gmean_values)
  optimal_threshold <- ROCR::performance(pred_obj, "tpr")@x.values[[1]][
    optimal_idx
  ]

  # Get predictions on validation fold
  val_pred_prob <- predict(model, newdata = val_fold, type = "prob")[, "1"]
  val_pred_class <- ifelse(val_pred_prob >= optimal_threshold, "1", "0")
  val_true <- as.character(val_fold$response)

  # Calculate G-mean on validation set with optimal threshold
  metrics <- calculate_gmean(val_pred_class, val_true)
  gmean <- metrics$gmean
  tpr <- metrics$tpr
  tnr <- metrics$tnr

  # Create result row
  result <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    model_type = "BRF",
    cv_fold = cv_fold,
    mtry = mtry,
    nodesize = nodesize,
    gmean = gmean,
    tpr = tpr,
    tnr = tnr,
    stringsAsFactors = FALSE
  )

  # Append to CSV immediately
  write_csv(result, progress_file, append = TRUE)
  cat("  Job", job_id, "completed - G-mean:", round(gmean, 4), "\n")

  return(result)
}

# Function to train a RFQ model (existing approach)
train_rfq_model <- function(
  cv_fold,
  mtry,
  nodesize,
  job_id,
  train_data,
  test_partition = 1
) {
  cat(
    "Job",
    job_id,
    "- RFQ Fold",
    cv_fold,
    "mtry:",
    mtry,
    "nodesize:",
    nodesize,
    "...\n"
  )

  # Split data by partition (exclude test_partition from all tuning)
  # train_fold: all partitions except cv_fold and test_partition
  # val_fold: only cv_fold partition
  train_fold <- train_data |>
    filter(partition != test_partition, partition != cv_fold) |>
    select(-partition) |>
    as.data.frame() # Convert tibble to data.frame for imbalanced()

  val_fold <- train_data |>
    filter(partition == cv_fold) |>
    select(-partition) |>
    as.data.frame() # Convert tibble to data.frame for imbalanced()

  # Convert response to factor with explicit levels (0, 1) for consistency
  train_fold$response <- factor(
    as.character(train_fold$response),
    levels = c("0", "1")
  )
  val_fold$response <- factor(
    as.character(val_fold$response),
    levels = c("0", "1")
  )

  # Ensure ar50 is properly handled as factor if it exists
  if ("ar50" %in% names(train_fold)) {
    train_fold$ar50 <- as.factor(as.character(train_fold$ar50))
    val_fold$ar50 <- as.factor(as.character(val_fold$ar50))
  }

  # Train model with current parameters
  model <- randomForestSRC::imbalanced(
    formula = response ~ .,
    data = train_fold,
    ntree = 3000,
    mtry = mtry,
    nodesize = nodesize,
    importance = FALSE, # Skip importance for speed during tuning
    do.trace = FALSE
  )

  # Predict on validation set
  pred <- predict(model, newdata = val_fold)

  # Calculate G-mean on validation set
  pred_class <- pred$class
  true_class <- val_fold$response

  # Calculate G-mean using common function
  metrics <- calculate_gmean(pred_class, true_class)
  gmean <- metrics$gmean
  tpr <- metrics$tpr
  tnr <- metrics$tnr

  # Create result row
  result <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    model_type = "RFQ",
    cv_fold = cv_fold,
    mtry = mtry,
    nodesize = nodesize,
    gmean = gmean,
    tpr = tpr,
    tnr = tnr,
    stringsAsFactors = FALSE
  )

  # Append to CSV immediately
  write_csv(result, progress_file, append = TRUE)
  cat("  Job", job_id, "completed - G-mean:", round(gmean, 4), "\n")

  return(result)
}

# Unified function to train either model type
train_single_model <- function(
  cv_fold,
  model_type,
  mtry,
  nodesize,
  job_id,
  train_data,
  test_partition = 1
) {
  if (model_type == "BRF") {
    return(train_brf_model(
      cv_fold,
      mtry,
      nodesize,
      job_id,
      train_data,
      test_partition
    ))
  } else if (model_type == "RFQ") {
    return(train_rfq_model(
      cv_fold,
      mtry,
      nodesize,
      job_id,
      train_data,
      test_partition
    ))
  } else {
    stop("Unknown model_type: ", model_type)
  }
}

# Presence-based partitioning along sorted feature values ####
# Creates k partitions by sorting presences along a feature and dividing into groups
# Non-presence rows are assigned to partitions based on feature value overlap
# Rows in gaps between partitions are dropped
#
# The response became 3-class with the EU integration (nonpeat / otherpeat / bog), so
# "presence" is now named explicitly rather than implied by `== 1`. Only the *anchor*
# class partitions the feature axis -- bog is both the rare class and the one the
# partitions exist to spread out -- and every other class is assigned by range overlap,
# exactly as the two-class version assigned absences. Class identity is untouched: this
# function decides which partition a row belongs to, never what it is.
partition_by_presence_sorting <- function(
  data,
  feature_name,
  k = 6,
  min_pres = 50,
  seed = NULL,
  presence_level = "bog"
) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Extract response and feature values
  Y <- as.character(data$response)
  X <- data[[feature_name]]

  if (!presence_level %in% Y) {
    stop(
      "presence_level '", presence_level, "' not found in the response.",
      "\n  Levels present: ", paste(sort(unique(Y)), collapse = ", ")
    )
  }

  # Identify anchor (presence) and non-anchor indices
  pres_idx <- which(Y == presence_level)
  abs_idx <- which(Y != presence_level)

  n_pres <- length(pres_idx)
  n_abs <- length(abs_idx)

  # Check minimum presence constraint
  if (n_pres < k * min_pres) {
    stop(
      "Insufficient presences for ", k, " partitions with min_pres = ", min_pres,
      "\n  Available presences: ", n_pres,
      "\n  Required: ", k * min_pres
    )
  }

  # STEP 1: Partition the presences ####
  # Sort presences by feature value
  pres_X <- X[pres_idx]
  pres_order <- order(pres_X)
  pres_sorted <- pres_idx[pres_order]

  # Divide into k groups of roughly equal size
  group_size <- floor(n_pres / k)
  remainder <- n_pres %% k

  pres_groups <- vector("list", k)
  start <- 1
  for (i in 1:k) {
    # Distribute remainder across first groups
    this_size <- group_size + ifelse(i <= remainder, 1, 0)
    pres_groups[[i]] <- pres_sorted[start:(start + this_size - 1)]
    start <- start + this_size
  }

  # STEP 2: Define partition boundaries ####
  # For each presence group, find X-range
  partitions <- vector("list", k)
  for (i in 1:k) {
    pres_in_group <- pres_groups[[i]]
    x_min <- min(X[pres_in_group])
    x_max <- max(X[pres_in_group])

    partitions[[i]] <- list(
      x_lower = x_min,
      x_upper = x_max,
      presence_idx = pres_in_group,
      absence_idx = integer(0)
    )
  }

  # STEP 3: Assign non-anchor rows ####
  # Each goes to the partition whose X-range contains it; rows in gaps are dropped.
  # Done as k vectorized passes rather than a per-row scan: the pooled frame has ~83k
  # non-anchor rows, and first-match-wins is preserved by only ever assigning rows that
  # are still unassigned -- the same rule as the earlier `break`.
  abs_assignment <- rep(NA_integer_, length(abs_idx))
  abs_X <- X[abs_idx]

  for (i in 1:k) {
    in_range <- is.na(abs_assignment) &
      abs_X >= partitions[[i]]$x_lower &
      abs_X <= partitions[[i]]$x_upper
    abs_assignment[in_range] <- i
  }

  for (i in 1:k) {
    partitions[[i]]$absence_idx <- abs_idx[which(abs_assignment == i)]
  }

  dropped_abs <- abs_idx[which(is.na(abs_assignment))]

  # Create partition assignment vector
  partition_assignment <- rep(NA_integer_, nrow(data))
  for (i in 1:k) {
    all_idx <- c(partitions[[i]]$presence_idx, partitions[[i]]$absence_idx)
    partition_assignment[all_idx] <- i
  }

  # Return results
  list(
    partitions = partition_assignment,
    partition_list = partitions,
    n_presences = sapply(partitions, function(p) length(p$presence_idx)),
    n_absences = sapply(partitions, function(p) length(p$absence_idx)),
    # Full class breakdown per partition: with three classes, "presences and absences"
    # no longer describes the table a reader needs to check the partitions against.
    n_by_class = table(partition = partition_assignment, class = Y, useNA = "no"),
    n_dropped = length(dropped_abs),
    dropped_idx = dropped_abs,
    dropped_by_class = table(Y[dropped_abs]),
    feature_name = feature_name,
    presence_level = presence_level,
    k = k
  )
}

# Shared absence sampling for the pooled EU + Norway block ####

# One rule, stated once, applied to both blocks (plan_EUintegration.md section 3 rule 11
# and section 5). Two properties the EU block was recruited for depend on it:
#
#  - Climate stratification, never uniform-random. Non-peat supply is effectively
#    unlimited in both domains, so a uniform draw puts almost every absence in the cool
#    wet core where Norway already has data and almost none on the warm flank. That exact
#    failure is already logged for pl0_modelGlobalScale.R.
#  - Symmetry between blocks. Norway and Europe are sampled by the same function with the
#    same bins, so "other-peat at 18 degrees" means the same thing on both sides and the
#    dataset-blocked comparison of section 1.5 compares like with like.

# Fixed 1-degree bio10 bins with a fixed origin, spanning both domains (Norway 1.8-17.9,
# EU 6.1-21.0). Fixed rather than per-block quantiles so a stratum index is comparable
# across blocks and across re-runs.
BIO10_BREAKS <- seq(0, 25, by = 1)

# Absences drawn per bog cell, per absence class. THE swept parameter of the sampling
# design (section 1.5: any weight that survives must be swept, not tuned). Equal targets
# per absence class deliberately depart from the natural class shares: raised bog is a
# climatic intermediate flanked by two ecologically opposite classes, so both flanks need
# real sample, and an area-proportional draw would leave other-peat too thin on the EU
# side to carry the contrast. This is stated, symmetric and swept -- not silent.
ABSENCES_PER_BOG_PER_CLASS <- 25

# Water-filling allocation across strata: an equal target per stratum, capped by what the
# stratum actually holds, with the remainder redistributed over strata that still have
# room. One-shot proportional allocation would let a scarce warm stratum quietly absorb
# the rounding error of the whole class budget.
allocate_across_strata <- function(available, total, floor_idx = integer(0),
                                   floor_frac = 0.5) {
  alloc <- rep(0L, length(available))
  remaining <- total
  open <- available > 0

  # Floor pass: secure the warm strata before general water-filling, because they are
  # the scarce ones and the reason the block exists.
  if (length(floor_idx) > 0 && any(open)) {
    want <- floor(floor(total / sum(open)) * floor_frac)
    for (i in intersect(floor_idx, which(open))) {
      take <- min(want, available[i], remaining)
      alloc[i] <- alloc[i] + take
      remaining <- remaining - take
    }
  }

  while (remaining > 0 && any(available - alloc > 0)) {
    room <- which(available - alloc > 0)
    per <- max(1, floor(remaining / length(room)))
    for (i in room) {
      if (remaining <= 0) break
      take <- min(per, available[i] - alloc[i], remaining)
      alloc[i] <- alloc[i] + take
      remaining <- remaining - take
    }
  }
  alloc
}

# pool: data frame with columns `response` (class label) and `stratum` (index into
#       BIO10_BREAKS), plus whatever else should travel with the drawn rows.
# targets: named vector of per-class draw targets.
# warm_strata: stratum indices above the warm-flank threshold, guaranteed a floor.
draw_stratified_absences <- function(pool, targets, warm_strata,
                                     floor_frac = 0.5) {
  n_strata <- length(BIO10_BREAKS) - 1L

  per_class <- lapply(names(targets), function(cls) {
    sub <- pool[pool$response == cls, , drop = FALSE]
    available <- tabulate(sub$stratum, nbins = n_strata)
    alloc <- allocate_across_strata(available, targets[[cls]], warm_strata, floor_frac)

    rows <- do.call(rbind, lapply(seq_len(n_strata), function(s) {
      if (alloc[s] == 0) {
        return(NULL)
      }
      idx <- which(sub$stratum == s)
      sub[sample(idx, alloc[s]), , drop = FALSE]
    }))

    list(
      rows = rows,
      allocation = data.frame(
        response = cls,
        stratum = seq_len(n_strata),
        bio10_from = BIO10_BREAKS[-(n_strata + 1L)],
        bio10_to = BIO10_BREAKS[-1L],
        warm = seq_len(n_strata) %in% warm_strata,
        available = available,
        drawn = alloc
      )
    )
  })

  list(
    rows = do.call(rbind, lapply(per_class, `[[`, "rows")),
    allocation = do.call(rbind, lapply(per_class, `[[`, "allocation"))
  )
}
