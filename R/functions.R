## Functions ####

# Input freshness ####

# Stops a script that is about to read an input older than the file it was derived from.
# On 2026-08-24 the two reliability scripts ran on a predictions file that predated the
# partitioned frame it should have been computed from, and produced a full set of numbers
# that read as results. A modification-time check is a crude dependency test, but it is
# exactly the test that would have caught that: `newer` must have been written after
# `older`. Both paths are named in the error so the fix (rerun the producer of `newer`)
# is obvious.
assert_fresher <- function(newer, older) {
  stopifnot(file.exists(newer), file.exists(older))
  t_new <- file.mtime(newer)
  t_old <- file.mtime(older)
  if (t_new < t_old) {
    stop(
      "Stale input: ", newer, " (", format(t_new, "%Y-%m-%d %H:%M:%S"), ") is older than ",
      older, " (", format(t_old, "%Y-%m-%d %H:%M:%S"), "). Rerun the script that writes ",
      basename(newer), " before this one."
    )
  }
  invisible(TRUE)
}

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

  # STEP 2: Define the presence envelope and the assignment boundaries ####
  #
  # These are two different things, and conflating them is what used to discard a quarter
  # of the frame. Each group's presence extremes describe where the ANCHORS of that group
  # lie; they are kept because they define the envelope and are what `envelope_side` is
  # measured against. They are NOT used to assign non-anchor rows, because outside the
  # outermost anchors there is no interval to fall into, and between two groups there is a
  # gap. Rows landing in either were previously dropped -- on the pooled frame that was
  # 20,447 of 82,359 rows, including every non-anchor beyond the warm end of the anchor
  # range, which is precisely the region a climate projection travels into.
  #
  # The assignment boundaries instead TILE the axis completely: half-open intervals cut at
  # the midpoint between adjacent groups, with the terminal intervals open to +/- Inf.
  # Every row therefore lands in exactly one partition and nothing is discarded.
  partitions <- vector("list", k)
  for (i in 1:k) {
    pres_in_group <- pres_groups[[i]]
    partitions[[i]] <- list(
      x_lower = min(X[pres_in_group]),   # presence extremes: the envelope, not the cut
      x_upper = max(X[pres_in_group]),
      presence_idx = pres_in_group,
      absence_idx = integer(0)
    )
  }

  # Interior cuts at midpoints between adjacent groups; ends open.
  interior_cuts <- vapply(
    seq_len(k - 1),
    function(i) mean(c(partitions[[i]]$x_upper, partitions[[i + 1]]$x_lower)),
    numeric(1)
  )
  boundaries <- c(-Inf, interior_cuts, Inf)
  for (i in 1:k) {
    partitions[[i]]$cut_lower <- boundaries[i]
    partitions[[i]]$cut_upper <- boundaries[i + 1]
  }

  # The presence envelope, against which position is reported.
  envelope_min <- min(pres_X)
  envelope_max <- max(pres_X)

  # STEP 3: Assign non-anchor rows ####
  # findInterval over the interior cuts rather than k range tests: the intervals are
  # disjoint and exhaustive by construction, so there is no first-match-wins rule left to
  # preserve and no row can fail to match. Anchors keep the group they were sorted into in
  # STEP 1 and are never reassigned here -- with tied feature values an anchor can sit on
  # the far side of a midpoint cut from its own group, and the sort is what defines the
  # equal-count property the partitions exist for.
  abs_X <- X[abs_idx]
  abs_assignment <- findInterval(abs_X, interior_cuts) + 1L
  abs_assignment[is.na(abs_X)] <- NA_integer_

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

  # Position relative to the presence envelope, named for the axis rather than for what
  # the axis measures -- this function is generic over `feature_name`. Every row carries
  # it, so any metric can be reported inside-envelope and tail-inclusive without refitting.
  envelope_side <- rep(NA_character_, nrow(data))
  envelope_side[!is.na(X)] <- "inside"
  envelope_side[!is.na(X) & X < envelope_min] <- "below"
  envelope_side[!is.na(X) & X > envelope_max] <- "above"

  # Return results
  list(
    partitions = partition_assignment,
    partition_list = partitions,
    envelope_side = envelope_side,
    envelope_range = c(lower = envelope_min, upper = envelope_max),
    boundaries = boundaries,
    n_presences = sapply(partitions, function(p) length(p$presence_idx)),
    n_absences = sapply(partitions, function(p) length(p$absence_idx)),
    # Full class breakdown per partition: with three classes, "presences and absences"
    # no longer describes the table a reader needs to check the partitions against.
    n_by_class = table(partition = partition_assignment, class = Y, useNA = "no"),
    n_by_envelope_side = table(
      partition = partition_assignment,
      side = envelope_side,
      useNA = "no"
    ),
    # Retained so callers that report them keep working. With complete tiling these are
    # empty unless the feature itself is NA, which makes a non-zero count a real signal.
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
