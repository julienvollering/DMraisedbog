# Compare feature space distances between observations and modeling domain

library(tidyverse)
library(terra)
library(sf)

# Install and load forked CAST

remotes::install_github("julienvollering/CAST", 
                        build_manual = TRUE,
                        build_vignettes = FALSE)
library(CAST)

# Streamlined feature distance function
featuredist <- function(training_data, prediction_data, variable_name,
                        cvfold_column = NULL, sample_size = 1000,
                        seed = NULL, normalize = FALSE) {
  
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
  
  # Normalize if requested
  if (normalize) {
    train_mean <- mean(train_clean)
    train_sd <- sd(train_clean)
    if (train_sd > 0) {
      train_clean <- (train_clean - train_mean) / train_sd
      pred_clean <- (pred_clean - train_mean) / train_sd
    }
  }
  
  # Convert to matrices for FNN
  train_matrix <- matrix(train_clean, ncol = 1)
  pred_matrix <- matrix(pred_clean, ncol = 1)
  
  # Calculate sample-to-sample distances
  s2s_dists <- numeric(length(train_clean))
  for (i in seq_along(train_clean)) {
    dists <- FNN::knnx.dist(train_matrix[i, , drop = FALSE], train_matrix, k = 1)
    dists[i] <- NA  # Exclude self-distance
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
    what = factor(c(rep("sample-to-sample", length(s2s_dists)),
                    rep("prediction-to-sample", length(s2p_dists))),
                  levels = c("sample-to-sample", "prediction-to-sample",
                             "CV-distances")),
    dist_type = "feature"
  )
  
  # Add CV distances if fold column is provided
  if (!is.null(cvfold_column)) {
    cv_folds <- train_sampled[[cvfold_column]]
    cv_folds <- cv_folds[!is.na(train_var)]  # Match cleaned training data
    
    cv_dists <- numeric(0)
    unique_folds <- unique(cv_folds)
    
    for (fold in unique_folds) {
      test_idx <- which(cv_folds == fold)
      train_idx <- which(cv_folds != fold)
      
      if (length(test_idx) > 0 && length(train_idx) > 0) {
        test_matrix <- matrix(train_clean[test_idx], ncol = 1)
        fold_train_matrix <- matrix(train_clean[train_idx], ncol = 1)
        
        for (i in seq_along(test_idx)) {
          dists <- FNN::knnx.dist(test_matrix[i, , drop = FALSE], fold_train_matrix, k = 1)
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

## Read modeling frame and domain raster ####
mf <- read_csv("output/pl2/modeling_frame_regional.csv")

raster_future <- rast("output/pl2/scenario_future.tif")

## Load feature weights ####
weights_features <- read_csv("output/pl2/weights_feature_data_partitioning.csv")
weights_features <- weights_features |>
  mutate(feature = forcats::fct_reorder(
    feature, median,
    .fun = first
  ))

weights <- weights_features |>
  select(feature, median) |>
  arrange(desc(median)) |>
  tibble::deframe()

## Feature spaces ####

mf_current <- mf |> 
  filter(scenario == "current") |>
  select(-scenario, -response) |> 
  st_as_sf(coords = c("x", "y"), crs = 3035)

mf_current_presence <- mf |> 
  filter(scenario == "current", response == 1) |>
  select(-scenario, -response) |> 
  st_as_sf(coords = c("x", "y"), crs = 3035)

mf_future <- mf |> 
  filter(scenario == "future") |>
  select(-scenario, -response) |> 
  st_as_sf(coords = c("x", "y"), crs = 3035)

toJoin <- mf |> 
  filter(scenario == "future") |> 
  select(-scenario, -response)
mf_future_presence <- mf |> 
  filter(response == 1) |> 
  select(x, y) |> 
  left_join(toJoin, by = c("x", "y")) |> 
  st_as_sf(coords = c("x", "y"), crs = 3035)

## Distances ####

#### Test when sample == prediction ####

dists <- CAST:::geodist.fast(
  sample = mf_current_presence[,"elevation"],
  prediction = mf_future_presence[,"elevation"],
  variables = "elevation",
  normalize = FALSE,
  weights = NULL,
  samplesize = 1e4) # Includes all
plot(dists, type = "simple")
dists |> 
  group_by(what) |>
  summarise(dist = median(dist))
# Sample-to-sample do not include NN to self, while prediction-to-sample do

### Domain = presences x future ####

#### 1 feature
x <- mf_current[,names(weights)[1]]
y <- mf_future_presence[,names(weights)[1]]
ggplot() +
  geom_density(aes(x = pull(x, names(weights)[1])), color = "blue") +
  geom_density(aes(x = pull(y, names(weights)[1])), color = "red") +
  labs(title = names(weights)[1], x = "Feature value", y = "Density") +
  theme_minimal()
dists1 <- CAST:::geodist.fast(
  sample = x, # Full geographic extent of current
  prediction = y,
  variables = names(weights)[1],
  normalize = FALSE,
  weights = NULL,
  samplesize = 1e4) # Increase samplesize for stability 
plot(dists1, type = "simple")
dists1 |> 
  group_by(what) |>
  summarise(dist = mean(dist))

#### 2 features
x <- mf_current[,names(weights)[1:2]]
y <- mf_future_presence[,names(weights)[1:2]]
f1 <- ggplot() +
  geom_density(aes(x = pull(x, names(weights)[1])), color = "blue") +
  geom_density(aes(x = pull(y, names(weights)[1])), color = "red") +
  labs(title = names(weights)[1], x = "Feature value", y = "Density") +
  theme_minimal()
f2 <- ggplot() +
  geom_density(aes(x = pull(x, names(weights)[2])), color = "blue") +
  geom_density(aes(x = pull(y, names(weights)[2])), color = "red") +
  labs(title = names(weights)[2], x = "Feature value", y = "Density") +
  theme_minimal()
patchwork::wrap_plots(f1, f2, ncol = 1)
dists2 <- CAST:::geodist.fast(
  sample = mf_current[,names(weights)[1:2]], # Full geographic extent of current
  prediction = mf_future_presence[,names(weights)[1:2]],
  variables = names(weights)[1:2],
  normalize = FALSE,
  weights = NULL,
  samplesize = 1e4) # Increase samplesize for stability 
plot(dists2, type = "simple")
dists2 |> 
  group_by(what) |>
  summarise(dist = median(dist))

#### 3 features
x <- mf_current[,names(weights)[1:3]]
y <- mf_future_presence[,names(weights)[1:3]]
f1 <- ggplot() +
  geom_density(aes(x = pull(x, names(weights)[1])), color = "blue" ) +
  geom_density(aes(x = pull(y, names(weights)[1])), color = "red") +
  labs(title = names(weights)[1], x = "Feature value", y = "Density") +
  theme_minimal()
f2 <- ggplot() +
  geom_density(aes(x = pull(x, names(weights)[2])), color = "blue") +
  geom_density(aes(x = pull(y, names(weights)[2])), color = "red") +
  labs(title = names(weights)[2], x = "Feature value", y = "Density") +
  theme_minimal()
f3 <- ggplot() +
  geom_density(aes(x = pull(x, names(weights)[3])), color = "blue") +
  geom_density(aes(x = pull(y, names(weights)[3])), color = "red") +
  labs(title = names(weights)[3], x = "Feature value", y = "Density") +
  theme_minimal()
patchwork::wrap_plots(f1, f2, f3, ncol = 1)
f3
dists3 <- CAST:::geodist.fast(
  sample = mf_current[,names(weights)[1:3]], # Full geographic extent of current
  prediction = mf_future_presence[,names(weights)[1:3]],
  variables = names(weights)[1:3],
  normalize = FALSE,
  weights = NULL,
  samplesize = 1e4) # Increase samplesize for stability 
plot(dists3, type = "simple")
dists3 |> 
  group_by(what) |>
  summarise(dist = median(dist))

### With CV ####  

##### Maximally separated by top feature ####
mf_current_cut <- mf_current |> 
  select(names(weights)[1]) |>
  mutate(assignment = cut(
    pull(mf_current, names(weights)[1]),
    breaks = quantile(pull(mf_current, names(weights)[1]), probs = seq(0, 1, 0.2)),
    include.lowest = TRUE
  ))
x <- mf_current_cut
y <- mf_future_presence[,names(weights)[1]]
g1 <- ggplot() +
  geom_density(aes(x = pull(x, names(weights)[1])), color = "blue") +
  geom_density(aes(x = pull(y, names(weights)[1])), color = "red") +
  labs(title = names(weights)[1], x = "Feature value", y = "Density") +
  coord_cartesian(xlim = c(200, 600)) +
  theme_minimal()
cv <- ggplot(mf_current_cut) +
  geom_density(aes(x = pull(mf_current_cut, names(weights)[1]), fill = assignment), alpha = 0.5) +
  labs(title = names(weights)[1], x = "Feature value", y = "Density") +
  coord_cartesian(xlim = c(200, 600)) +
  theme_minimal()
patchwork::wrap_plots(g1, cv, ncol = 1)
dists <- CAST:::geodist.fast(
  sample = mf_current[,names(weights)[1]],
  prediction = mf_future_presence[,names(weights)[1]],
  cvfolds = as.numeric(mf_current_cut$assignment),
  variables = names(weights)[1],
  normalize = FALSE,
  weights = NULL,
  samplesize = 1e4) # Increase samplesize for stability 
dists |> 
  group_by(what) |>
  summarise(dist = median(dist))
ggplot(dists) +
  geom_density(aes(x = dist, color = what)) +
  theme_minimal()

##### Unseparated CV folds ####
dists <- CAST:::geodist.fast(
  sample = mf_current[,names(weights)[1]],
  prediction = mf_future_presence[,names(weights)[1]],
  cvfolds = rep(1:5, length.out = nrow(mf_current)),
  variables = names(weights)[1],
  normalize = FALSE,
  weights = NULL,
  samplesize = 1e4) # Increase samplesize for stability 
dists |> 
  group_by(what) |>
  summarise(dist = median(dist))

##### Maximally separated by 2rd feature ####
mf_current_cut <- mf_current |> 
  select(names(weights)[2]) |>
  mutate(assignment = cut(
    pull(mf_current, names(weights)[2]),
    breaks = quantile(pull(mf_current, names(weights)[2]), probs = seq(0, 1, 0.2)),
    include.lowest = TRUE
  ))
x <- mf_current_cut
y <- mf_future_presence[,names(weights)[2]]
g1 <- ggplot() +
  geom_density(aes(x = pull(x, names(weights)[2])), color = "blue") +
  geom_density(aes(x = pull(y, names(weights)[2])), color = "red") +
  labs(title = names(weights)[2], x = "Feature value", y = "Density") +
  coord_cartesian(xlim = c(200, 1100)) +
  theme_minimal()
cv <- ggplot(mf_current_cut) +
  geom_density(aes(x = pull(mf_current_cut, names(weights)[2]), fill = assignment), alpha = 0.5) +
  labs(title = names(weights)[2], x = "Feature value", y = "Density") +
  coord_cartesian(xlim = c(200, 1100)) +
  theme_minimal()
patchwork::wrap_plots(g1, cv, ncol = 1)
dists <- CAST:::geodist.fast(
  sample = mf_current[,names(weights)[2]],
  prediction = mf_future_presence[,names(weights)[2]],
  cvfolds = as.numeric(mf_current_cut$assignment),
  variables = names(weights)[2],
  normalize = FALSE,
  weights = NULL,
  samplesize = 1e4) # Increase samplesize for stability 
dists |> 
  group_by(what) |>
  summarise(dist = median(dist))
ggplot(dists) +
  geom_density(aes(x = dist, color = what)) +
  theme_minimal()

##### Maximally separated by 3rd feature ####
mf_current_cut <- mf_current |> 
  select(names(weights)[3]) |>
  mutate(assignment = cut(
    pull(mf_current, names(weights)[3]),
    breaks = quantile(pull(mf_current, names(weights)[3]), probs = seq(0, 1, 0.2)),
    include.lowest = TRUE
  ))
x <- mf_current_cut
y <- mf_future_presence[,names(weights)[3]]
g1 <- ggplot() +
  geom_density(aes(x = pull(x, names(weights)[3])), color = "blue") +
  geom_density(aes(x = pull(y, names(weights)[3])), color = "red") +
  labs(title = names(weights)[3], x = "Feature value", y = "Density") +
  coord_cartesian(xlim = c(0, 30)) +
  theme_minimal()
cv <- ggplot(mf_current_cut) +
  geom_density(aes(x = pull(mf_current_cut, names(weights)[3]), fill = assignment), alpha = 0.5) +
  labs(title = names(weights)[3], x = "Feature value", y = "Density") +
  coord_cartesian(xlim = c(0, 30)) +
  theme_minimal()
patchwork::wrap_plots(g1, cv, ncol = 1)
dists <- CAST:::geodist.fast(
  sample = mf_current[,names(weights)[3]],
  prediction = mf_future_presence[,names(weights)[3]],
  cvfolds = as.numeric(mf_current_cut$assignment),
  variables = names(weights)[3],
  normalize = FALSE,
  weights = NULL,
  samplesize = 1e4) # Increase samplesize for stability 
dists |> 
  group_by(what) |>
  summarise(dist = median(dist))
ggplot(dists) +
  geom_density(aes(x = dist, color = what)) +
  theme_minimal()

## Distances w/ `featuredist` ####

#### Test when sample == prediction ####

dists <- featuredist(
  mf_current_presence[,"elevation"],
  mf_future_presence[,"elevation"],
  variable_name = "elevation",
  normalize = FALSE,
  sample_size = 1e4)
plot(dists, type = "simple")
dists |> 
  group_by(what) |>
  summarise(dist = median(dist))
# Sample-to-sample do not include NN to self, while prediction-to-sample do

### With CV ####  

##### Maximally separated by top feature ####
mf_current_cut <- mf_current |> 
  select(names(weights)[1]) |>
  mutate(assignment = cut(
    pull(mf_current, names(weights)[1]),
    breaks = quantile(pull(mf_current, names(weights)[1]), probs = seq(0, 1, 0.2)),
    include.lowest = TRUE
  ))
x <- mf_current_cut
y <- mf_future_presence[,names(weights)[1]]
g1 <- ggplot() +
  geom_density(aes(x = pull(x, names(weights)[1])), color = "blue") +
  geom_density(aes(x = pull(y, names(weights)[1])), color = "red") +
  labs(title = names(weights)[1], x = "Feature value", y = "Density") +
  coord_cartesian(xlim = c(200, 600)) +
  theme_minimal()
cv <- ggplot(mf_current_cut) +
  geom_density(aes(x = pull(mf_current_cut, names(weights)[1]), fill = assignment), alpha = 0.5) +
  labs(title = names(weights)[1], x = "Feature value", y = "Density") +
  coord_cartesian(xlim = c(200, 600)) +
  theme_minimal()
patchwork::wrap_plots(g1, cv, ncol = 1)
dists <- featuredist(x, y, variable_name = names(weights)[1], 
                     cvfold_column = "assignment")
dists |> 
  group_by(what) |>
  summarise(dist = median(dist))
ggplot(dists) +
  geom_density(aes(x = dist, color = what)) +
  theme_minimal()
plot(dists, stat = "ecdf")
attr(dists, "W_CV")

##### Unseparated CV folds ####
x <- mf_current_cut |> 
  mutate(assignment2 = rep(1:5, length.out = nrow(mf_current)))
dists <- featuredist(x, y, variable_name = names(weights)[1], 
                     cvfold_column = "assignment2")
dists |> 
  group_by(what) |>
  summarise(dist = median(dist))
plot(dists, stat = "ecdf")
attr(dists, "W_CV")

# sessionInfo ####

sessioninfo::session_info()
