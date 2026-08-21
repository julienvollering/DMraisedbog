library(tidyverse)
library(terra)
library(sf)

# Install and load forked CAST

remotes::install_github("julienvollering/CAST", 
                        build_manual = TRUE,
                        build_vignettes = FALSE)
library(CAST)
?geodist
?knndm

# Load the training data ####
rf_global_current <- rast("output/rf_global_pred_regional_current.tif") 
names(rf_global_current) <- "rf_global"
chelsa_nor_250m_cur <- rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")[[1:19]]
ar50_nor_250m <- rast("output/ar50_250m_EPSG3035.tif")
names(ar50_nor_250m) <- "ar50"
ar50_nor_250m <- as.factor(ar50_nor_250m) # Convert to factor
preds_train <- c(rf_global_current, chelsa_nor_250m_cur, ar50_nor_250m)

# Load the projection data
rf_global_future <- rast("output/rf_global_pred_regional_future.tif")
names(rf_global_future) <- "rf_global"
chelsa_nor_250m_fut <- rast("output/predictors_regional_250m_Norway_future_EPSG3035.tif")[[1:19]]
preds_proj <- c(rf_global_future, chelsa_nor_250m_fut, ar50_nor_250m)

presence <- read_csv("output/presence_coords_regional.csv")
absence <- read_csv("output/absence_coords_regional.csv")

# Create training dataset ####

# Extract predictor values for presence points
presence_df <- terra::extract(preds_train, presence[c("x", "y")], xy=TRUE, ID=FALSE) |> 
  mutate(response = 1, .before = 1) |>
  drop_na() |> 
  tibble()

# Extract predictor values for absence points
absence_df <- terra::extract(preds_train, absence[c("x", "y")], xy=TRUE, ID=FALSE) |> 
  mutate(response = 0, .before = 1) |>
  drop_na() |> 
  tibble()

# Combine into training dataset
training_data <- bind_rows(presence_df, absence_df) |>
  mutate(response = factor(response), ar50 = factor(ar50)) |> 
  select(-x, -y) # Remove x and y coordinates

# Get final counts
n_final_presence <- sum(training_data$response == "1")
n_final_background <- sum(training_data$response == "0")

cat("Final training data - Presences:", n_final_presence, "Background:", n_final_background, "\n")
cat("Prevalence in training data:", 
    round(n_final_presence / (n_final_presence + n_final_background) * 100, 2), "%\n")

# Feature weights for feature distances ####

# library(ranger)  # Faster for repeated runs
# 
# get_stable_weights <- function(outcome, predictors, n_reps = 10) {
#   importance_matrix <- matrix(0, n_reps, ncol(predictors))
#   
#   # Calculate the fraction needed for balanced sampling
#   n_pos <- sum(outcome == 1)
#   n_neg <- sum(outcome == 0)
#   
#   # For balanced: sample ALL positives, and equal number of negatives
#   frac_pos <- 1.0  # Use all positive cases
#   frac_neg <- n_pos / n_neg  # Downsample negatives to match
#   
#   for(i in 1:n_reps) {
#     rf <- ranger(outcome ~ ., 
#                  data = cbind(predictors, outcome = as.factor(outcome)),
#                  num.trees = 500,
#                  importance = "impurity",
#                  sample.fraction = c(frac_neg, frac_pos),  # Order matters! 
#                  replace = TRUE,
#                  oob.error = TRUE)
#     
#     importance_matrix[i, ] <- rf$variable.importance
#   }
#   
#   weights <- apply(importance_matrix, 2, median)
#   return(weights / max(weights))
# }
# 
# tictoc::tic()
# weights_gini <- get_stable_weights(training_data$response, 
#                                    training_data[, -1], 
#                                    n_reps=2)
# tictoc::toc() #2525 sec
# names(weights_gini) <- names(training_data)[-1] # Exclude response column
# as.data.frame(weights_gini) |> 
#   rownames_to_column("feature") |> 
#   write_csv("output/feature_weights_gini.csv")

weights_gini <- read_csv("output/feature_weights_gini.csv")
weights_gini <- weights_gini |> 
  pivot_wider(names_from = feature, values_from = weights_gini) #argument format

# CAST functions ####

training_data <- bind_rows(presence_df, absence_df) |>
  mutate(response = factor(response), ar50 = factor(ar50))

training_points <- 
  training_data |> 
  select(-response) |> 
  st_as_sf(coords = c("x", "y"), crs = 3035)

## geodist() ####

?geodist
dist_train_proj_2features <- geodist(
  x = training_points[sample.int(nrow(training_points), 1000),1:2], # From docs: Predictor values for x, testdata and preddata are optional if modeldomain is a raster. If not provided they are extracted from the modeldomain rasterStack.
  modeldomain = preds_proj, # From docs: If type = "feature", the argument modeldomain...has to include predictors
  type = "feature",
  samplesize = 1000
  )
plot(dist_train_proj_2features, stat = "ecdf", type = "simple")
plot(dist_train_proj_2features, type = "simple")

dist_train_proj_unweighted <- geodist(
  x = training_points[sample.int(nrow(training_points), 1000),], # From docs: Predictor values for x, testdata and preddata are optional if modeldomain is a raster. If not provided they are extracted from the modeldomain rasterStack.
  modeldomain = preds_proj, # From docs: If type = "feature", the argument modeldomain...has to include predictors
  type = "feature",
  samplesize = 1000
)
plot(dist_train_proj_unweighted, stat = "ecdf", type = "simple")
plot(dist_train_proj_unweighted, type = "simple")

dist_train_proj_weighted <- geodist(
  x = training_points[sample.int(nrow(training_points), 1000),], # From docs: Predictor values for x, testdata and preddata are optional if modeldomain is a raster. If not provided they are extracted from the modeldomain rasterStack.
  modeldomain = preds_proj, # From docs: If type = "feature", the argument modeldomain...has to include predictors
  type = "feature",
  weight = weights_gini,
  samplesize = 1000
)
plot(dist_train_proj_weighted, stat = "ecdf", type = "simple")
plot(dist_train_proj_weighted, type = "simple")

## knndm() ####

### Debugging ####

# Error in check_knndm_feature(tpoints, predpoints, space, k, maxp, clustering,  : 
# Some values of factorar50are only present in training / prediction points.
# All factor values in the prediction points must be present in the training points.
# saveRDS(st_drop_geometry(training_points), "output/tpoints.rds")
# saveRDS(ppoints, "output/predpoints.rds")

### Unweighted knndm ####

?knndm
training_points_df <- training_points |> 
  st_drop_geometry() |> 
  as.data.frame()

set.seed(123)
ppoints <- preds_proj |> 
  as.data.frame(na.rm = TRUE) |>
  sample_n(1000) # Sample 1000 points, equivalent to default `samplesize` in `knndm()`

set.seed(123)
training_points_sample <- training_points |> 
  sample_n(1e2)  # Sample n points for shorter runtime
training_points_sample_df <- training_points_sample |> 
  st_drop_geometry() |> 
  as.data.frame()

tictoc::tic()
knndm_train_proj_unweighted <- knndm(
  tpoints = training_points_sample_df, # tpoints must be a `data.frame` (NOT tibble!) because of check_knndm_feature() with factor levels
  predpoints = ppoints,
  space = "feature",
  k = 10,
  samplesize = 1e3,
  weight = NULL 
)
tictoc::toc()

runtime <- tribble(
  ~ntpoints, ~elapsed,
  100, 6.8,
  1e3, 40.6,
  1e4, 1234
) 
plot(runtime)

knndm_train_proj_unweighted
plot(knndm_train_proj_unweighted, stat = "ecdf", type = "simple")
ggplot() + geom_sf(
  data = mutate(training_points_sample, 
                cvfold = as.factor(knndm_train_proj_unweighted$clusters)), 
  aes(color=cvfold))

# From knndm() docs: 
# If training data points are very clustered with respect to the prediction area
# and the presented `knndm` configuration still show signs of Gj * > Gij, there
# are several things that can be tried. First, increase the `maxp` parameter
# this may help to control for strong clustering (at the cost of having
# unbalanced folds). Secondly, decrease the number of final folds `k`, which may
# help to have larger clusters.

tictoc::tic()
knndm_train_proj_unweighted <- knndm(
  tpoints = training_points_sample_df, # tpoints must be a `data.frame` (NOT tibble!) because of check_knndm_feature() with factor levels
  predpoints = ppoints,
  space = "feature",
  k = 5, # Decrease k
  maxp = 0.5, # Increase maxp
  samplesize = 1e3,
  weight = NULL 
)
tictoc::toc()
plot(knndm_train_proj_unweighted, stat = "ecdf", type = "simple")

### Weighted knndm() ####

knndm_train_proj_weighted <- knndm(
  tpoints = training_points_sample_df, # tpoints must be a `data.frame` (NOT tibble!) because of check_knndm_feature() with factor levels
  predpoints = ppoints,
  space = "feature",
  k = 5, # Decrease k
  maxp = 0.5, # Increase maxp
  samplesize = 1e3,
  weight = weights_gini 
)
plot(knndm_train_proj_weighted, stat = "ecdf", type = "simple")

# sessionInfo ####

sessioninfo::session_info()

