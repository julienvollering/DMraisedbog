# Compare feature space distances between observations and modeling domain

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

## Read modeling frame and domain raster ####
mf <- read_csv("output/pl2/modeling_frame_regional.csv")

raster_future <- rast("output/pl2/scenario_future.tif")

## Load feature weights ####
weights_features <- read_csv("output/pl2/weights_feature_data_partitioning.csv")
weights_features <- weights_features |>
  mutate(feature = forcats::fct_reorder(
    feature, median_normalized,
    .fun = first
  ))

weights <- weights_features |>
  filter(grepl("Balanced", method)) |>
  select(feature, median_normalized) |>
  arrange(desc(median_normalized)) |>
  tibble::deframe()

## Distances ####

mf_current <- mf |> 
  filter(scenario == "current") |>
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

### With domain = all future ####

#### All features, unweighted ####
dists <- CAST:::geodist.fast(
  sample = mf_current,
  prediction = mf_future,
  weights = NULL,
  samplesize = 1e4) # Increase samplesize for stability 
plot(dists, type = "simple")

#### All features, weighted ####
dists <- CAST:::geodist.fast(
  sample = mf_current,
  prediction = mf_future,
  weights = weights,
  samplesize = 1e4) # Increase samplesize for stability 
plot(dists, type = "simple")

#### 5 features, unweighted ####
dists <- CAST:::geodist.fast(
  sample = mf_current[,names(weights)[1:5]],
  prediction = mf_future[,names(weights)[1:5]],
  weights = NULL,
  samplesize = 1e4) # Increase samplesize for stability 
plot(dists, type = "simple")

#### 5 features, weighted ####
dists <- CAST:::geodist.fast(
  sample = mf_current[,names(weights)[1:5]],
  prediction = mf_future[,names(weights)[1:5]],
  weights = weights,
  samplesize = 1e4) # Increase samplesize for stability 
plot(dists, type = "simple")

### With domain = presences x future ####

#### All features, unweighted ####
dists <- CAST:::geodist.fast(
  sample = mf_current,
  prediction = mf_future_presence,
  weights = NULL,
  samplesize = 1e4) # Increase samplesize for stability 
plot(dists, type = "simple")

#### All features, weighted ####
dists <- CAST:::geodist.fast(
  sample = mf_current,
  prediction = mf_future_presence,
  weights = weights,
  samplesize = 1e4) # Increase samplesize for stability 
plot(dists, type = "simple")

#### 5 features, unweighted ####
dists <- CAST:::geodist.fast(
  sample = mf_current[,names(weights)[1:5]],
  prediction = mf_future_presence[,names(weights)[1:5]],
  weights = NULL,
  samplesize = 1e4) # Increase samplesize for stability 
plot(dists, type = "simple")

#### 5 features, weighted ####
dists <- CAST:::geodist.fast(
  sample = mf_current[,names(weights)[1:5]],
  prediction = mf_future_presence[,names(weights)[1:5]],
  weights = weights,
  samplesize = 1e4) # Increase samplesize for stability 
plot(dists, type = "simple")

#### Test when sample == prediction ####

mf_current_presence <- mf |> 
  filter(scenario == "current", response == 1) |>
  select(-scenario, -response) |> 
  st_as_sf(coords = c("x", "y"), crs = 3035)

dists <- CAST:::geodist.fast(
  sample = mf_current_presence[,"elevation"],
  prediction = mf_future_presence[,"elevation"],
  weights = NULL,
  samplesize = 1e2) # Increase samplesize for stability
plot(dists, type = "simple")
dists |> 
  filter(what == "prediction-to-sample") |>
  pull(dist) |> 
  unique()
# Sample-to-sample do not include NN to self, while prediction-to-sample do

# sessionInfo ####

sessioninfo::session_info()
