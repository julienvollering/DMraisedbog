library(tidyverse)
library(terra)
library(cluster)
library(sf)

## Load and prepare data ####

rf_global_current <- rast("output/rf_global_pred_regional_current.tif") 
names(rf_global_current) <- "rf_global"
chelsa_nor_250m_cur <- rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")[[1:19]]
ar50_nor_250m <- rast("output/ar50_250m_EPSG3035.tif")
names(ar50_nor_250m) <- "ar50"
ar50_nor_250m <- as.factor(ar50_nor_250m)
preds_train <- c(rf_global_current, chelsa_nor_250m_cur, ar50_nor_250m)

presence <- read_csv("output/presence_coords_regional.csv")
absence <- read_csv("output/absence_coords_regional.csv")

presence_df <- terra::extract(preds_train, presence[c("x", "y")], xy=TRUE, ID=FALSE) |> 
  mutate(response = 1, .before = 1) |>
  drop_na() |> 
  tibble()

absence_df <- terra::extract(preds_train, absence[c("x", "y")], xy=TRUE, ID=FALSE) |> 
  mutate(response = 0, .before = 1) |>
  drop_na() |> 
  tibble()

training_data_with_coords <- bind_rows(presence_df, absence_df) |>
  mutate(response = factor(response), ar50 = factor(ar50))

n_final_presence <- sum(training_data_with_coords$response == "1")
n_final_background <- sum(training_data_with_coords$response == "0")
n_final_prevalence <- n_final_presence / (n_final_presence + n_final_background)

cat("Final training data - Presences:", n_final_presence, "Background:", n_final_background, "\n")
cat("Prevalence in training data:", round(n_final_prevalence * 100, 3), "%\n")

weights_gini <- read_csv("output/feature_weights_gini.csv")
weights <- pull(weights_gini, weights_gini)
names(weights) <- pull(weights_gini, feature)

## Create clustering sample ####

sample_size <- 1e4
set.seed(42)
training_sample_with_coords <- training_data_with_coords |> 
  sample_n(sample_size)

cat("Sample prevalence:", 
    round(sum(training_sample_with_coords$response == "1") / nrow(training_sample_with_coords) * 100, 2), "%\n")

## Run clustering ####

# Select top 3 features by weight
top_features <- weights |> 
  sort(decreasing = TRUE) |> 
  head(3) |>
  names()

cat("Using top 3 features:", paste(top_features, collapse = ", "), "\n")

# Extract feature data for clustering (numeric only)
clustering_data <- training_sample_with_coords |>
  select(all_of(top_features)) |>
  select(where(is.numeric))  # Ensure only numeric features

cat("Clustering data dimensions:", nrow(clustering_data), "x", ncol(clustering_data), "\n")

start_time <- Sys.time()
# Run PAM clustering with Euclidean distance
set.seed(42)
pam_result <- pam(
  x = clustering_data,
  k = 5,  # Start with fewer clusters for simplicity
  metric = "euclidean",
  stand = TRUE,  # Standardize features
  pamonce = 2
)
end_time <- Sys.time()
cat("PAM clustering completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

# Add cluster assignments to training sample
training_sample_with_clusters <- training_sample_with_coords |>
  mutate(cluster = factor(pam_result$clustering))

# Show cluster summary
cluster_summary <- training_sample_with_clusters |>
  group_by(cluster) |>
  summarise(
    n_points = n(),
    n_presence = sum(as.numeric(as.character(response))),
    prevalence = n_presence / n_points,
    .groups = 'drop'
  )

cat("\nCluster summary:\n")
print(cluster_summary)

# Geographic visualization
training_sample_with_clusters |> 
  select(cluster, x, y) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(preds_train)) |>
  ggplot() +
  geom_sf(aes(color = cluster), size = 0.5) +
  labs(title = "PAM Clusters - Geographic Distribution")

## Assign new data to clusters ####

# Create new data sample
set.seed(123)
new_data_coords <- training_data_with_coords |> 
  sample_n(2000)  # Smaller sample for debugging

# Extract the same features for new data
new_data_features <- new_data_coords |>
  select(all_of(top_features)) |>
  select(where(is.numeric))

# Get medoids from PAM result
medoids <- pam_result$medoids |> 
  as_tibble()

cat("\nMedoid indices:", pam_result$id.med, "\n")
cat("Number of medoids:", nrow(medoids), "\n")

# Simple assignment function
assign_to_clusters <- function(new_data, medoids) {
  
  n_new <- nrow(new_data)
  assignments <- numeric(n_new)
  
  for (i in 1:n_new) {
    # Calculate Euclidean distances to each medoid
    distances <- numeric(nrow(medoids))
    
    for (j in 1:nrow(medoids)) {
      distances[j] <- sqrt(sum((new_data[i, ] - medoids[j, ])^2))
    }
    
    # Assign to closest medoid
    assignments[i] <- which.min(distances)
  }
  
  return(assignments)
}

# Assign clusters
cat("Assigning clusters to", nrow(new_data_features), "new points...\n")
new_assignments <- assign_to_clusters(new_data_features, medoids)

# Add assignments to new data
new_data_with_clusters <- new_data_coords |>
  mutate(cluster = factor(new_assignments))

# Compare cluster distributions
cat("\nCluster distribution comparison:\n")
cat("Training sample:\n")
print(table(training_sample_with_clusters$cluster))
cat("\nNew data:\n")
print(table(new_data_with_clusters$cluster))

## Diagnostic visualization ####

# Geographic plot of new assignments
new_data_with_clusters |>
  select(cluster, x, y) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(preds_train)) |>
  ggplot() +
  geom_sf(aes(color = cluster), size = 0.5) +
  labs(title = "New Data - Cluster Assignments")

# Feature space visualization (if 2D or 3D)
if (ncol(clustering_data) >= 2) {
  # Plot first two features
  feature_plot_data <- bind_rows(
    training_sample_with_clusters |>
      select(all_of(names(clustering_data)[1:2]), cluster) |>
      mutate(data_type = "Training"),
    new_data_with_clusters |>
      select(all_of(names(clustering_data)[1:2]), cluster) |>
      mutate(data_type = "New")
  )
  
  # Add medoids to plot
  medoid_plot_data <- medoids[1:2] |>
    mutate(cluster = factor(1:nrow(medoids)), data_type = "Medoid")
  
  ggplot() +
    geom_point(data = feature_plot_data, 
               aes_string(x = names(clustering_data)[1], 
                         y = names(clustering_data)[2], 
                         color = "cluster", 
                         shape = "data_type"), 
               alpha = 0.6, size = 1) +
    geom_point(data = medoid_plot_data,
               aes_string(x = names(clustering_data)[1], 
                         y = names(clustering_data)[2], 
                         color = "cluster"),
               shape = "X", size = 4, stroke = 2) +
    labs(title = "Feature Space: Training vs New Data",
         subtitle = "X marks show medoid locations") +
    theme_minimal()
}

# Detailed diagnostic for first few points
cat("\n=== ASSIGNMENT DIAGNOSTICS ===\n")
for (i in 1:5) {
  cat("\nTest point", i, ":\n")
  
  # Show feature values
  point_features <- new_data_features[i, ]
  cat("Features:", paste(names(point_features), "=", round(as.numeric(point_features), 3), collapse = ", "), "\n")
  
  # Calculate distances to all medoids
  distances <- numeric(nrow(medoids))
  for (j in 1:nrow(medoids)) {
    distances[j] <- sqrt(sum((point_features - medoids[j, ])^2))
  }
  
  cat("Distances to medoids:", paste(round(distances, 4), collapse = ", "), "\n")
  cat("Assigned to cluster:", which.min(distances), "\n")
  cat("Minimum distance:", round(min(distances), 4), "\n")
}

# sessionInfo ####

sessioninfo::session_info()
