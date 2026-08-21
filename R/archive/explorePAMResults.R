library(tidyverse)
library(terra)
library(fpc)
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

training_data_with_coords |> 
  select(response, x, y) |> 
  st_as_sf(coords = c("x", "y"), crs = st_crs(preds_train)) |>
  sample_n(1e5) |> 
  arrange(response) |>
  ggplot() +
  geom_sf(aes(color = response), size = 0.5)

n_final_presence <- sum(training_data_with_coords$response == "1")
n_final_background <- sum(training_data_with_coords$response == "0")
n_final_prevalence <- n_final_presence / (n_final_presence + n_final_background)

cat("Final training data - Presences:", n_final_presence, "Background:", n_final_background, "\n")
cat("Prevalence in training data:", round(n_final_prevalence * 100, 3), "%\n")

weights_gini <- read_csv("output/feature_weights_gini.csv")
weights <- pull(weights_gini, weights_gini)
names(weights) <- pull(weights_gini, feature)

## Create clustering sample ####

# Use random sample instead of stratified to see if this helps cluster 1 issue
sample_size <- 1e4 # 1e5 produces error in call to cluster::daisy()
set.seed(42)
training_sample_with_coords <- training_data_with_coords |> 
  sample_n(sample_size)

sum(training_sample_with_coords$response == "1") / 
  (sum(training_sample_with_coords$response == "1") + 
     sum(training_sample_with_coords$response == "0"))

training_sample_with_coords |>
  select(response, x, y) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(preds_train)) |>
  arrange(response) |>
  ggplot() +
  geom_sf(aes(color = response), size = 0.5)

training_sample <- training_sample_with_coords |>
  select(-x, -y)

## Run clustering ####

### With unweighted distance matrix ####

# feature_cols <- setdiff(names(training_sample), "response")
# 
# gower_dist_unweighted <- cluster::daisy(
#   x = training_sample[, feature_cols],
#   metric = "gower",
#   weights = 1 # alt: weights
# )
# 
# set.seed(42)
# pam_result <- cluster::pam(
#   x = gower_dist_unweighted,
#   k = 15,  # Number of clusters
#   diss = TRUE,
#   cluster.only = TRUE,
#   variant = "faster"
# )
# 
# training_sample_with_clusters <- training_sample_with_coords |>
#   mutate(cluster = factor(pam_result))
# 
# training_sample_with_clusters |>
#   group_by(cluster) |>
#   summarise(
#     n_points = n(),
#     n_presence = sum(as.numeric(as.character(response))),
#     prevalence = n_presence / n_points,
#     .groups = 'drop'
#   ) |>
#   arrange(cluster) |> 
#   print(n=50)
# 
# # PCA for feature space
# numeric_features <- training_sample |> 
#   select(where(is.numeric)) |>
#   na.omit()
# 
# pca_result <- prcomp(numeric_features, scale. = TRUE)
# pca_data <- data.frame(
#   PC1 = pca_result$x[,1],
#   PC2 = pca_result$x[,2],
#   cluster = factor(pam_result),
#   response = training_sample$response[1:nrow(numeric_features)]
# )
# plot(pca_data$PC1, pca_data$PC2, col = as.numeric(pca_data$cluster), 
#      main = "Unweighted PAM in PCA Space", xlab = "PC1", ylab = "PC2")
# 
# # Geographic plot using sf
# training_sample_with_clusters |> 
#   select(cluster, x, y) |>
#   st_as_sf(coords = c("x", "y"), crs = st_crs(preds_train)) |>
#   arrange(cluster) |>
#   ggplot() +
#   geom_sf(aes(color = cluster), size = 0.5) +
#   labs(title = "Unweighted PAM")

### With weighted distance matrix ####

feature_cols <- setdiff(names(training_sample), "response")
weights |> 
  as_tibble(rownames = "feature") |> 
  mutate(feature = fct_reorder(feature, desc(value))) |>
  ggplot() + geom_col(aes(x = feature, y = value))

# Add progress monitoring for clustering
cat("\n=== CLUSTERING PROGRESS ===\n")
cat("Computing Gower distance matrix...\n")
cat("Dataset size:", nrow(training_sample), "x", length(feature_cols), "features\n")
cat("Start time:", format(Sys.time()), "\n")

gower_start <- Sys.time()
gower_dist_weighted <- cluster::daisy(
  x = training_sample[, feature_cols],
  metric = "gower",
  weights = weights
)
gower_end <- Sys.time()
cat("Gower distance completed in", round(difftime(gower_end, gower_start, units = "mins"), 2), "minutes\n")

cat("\nStarting PAM clustering (k=15)...\n")
cat("This may take several minutes for", nrow(training_sample), "points...\n")
pam_start <- Sys.time()

set.seed(42)
pam_result <- cluster::pam(
  x = gower_dist_weighted,
  k = 15,  # Number of clusters
  diss = TRUE,
  cluster.only = FALSE,
  nstart = 10,
  pamonce = 1
)

pam_end <- Sys.time()
total_time <- difftime(pam_end, gower_start, units = "mins")
pam_time <- difftime(pam_end, pam_start, units = "mins")

cat("PAM clustering completed!\n")
cat("PAM time:", round(pam_time, 2), "minutes\n")  
cat("Total clustering time:", round(total_time, 2), "minutes\n") # 37 min
cat("=== CLUSTERING COMPLETE ===\n\n") 

saveRDS(pam_result, "output/pam_result_weighted.rds")
# pam_result <- readRDS("output/pam_result_weighted.rds")

training_sample_with_clusters <- training_sample_with_coords |>
  mutate(cluster = factor(pam_result$clustering))

training_sample_with_clusters |>
  group_by(cluster) |>
  summarise(
    n_points = n(),
    n_presence = sum(as.numeric(as.character(response))),
    prevalence = n_presence / n_points,
    .groups = 'drop'
  ) |>
  arrange(cluster) |> 
  print(n=50)

# PCA for feature space
numeric_features <- training_sample |> 
  select(where(is.numeric)) |>
  na.omit()

pca_result <- prcomp(numeric_features, scale. = TRUE)
pca_data <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  cluster = factor(pam_result$clustering),
  response = training_sample$response[1:nrow(numeric_features)]
)
plot(pca_data$PC1, pca_data$PC2, col = as.numeric(pca_data$cluster), 
     main = "Weighted PAM in PCA Space", xlab = "PC1", ylab = "PC2")

# Geographic plot using sf
training_sample_with_clusters |> 
  select(cluster, x, y) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(preds_train)) |>
  arrange(cluster) |>
  ggplot() +
  geom_sf(aes(color = cluster), size = 0.5) +
  labs(title = "Weighted PAM")

### With selected features ####

feature_cols <- setdiff(names(training_sample), "response")
weights |> 
  as_tibble(rownames = "feature") |> 
  mutate(feature = fct_reorder(feature, desc(value))) |>
  arrange(desc(value))
selected_features <- c("bio10", "rf_global", "bio1")

# Add progress monitoring for clustering
cat("\n=== CLUSTERING PROGRESS ===\n")
cat("Computing Gower distance matrix...\n")
cat("Dataset size:", nrow(training_sample), "x", length(selected_features), "features\n")
cat("Start time:", format(Sys.time()), "\n")

gower_start <- Sys.time()
gower_dist_weighted <- cluster::daisy(
  x = training_sample[, selected_features],
  metric = "gower",
  weights = weights
)
gower_end <- Sys.time()
cat("Gower distance completed in", round(difftime(gower_end, gower_start, units = "mins"), 2), "minutes\n")

cat("\nStarting PAM clustering (k=15)...\n")
cat("This may take several minutes for", nrow(training_sample), "points...\n")
pam_start <- Sys.time()

set.seed(42)
pam_result <- cluster::pam(
  x = gower_dist_weighted,
  k = 15,  # Number of clusters
  diss = TRUE,
  cluster.only = FALSE,
  nstart = 10,
  pamonce = 6
)

pam_end <- Sys.time()
total_time <- difftime(pam_end, gower_start, units = "mins")
pam_time <- difftime(pam_end, pam_start, units = "mins")

cat("PAM clustering completed!\n")
cat("PAM time:", round(pam_time, 2), "minutes\n")  
cat("Total clustering time:", round(total_time, 2), "minutes\n") # 37 min
cat("=== CLUSTERING COMPLETE ===\n\n") 

training_sample_with_clusters <- training_sample_with_coords |>
  mutate(cluster = factor(pam_result$clustering))

training_sample_with_clusters |>
  group_by(cluster) |>
  summarise(
    n_points = n(),
    n_presence = sum(as.numeric(as.character(response))),
    prevalence = n_presence / n_points,
    .groups = 'drop'
  ) |>
  arrange(cluster) |> 
  print(n=50)

# PCA for feature space
numeric_features <- training_sample |> 
  select(where(is.numeric)) |>
  na.omit()

pca_result <- prcomp(numeric_features, scale. = TRUE)
pca_data <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  cluster = factor(pam_result$clustering),
  response = training_sample$response[1:nrow(numeric_features)]
)
plot(pca_data$PC1, pca_data$PC2, col = as.numeric(pca_data$cluster), 
     main = "Weighted PAM in PCA Space", xlab = "PC1", ylab = "PC2")

# Geographic plot using sf
training_sample_with_clusters |> 
  select(cluster, x, y) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(preds_train)) |>
  arrange(cluster) |>
  ggplot() +
  geom_sf(aes(color = cluster), size = 0.5) +
  labs(title = "Weighted PAM")

## Assign new data to established clusters ####

# Create new data subset for cluster assignment (excluding the training sample)
set.seed(123)
new_data <- training_data_with_coords |> 
  sample_n(1e4) # Sample for testing runtime

# Batch assignment function with progress bar
assign_clusters_batch <- function(
    new_data, 
    medoids_data, 
    feature_cols, 
    weights, 
    batch_size = 500) {
  library(progress)
  
  n_new <- nrow(new_data)
  n_batches <- ceiling(n_new / batch_size)
  
  # Initialize progress bar
  pb <- progress_bar$new(
    format = "  Assigning clusters [:bar] :percent in :elapsed",
    total = n_batches, clear = FALSE, width = 60
  )
  
  # Initialize results vector
  cluster_assignments <- numeric(n_new)
  
  # Process in batches
  for (i in 1:n_batches) {
    pb$tick()
    
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, n_new)
    batch_indices <- start_idx:end_idx
    
    # Get batch data
    batch_data <- new_data[batch_indices, feature_cols, drop = FALSE]
    
    # Calculate distances from batch points to each medoid
    batch_assignments <- numeric(length(batch_indices))
    
    for (j in seq_along(batch_indices)) {
      point_data <- batch_data[j, , drop = FALSE]
      
      # Calculate Gower distance to each medoid
      distances_to_medoids <- numeric(nrow(medoids_data))
      
      for (k in 1:nrow(medoids_data)) {
        medoid_data <- medoids_data[k, , drop = FALSE]
        combined_data <- rbind(point_data, medoid_data)
        
        # Calculate weighted Gower distance
        dist_matrix <- cluster::daisy(
          x = combined_data,
          metric = "gower",
          weights = weights
        )
        distances_to_medoids[k] <- as.matrix(dist_matrix)[1, 2]
      }
      
      # Assign to closest medoid
      batch_assignments[j] <- which.min(distances_to_medoids)
    }
    
    cluster_assignments[batch_indices] <- batch_assignments
  }
  
  pb$terminate()
  return(cluster_assignments)
}

# Extract medoid data from the PAM result
medoid_indices <- pam_result$medoids
medoids_data <- training_sample[medoid_indices, selected_features, drop = FALSE]

cat("Starting batch cluster assignment with", nrow(medoids_data), "medoids...\n")
start_time <- Sys.time()

# Assign clusters to new data
new_cluster_assignments <- assign_clusters_batch(
  new_data = new_data,
  medoids_data = medoids_data,
  feature_cols = selected_features,
  weights = weights[selected_features],
  batch_size = 500
)

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat("Cluster assignment completed in", round(runtime, 2), "seconds\n")
cat("Average time per point:", round(runtime / nrow(new_data) * 1000, 3), "milliseconds\n")

# Estimate time for full dataset
full_data_size <- nrow(training_data_with_coords) - nrow(training_sample_with_coords)
estimated_full_runtime <- (runtime / nrow(new_data)) * full_data_size
cat("Estimated runtime for full dataset (", full_data_size, "points):", 
    round(estimated_full_runtime / 60, 1), "minutes\n")

# Add cluster assignments to new data
new_data_with_clusters <- new_data |>
  mutate(cluster = factor(new_cluster_assignments))

# Summary of cluster assignments
cat("\nCluster assignment summary:\n")
new_data_with_clusters |>
  count(cluster, sort = TRUE) |>
  print()

# Quick validation plot
new_data_with_clusters |>
  select(cluster, x, y) |>
  sample_n(10000) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(preds_train)) |>
  ggplot() +
  geom_sf(aes(color = cluster), size = 0.5) +
  labs(title = "New Data Cluster Assignments")

## DIAGNOSTICS: Investigating cluster assignment issues ####

cat("\n=== DIAGNOSTIC ANALYSIS ===\n")

# 1. Check medoid extraction and characteristics
cat("\n1. MEDOID ANALYSIS:\n")
cat("PAM medoid indices:", medoid_indices, "\n")
cat("Number of medoids extracted:", nrow(medoids_data), "\n")
cat("Expected number of clusters:", length(unique(pam_result$clustering)), "\n")

# Compare medoids to cluster centers
medoid_cluster_mapping <- data.frame(
  medoid_idx = medoid_indices,
  cluster_id = 1:length(medoid_indices),
  original_cluster = pam_result$clustering[medoid_indices]
)
print("Medoid to cluster mapping:")
print(medoid_cluster_mapping)

# Check if medoids actually represent their assigned clusters
cat("\n2. MEDOID VALIDATION:\n")
for(i in 1:length(medoid_indices)) {
  original_cluster <- pam_result$clustering[medoid_indices[i]]
  cat("Medoid", i, "should represent cluster", original_cluster, 
      "but is assigned as cluster", i, "\n")
}

# 3. Distance distribution diagnostics
cat("\n3. DISTANCE DISTRIBUTION ANALYSIS:\n")

# Sample a subset of new data for detailed analysis
set.seed(456)
diagnostic_sample <- new_data_with_clusters |> sample_n(min(1000, nrow(new_data_with_clusters)))

# Calculate distances from diagnostic sample to all medoids
distance_diagnostics <- list()

for(i in 1:nrow(diagnostic_sample)) {
  point_data <- diagnostic_sample[i, selected_features, drop = FALSE]
  distances <- numeric(nrow(medoids_data))
  
  for(k in 1:nrow(medoids_data)) {
    medoid_data <- medoids_data[k, , drop = FALSE]
    combined_data <- rbind(point_data, medoid_data)
    
    dist_matrix <- cluster::daisy(
      x = combined_data,
      metric = "gower", 
      weights = weights[selected_features]
    )
    distances[k] <- as.matrix(dist_matrix)[1, 2]
  }
  
  distance_diagnostics[[i]] <- data.frame(
    point_id = i,
    assigned_cluster = diagnostic_sample$cluster[i],
    medoid_1_dist = distances[1],
    medoid_2_dist = distances[2], 
    medoid_3_dist = distances[3],
    min_dist = min(distances),
    min_dist_cluster = which.min(distances),
    second_min_dist = sort(distances)[2],
    dist_ratio = sort(distances)[2] / min(distances)
  )
}

distance_df <- do.call(rbind, distance_diagnostics)

# Check for assignment mismatches
assignment_check <- distance_df |>
  mutate(
    correct_assignment = as.numeric(assigned_cluster) == min_dist_cluster,
    assignment_error = !correct_assignment
  )

cat("Assignment accuracy check:\n")
cat("Correct assignments:", sum(assignment_check$correct_assignment), "out of", nrow(assignment_check), "\n")
cat("Assignment errors:", sum(assignment_check$assignment_error), "\n")
cat("Assignment accuracy:", round(mean(assignment_check$correct_assignment) * 100, 2), "%\n")

if(sum(assignment_check$assignment_error) > 0) {
  cat("\nERRORS DETECTED! Showing examples of misassigned points:\n")
  error_examples <- assignment_check |> 
    filter(assignment_error) |>
    head(10)
  print(error_examples)
}

# 4. Compare training vs new data cluster distributions
cat("\n4. CLUSTER DISTRIBUTION COMPARISON:\n")

training_cluster_dist <- training_sample_with_clusters |>
  count(cluster, name = "training_count") |>
  mutate(training_prop = training_count / sum(training_count))

new_data_cluster_dist <- new_data_with_clusters |>
  count(cluster, name = "new_count") |>
  mutate(new_prop = new_count / sum(new_count))

cluster_comparison <- training_cluster_dist |>
  full_join(new_data_cluster_dist, by = "cluster") |>
  replace_na(list(training_count = 0, new_count = 0, training_prop = 0, new_prop = 0)) |>
  mutate(
    prop_diff = new_prop - training_prop,
    fold_change = ifelse(training_prop > 0, new_prop / training_prop, Inf)
  )

print("Cluster distribution comparison (training vs new):")
print(cluster_comparison)

# 5. Feature space diagnostics 
cat("\n5. FEATURE SPACE ANALYSIS:\n")

# Compare feature distributions between training clusters and new assignments
feature_comparison <- list()

for(cluster_id in unique(training_sample_with_clusters$cluster)) {
  
  # Training data for this cluster
  training_cluster_data <- training_sample_with_clusters |>
    filter(cluster == cluster_id) |>
    select(all_of(feature_cols))
  
  # New data assigned to this cluster  
  new_cluster_data <- new_data_with_clusters |>
    filter(cluster == cluster_id) |>
    select(all_of(feature_cols))
  
  if(nrow(new_cluster_data) > 0 && nrow(training_cluster_data) > 0) {
    
    # Compare means for numeric features
    numeric_cols <- names(select(training_cluster_data, where(is.numeric)))
    
    if(length(numeric_cols) > 0) {
      training_means <- training_cluster_data |>
        select(all_of(numeric_cols)) |>
        summarise(across(everything(), mean, na.rm = TRUE)) |>
        pivot_longer(everything(), names_to = "feature", values_to = "training_mean")
      
      new_means <- new_cluster_data |>
        select(all_of(numeric_cols)) |>
        summarise(across(everything(), mean, na.rm = TRUE)) |>
        pivot_longer(everything(), names_to = "feature", values_to = "new_mean")
      
      feature_comparison[[as.character(cluster_id)]] <- training_means |>
        left_join(new_means, by = "feature") |>
        mutate(
          cluster = cluster_id,
          mean_diff = new_mean - training_mean,
          relative_diff = ifelse(training_mean != 0, mean_diff / abs(training_mean), NA)
        )
    }
  }
}

if(length(feature_comparison) > 0) {
  feature_comparison_df <- do.call(rbind, feature_comparison)
  
  cat("Feature mean differences (new - training) by cluster:\n")
  feature_comparison_df |>
    filter(abs(relative_diff) > 0.1, !is.na(relative_diff)) |>  # Show substantial differences
    arrange(cluster, desc(abs(relative_diff))) |>
    print(n = 50)
}

# 6. Visualization of potential issues
cat("\n6. CREATING DIAGNOSTIC PLOTS:\n")

# Plot distance ratios (lower ratios = more ambiguous assignments)
distance_df |>
  ggplot(aes(x = dist_ratio)) +
  geom_histogram(bins = 30, alpha = 0.7) +
  geom_vline(xintercept = 1.5, color = "red", linetype = "dashed") +
  labs(
    title = "Distance Ratio Distribution", 
    subtitle = "Ratio of 2nd closest / closest medoid distance",
    x = "Distance Ratio", 
    y = "Count"
  ) +
  annotate("text", x = 2, y = Inf, vjust = 2, 
           label = "Higher ratios = more confident assignments")

# Plot assignment errors by cluster
if(sum(assignment_check$assignment_error) > 0) {
  assignment_check |>
    ggplot(aes(x = assigned_cluster, fill = assignment_error)) +
    geom_bar(position = "fill") +
    labs(
      title = "Assignment Error Rate by Cluster",
      x = "Assigned Cluster", 
      y = "Proportion",
      fill = "Assignment Error"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Deep diagnostic - let's see EXACTLY what's going on
debug_assignments <- function(new_data, train_data, pam_model) {
  
  medoids <- train_data[pam_model$id.med, ]
  
  # Take just 5 new points for detailed inspection
  test_points <- new_data |> 
    sample_n(5)
  
  for (i in 1:nrow(test_points)) {
    cat("\n=== Test Point", i, "===\n")
    
    # Calculate distances to ALL medoids
    combined <- rbind(test_points[i, , drop = FALSE], medoids)
    dist_matrix <- as.matrix(cluster::daisy(combined, metric = "gower"))
    dists_to_medoids <- dist_matrix[1, -1]
    
    # Show all distances
    cat("Distances to medoids:\n")
    for (j in 1:length(dists_to_medoids)) {
      cat(sprintf("  Medoid %d: %.10f\n", j, dists_to_medoids[j]))
    }
    
    # Which is minimum?
    min_idx <- which.min(dists_to_medoids)
    cat("which.min() returns:", min_idx, "\n")
    
    # Are there near-ties?
    sorted_dists <- sort(dists_to_medoids)
    cat("Sorted distances:", sorted_dists[1:3], "\n")
    cat("Difference between 1st and 2nd:", sorted_dists[2] - sorted_dists[1], "\n")
  }
}

# Run this first
debug_assignments(new_data[, selected_features], 
                  training_sample[, selected_features], 
                  pam_result)
