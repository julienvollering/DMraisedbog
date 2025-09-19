# Explore Quantile Random Forests with Imbalanced Data
# Short script for spatial predictions and G-mean evaluation

library(tidyverse)
library(sf)
library(terra)
library(randomForestSRC)

# Load partitioned data
data <- read_csv("output/pl1/data_partitioned.csv",
                          col_types = cols(
                            response = col_factor(),
                            ar50 = col_factor(),
                            outer = col_character(),
                            inner = col_integer(),
                            x = col_double(),
                            y = col_double()
                          ))
str(data)

## Unpartitioned ####

# Subset to small sample for fast training (adjust n as needed)
set.seed(123)
sample_data <- data |>
  slice_sample(n = 1e5) |>
  select(-outer, -inner) |>  # Remove partition columns
  select(-x, -y) |>  # Remove coordinates for model training
  as.data.frame()
str(sample_data)
  
# Prepare response as factor
sum(sample_data$response == 1)  # Check presence count

# Quick QRF model with imbalanced setting
qrf_model <- imbalanced(
  response ~ .,
  data = sample_data,
  ntree = 3000, 
  mtry = 3,
  nodesize = 10,
  rfq = TRUE
)

qrf_model
get.imbalanced.performance(qrf_model) |> 
  enframe()

# Make predictions on training sample
predictions <- predict(qrf_model, sample_data)

# Calculate G-mean (reproduce tuning CV calculation)
conf_matrix <- table(sample_data$response, predictions$class)
sensitivity <- conf_matrix[2,2] / sum(conf_matrix[2,])  # TPR
specificity <- conf_matrix[1,1] / sum(conf_matrix[1,])  # TNR
gmean <- sqrt(sensitivity * specificity)

cat("Model Performance:\n")
cat("Sensitivity (TPR):", round(sensitivity, 3), "\n")
cat("Specificity (TNR):", round(specificity, 3), "\n")
cat("G-mean:", round(gmean, 3), "\n")

### Create spatial dataset ####
rf_global_current <- rast("output/rf_global_pred_regional_current.tif") 
names(rf_global_current) <- "rf_global"
preds_nor_250m_cur <- rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")
predictors <- c(rf_global_current, preds_nor_250m_cur)
names(predictors) %in% names(sample_data)  # Check names match
names(sample_data) %in% names(predictors)  # Check names match

plot(predictors$bio1)
predictors_se <- crop(predictors, ext(435e4, 447e4, 410e4, 420e4))
plot(predictors_se$bio1)
  
# Make spatial predictions
predfun <- function(model, data) {
  v <- predict.rfsrc(model, data, importance = FALSE)
  cbind(p = v$predicted[,1], class = v$class) # Taking the first of two classes to be the positive class
}
spatial_preds <- terra::predict(predictors_se, qrf_model, 
                                fun = predfun, na.rm = TRUE)
plot(spatial_preds$p)
plot(spatial_preds$class)
prevalence_training <- mean(sample_data$response == 1)
plot(spatial_preds$p < prevalence_training)

## Partitioned ####

data_training_fold12 <- data |> 
  filter(outer == "train", inner == c(1, 2)) |> 
  slice_sample(n = 1e5) |>  # Subsample for fast training
  select(-outer, -inner) |>  # Remove partition columns
  select(-x, -y) |>  # Remove coordinates for model training
  as.data.frame()

data_training_fold3 <- data |> 
  filter(outer == "train", inner == 3) |> 
  slice_sample(n = 1e5) |>  # Subsample for fast evaluation
  select(-outer, -inner) |>  # Remove partition columns
  select(-x, -y) |>  # Remove coordinates for model training
  as.data.frame()

# Prepare response as factor
sum(data_training_fold12$response == 1)  # Check presence count
mean(data_training_fold12$response == 1)  # Check prevalence

# Quick QRF model with imbalanced setting
qrf_model <- imbalanced(
  response ~ .,
  data = data_training_fold12,
  ntree = 3000, 
  mtry = 3,
  nodesize = 10,
  rfq = TRUE
)

qrf_model
get.imbalanced.performance(qrf_model) |> 
  enframe()

# Make predictions on leave out fold
predictions <- predict(qrf_model, data_training_fold3)

# Calculate G-mean (reproduce tuning CV calculation)
conf_matrix <- table(data_training_fold3$response, predictions$class)
sensitivity <- conf_matrix[2,2] / sum(conf_matrix[2,])  # TPR
specificity <- conf_matrix[1,1] / sum(conf_matrix[1,])  # TNR
gmean <- sqrt(sensitivity * specificity)

cat("Model Performance:\n")
cat("Sensitivity (TPR):", round(sensitivity, 3), "\n")
cat("Specificity (TNR):", round(specificity, 3), "\n")
cat("G-mean:", round(gmean, 3), "\n")

# sessionInfo

sessioninfo::session_info()
