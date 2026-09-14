# Spatial Interpretation of Current and Future Predictions ####

# PURPOSE: Summarises the projection across Lyngstad polygons as class transitions and the change in P(bog).

# Reads the prediction rasters for current and future scenarios and summarizes them
# spatially across Lyngstad survey polygons (raised bogs).
#
# The 3-class response changes what "change" means here, and for the better. The binary
# version could only ask whether P(presence) crossed a prevalence threshold. The simplex
# says what a bog polygon is predicted to become: the interesting transition is
# bog -> otherpeat (the site stays peatland but stops being an active raised bog), and
# that was previously invisible -- both classes were the same `0`.
#
# Two readings are reported, because they answer different questions:
#  - the argmax class transition, the operating point implied by the balanced forest;
#  - the change in P(bog) itself, which is continuous and does not hide small shifts
#    behind a class boundary.

library(readr)
library(dplyr)
library(tidyr)
library(sf)
library(terra)

## Configuration ####

RESPONSE_LEVELS <- c("nonpeat", "otherpeat", "bog")

## Load spatial data ####

cat("Loading Lyngstad raised bog polygons...\n")
lyngstad <- st_read(
  "data/DMraisedbog.gpkg",
  layer = "lyngstad-MTYPE_A",
  quiet = TRUE
)
cat("Loaded", nrow(lyngstad), "polygons\n")

lyngstad_proj <- st_transform(lyngstad, crs = 3035)
lyngstad_vect <- vect(lyngstad_proj)

## Load predictions ####

pred_file <- "output/pl2/rf_local_pred_brf.tif"
pred_rasters <- rast(pred_file)
cat("Loaded", nlyr(pred_rasters), "prediction layers:", names(pred_rasters), "\n\n")

expected_layers <- c(
  paste0("current_", RESPONSE_LEVELS),
  paste0("future_", RESPONSE_LEVELS)
)
stopifnot(identical(names(pred_rasters), expected_layers))

## Extract zonal statistics ####

# Polygon means of each class probability. The three current layers still sum to 1 after
# averaging, so the polygon-level values remain a simplex and the argmax stays meaningful.
zonal <- terra::extract(
  pred_rasters,
  lyngstad_vect,
  fun = mean,
  na.rm = TRUE
) |>
  as_tibble()

# Polygons too small or entirely off the prediction grid return NA and cannot be
# classified either way; count them rather than letting them fall into a transition class.
n_missing <- sum(!complete.cases(zonal))
cat("Polygons with no prediction coverage:", n_missing, "of", nrow(zonal), "\n\n")

## Analyze change ####

current_mat <- as.matrix(zonal[, paste0("current_", RESPONSE_LEVELS)])
future_mat <- as.matrix(zonal[, paste0("future_", RESPONSE_LEVELS)])

argmax_class <- function(m) {
  out <- rep(NA_character_, nrow(m))
  ok <- complete.cases(m)
  out[ok] <- RESPONSE_LEVELS[max.col(m[ok, , drop = FALSE], ties.method = "first")]
  out
}

change <- zonal |>
  mutate(
    class_current = argmax_class(current_mat),
    class_future = argmax_class(future_mat),
    prob_bog_current = current_bog,
    prob_bog_future = future_bog,
    change_bog = future_bog - current_bog,
    transition = if_else(
      is.na(class_current) | is.na(class_future),
      NA_character_,
      paste(class_current, class_future, sep = " -> ")
    )
  ) |>
  select(
    ID,
    starts_with("current_"),
    starts_with("future_"),
    class_current,
    class_future,
    transition,
    prob_bog_current,
    prob_bog_future,
    change_bog
  )

## Create output polygon layer ####

lyngstad_predicted <- bind_cols(lyngstad_proj, select(change, -ID))

st_write(
  lyngstad_predicted,
  dsn = "output/pl2/lyngstad_predictions.gpkg",
  layer = "brf",
  delete_layer = TRUE,
  quiet = TRUE
)
cat("Saved polygon predictions to output/pl2/lyngstad_predictions.gpkg\n\n")

## Summary statistics ####

n_classified <- sum(!is.na(change$transition))

cat("Predicted class at Lyngstad raised bogs, current vs future:\n")
table(
  current = change$class_current,
  future = change$class_future,
  useNA = "ifany"
) |>
  print()

cat("\nTransitions (of", n_classified, "classified polygons):\n")
transitions <- change |>
  filter(!is.na(transition)) |>
  count(class_current, class_future, transition, sort = TRUE) |>
  mutate(prop = round(n / n_classified, 4))
print(as.data.frame(transitions), row.names = FALSE)
write_csv(transitions, "output/pl2/lyngstad_transitions.csv", append = FALSE)

cat("\nChange in P(bog):\n")
change |>
  summarise(
    n = sum(!is.na(change_bog)),
    mean_current = mean(prob_bog_current, na.rm = TRUE),
    mean_future = mean(prob_bog_future, na.rm = TRUE),
    mean_change = mean(change_bog, na.rm = TRUE),
    q05_change = quantile(change_bog, 0.05, na.rm = TRUE),
    median_change = median(change_bog, na.rm = TRUE),
    q95_change = quantile(change_bog, 0.95, na.rm = TRUE),
    prop_declining = mean(change_bog < 0, na.rm = TRUE)
  ) |>
  as.data.frame() |>
  print(row.names = FALSE)

hist(
  change$change_bog,
  breaks = 40,
  main = "Change in P(bog) at Lyngstad raised bogs",
  xlab = "P(bog) future - P(bog) current"
)
abline(v = 0, lty = 2)

# sessionInfo ####

sessioninfo::session_info()
