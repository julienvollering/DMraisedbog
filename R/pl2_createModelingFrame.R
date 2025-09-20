# Create full modeling frame with all locations x scenarios #### 

library(tidyverse)
library(terra)
library(sf)

## Load data ####

### Current ####
rf_global_current <- rast("output/rf_global_pred_regional_current.tif")
names(rf_global_current) <- "rf_global"
preds_nor_250m_cur <- rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")
raster_current <- c(rf_global_current, preds_nor_250m_cur)

presence <- read_csv("output/presence_coords_regional.csv")
absence <- read_csv("output/absence_coords_regional.csv")

### Future ####
rf_global_future <- rast("output/rf_global_pred_regional_future.tif")
names(rf_global_future) <- "rf_global"
preds_nor_250m_fut <- rast("output/predictors_regional_250m_Norway_future_EPSG3035.tif")
raster_future <- c(rf_global_future, preds_nor_250m_fut)

## Define modeling population N ####

# Mask by artype_60 >= 1 | presence (i.e., peatland areas only)
mask_artype60 <- preds_nor_250m_cur[["artype_60"]] >= 0.5
mask_presence <- rasterize(presence, preds_nor_250m_cur, field = 1, background = 0)
mask <- mask_artype60 | (mask_presence == 1)
plot(mask)

writeRaster(mask, "output/pl2/mask_artype60PlusPresence.tif", overwrite=TRUE)

raster_current <- mask(raster_current, mask, maskvalue = 0, updatevalue = NA)
raster_future <- mask(raster_future, mask, maskvalue = 0, updatevalue = NA)

## Define modeling features p ####

# Exclude artype features 
raster_current <- raster_current[[!grepl("artype", names(raster_current))]]
raster_future <- raster_future[[!grepl("artype", names(raster_future))]]

writeRaster(raster_current, "output/pl2/scenario_current.tif", overwrite=TRUE)
writeRaster(raster_future, "output/pl2/scenario_future.tif", overwrite=TRUE)

## Create frames ####

### Current ####
# Extract at presence points
df_presence <- terra::extract(
  raster_current, presence[c("x", "y")],
  xy = TRUE,
  ID = FALSE) |>
  mutate(response = 1, .before = 1) |>
  drop_na() |>
  tibble()

# Extract at absence points
df_absence <- terra::extract(
  raster_current, absence[c("x", "y")],
  xy = TRUE,
  ID = FALSE) |>
  mutate(response = 0, .before = 1) |>
  drop_na() |>
  tibble()

# Combine presence and absence
current <- bind_rows(df_presence, df_absence) |>
  mutate(response = factor(response))
nrow(current)

### Future ####
# Extract at all non-NA locations in future raster
# Includes areas not surveyed by Lyngstad; more rows than current
future <- terra::as.data.frame(raster_future, xy = TRUE, na.rm = TRUE) |>
  mutate(response = NA, .before = 1) |>
  tibble()

### Combined
mf <- bind_rows(current = current, future = future, .id = "scenario") |> 
  mutate(scenario = factor(scenario, levels = c("current", "future")))
head(mf)
tail(mf)

mf |> 
  group_by(scenario, response) |> 
  count()
mf |> 
  group_by(scenario) |> 
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = FALSE))) |> 
  pivot_longer(-scenario) |>
  pivot_wider(names_from = scenario, values_from = value) |> 
  print(n = Inf)

## Save ####
write_csv(mf, "output/pl2/modeling_frame_regional.csv", append = FALSE)

# sessionInfo ####

sessioninfo::session_info()

