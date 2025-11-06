library(tidyverse)
library(sf)

# Read data ####
# Imagery listed in reports
lyngstad1 <- read_csv("data/NiB/lyngstad-tabell3.csv", col_types = 'ccdcc')
# Imagery not listed in reports but included in the dataset
lyngstad2 <- read_csv(
  "data/NiB/lyngstad-tabell3mangler.csv",
  col_types = 'ccdcc'
)
lyngstad <- bind_rows(lyngstad1, lyngstad2)

gpkgfiles <- list.files(
  "data/NiB/",
  pattern = '.gpkg',
  full.names = TRUE,
  recursive = TRUE
)
polygons <- gpkgfiles |>
  map(st_read) |>
  bind_rows()

# Clean and standardize the stripe numbers in reference and polygons datasets ####
# We'll preserve the original values and create standardized versions for matching

# First, clean the polygons dataset
polygons_clean <- polygons %>%
  mutate(
    # For regular numeric stripes, convert to character
    stripenummer_clean = case_when(
      str_detect(stripenummer, "^\\d+$") ~ as.character(as.numeric(
        stripenummer
      )), # Normalize numeric values
      TRUE ~ as.character(stripenummer) # Keep other formats as is
    )
  )

# Warning message:
#   There was 1 warning in `stopifnot()`.
# ℹ In argument: `stripenummer_clean = case_when(...)`.
# Caused by warning:
#   ! NAs introduced by coercion

# filter(polygons_clean, is.na(stripenummer_clean))
# Ignore warning

# Next, clean the lyngstad dataset
lyngstad_clean <- lyngstad %>%
  mutate(
    # Standardize stripe numbers similarly to polygons
    stripe_clean = case_when(
      str_detect(stripe, "^\\d+$") ~ as.character(as.numeric(stripe)), # Normalize numeric values
      TRUE ~ as.character(stripe) # Keep other formats as is
    ),
    # Check for ranges in number column
    is_range = str_detect(number, "-")
  )

# Warning message:
#   There was 1 warning in `stopifnot()`.
# ℹ In argument: `stripe_clean = case_when(...)`.
# Caused by warning:
#   ! NAs introduced by coercion

# filter(lyngstad_clean, is.na(stripe_clean))
# Ignore warning

# Function to expand ranges
expand_number_ranges <- function(df) {
  # Extract special "alle" cases - these will be handled separately
  alle_cases <- df %>%
    filter(number == "alle") %>%
    mutate(is_alle = TRUE, bildenummer = "match_all") # Placeholder value

  # For non-ranges (and not "alle"), keep as is
  non_ranges <- df %>%
    filter(!is_range & number != "alle") %>%
    mutate(bildenummer = number, is_alle = FALSE)

  # For ranges, expand into individual numbers
  ranges <- df %>%
    filter(is_range & number != "alle") %>%
    mutate(
      start_num = as.numeric(str_extract(number, "^\\d+")),
      end_num = as.numeric(str_extract(number, "\\d+$"))
    ) %>%
    rowwise() %>%
    mutate(
      expanded = list(seq(start_num, end_num))
    ) %>%
    unnest(expanded) %>%
    mutate(bildenummer = as.character(expanded), is_alle = FALSE) %>%
    select(-expanded)

  # Combine the datasets
  bind_rows(alle_cases, non_ranges, ranges)
}

lyngstad_expanded <- expand_number_ranges(lyngstad_clean)

# Filter the original polygons to match the reference data ####
# Split reference data into "alle" records and specific records
lyngstad_specific <- lyngstad_expanded %>%
  filter(is_alle == FALSE) %>%
  select(nib_project_id, stripe_clean, bildenummer) %>%
  # Convert bildenummer to character to ensure matching
  mutate(bildenummer = as.character(bildenummer))

lyngstad_alle <- lyngstad_expanded %>%
  filter(is_alle == TRUE, stripe_clean != "alle") %>%
  select(nib_project_id, stripe_clean)

lyngstad_alle_alle <- lyngstad_expanded |>
  filter(is_alle == TRUE, stripe_clean == "alle") |>
  select(nib_project_id, stripe_clean)

# Process specific matches (join on project, stripe, number)
matches_specific <- polygons_clean %>%
  semi_join(
    lyngstad_specific,
    by = c(
      "nib_project_id" = "nib_project_id",
      "stripenummer_clean" = "stripe_clean",
      "bildenummer" = "bildenummer"
    )
  )

# Process "alle" matches (join on project and stripe only)
matches_alle <- polygons_clean %>%
  semi_join(
    lyngstad_alle,
    by = c(
      "nib_project_id" = "nib_project_id",
      "stripenummer_clean" = "stripe_clean"
    )
  )

# Additionally, include all stripes for projects where "alle" is specified
matches_alle_alle <- polygons_clean %>%
  semi_join(
    lyngstad_alle_alle,
    by = c(
      "nib_project_id" = "nib_project_id"
    )
  )

# Combine sets of matches
polygons_filtered <- bind_rows(
  matches_specific,
  matches_alle,
  matches_alle_alle
) %>%
  # Remove duplicates that might occur if a row matches both criteria
  distinct(nib_project_id, stripenummer, bildenummer, .keep_all = TRUE) %>%
  # Remove temporary columns used for matching
  select(
    nib_project_id,
    prosjektnavn,
    pixelstorrelse,
    stripenummer,
    stripenummer_clean,
    bildenummer,
    shape
  )

# Result: polygons_filtered now contains all rows that match the reference data,
# including proper handling of "alle" entries and proper standardization of stripe numbers

# Display summary statistics ####
polygons_filtered |>
  st_drop_geometry() |>
  group_by(prosjektnavn, nib_project_id) |>
  count()

# Add orthofoto coverage in Lyngstad 2016 ####

orthofoto_2016 <- st_read(
  "data/Lyngstad/Lyngstad2016-Fig1.gpkg",
  layer = "searched-orthophoto"
) |>
  st_transform(st_crs(polygons_filtered)) |>
  st_set_geometry(attr(polygons_filtered, "sf_column"))

polygons_filtered <- bind_rows(
  polygons_filtered,
  orthofoto_2016
)
tail(polygons_filtered)

# Plot the geometries ####

# polygons_filtered %>%
#   st_union() %>%
#   st_cast("POLYGON") %>%  # Cast to individual polygons
#   st_combine() %>%        # Recombine to a geometry collection
#   plot()
# Above approach doesn't resolve internal boundaries/overlaps, maybe because invalid geometry?

# polygons_filtered %>%
#   st_make_valid() %>%     # Ensure geometries are valid
#   st_buffer(0) %>%        # Minor buffer to fix tiny gaps/overlaps
#   st_union(is_coverage = TRUE)  |>
#   st_geometry() |>
#   plot()
# Above approach maintains overlaps but not differentiated by project

footprint_projects <- polygons_filtered |>
  select(nib_project_id, prosjektnavn) |>
  group_by(nib_project_id, prosjektnavn) |>
  summarize()
footprint_projects |>
  st_geometry() |>
  plot()

footprint_total <- footprint_projects |>
  st_union() |>
  st_geometry()
footprint_total |>
  plot()

# Write to file ####
st_write(
  polygons_filtered,
  "data/DMraisedbog.gpkg",
  layer = "nib-lyngstad-polygons",
  append = FALSE
)

st_write(
  footprint_projects,
  "data/DMraisedbog.gpkg",
  layer = "nib-lyngstad-projects",
  append = FALSE
)

st_write(
  footprint_total,
  "data/DMraisedbog.gpkg",
  layer = "nib-lyngstad-footprint",
  append = FALSE
)

# sessionInfo ####

sessioninfo::session_info()
