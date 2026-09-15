# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> Scope: durable conventions and architecture orientation only. Project *state*
> (what's being tried, what was superseded, active direction) lives in `notebook.md`
> (dated research log, newest entry on top) — not here.

## Project Overview

R-based research project modelling the distribution of raised bogs in Norway under
climate-change scenarios. One balanced 3-class random forest is trained on pooled
Norwegian and European data and projected over Norway, with a focus on feature-space
data partitioning for honest cross-validation, skill-versus-novelty reliability curves,
and extrapolation (area-of-applicability) analysis.

## Development Environment

- **Language**: R (prefer the native pipe `|>`)
- **Project Type**: RStudio project (DMraisedbog.Rproj)
- **Coordinate systems**: EPSG:25833 (UTM 33N) for Norwegian data; EPSG:3035 for modelling
- Core install: `install.packages(c('tidyverse', 'sf', 'terra'))`

## Key Dependencies

- **Spatial analysis**: `sf`, `terra`, `rnaturalearth`, `units`
- **Data manipulation**: `tidyverse` (readr, dplyr, tidyr, purrr, tibble, stringr, forcats), `foreign`
- **Learner**: `randomForestSRC` (production); `randomForest` and `glmnet` only as
  cross-implementation checks on the variable-importance ruler; `MIAmaxent` for FOP plots
- **Distances**: `FNN` (nearest neighbours behind the dissimilarity index; the metric
  itself is implemented in `R/functions.R`, not taken from `CAST`)
- **Elevation (outside Norway)**: `elevatr`
- **Visualization**: `ggplot2`, `ggrepel`, `scales`, `patchwork` (manuscript figures only)
- **Run infrastructure**: `rmarkdown` (RUNALL knits), `sessioninfo`, `tictoc`
- **Retired (scripts in `R/archive/`)**: `sabinaNSDM`, `biomod2`, `GeoThinneR`, `cluster`,
  `dbscan`, `RANN`, `ROCR`, `twosamples`, `furrr`

## Data Architecture

### Data Sources
- **Lyngstad dataset**: Historical raised bog mapping from Norwegian regions (2012-2023)
- **NiB (Norge i Bilder)**: Orthophoto database; defines the spatial extent of the Lyngstad survey
- **CHELSA**: Climate data (past: 1981-2010, future: 2071-2100)
- **EUNIS**: European habitat classification data
- **Natura 2000**: Protected-areas database
- **DTM50**: Norwegian digital terrain model at 50 m resolution (Norway only)

### Pipeline
`R/RUNALL.R` is the table of contents: it lists every script in dependency order, with a
one-line description each, and knits them to HTML. **If a script is not listed there it is
not part of the analysis** — retired and exploratory scripts live in `R/archive/`, and
RUNALL prints a warning for anything in `R/` that is neither listed nor archived.

Stages: `pl0_` shared upstream (raw data, predictor stacks, the two labelled training
blocks), `pl2_` the production pipeline, `pl3_` reliability mapping (parallel to pl2, not
nested). **There is no `pl1_`** — those were the retired hierarchical cascade; the
pl2/pl3 names are kept so existing references in `notebook.md` and
`plan_EUintegration.md` stay valid.

### Key Data Files
- `data/DMraisedbog.gpkg`: Main GeoPackage of processed spatial layers
- `data/CHELSA/`: Climate predictor variables (past and future)
- `data/Lyngstad/`: Raw shapefiles from bog mapping projects
- `data/NiB/`: Orthophotos defining survey extent
- `output/predictors_*.tif`: Processed predictor rasters at different scales/projections
- `output/*_coords_*.csv`: Extracted presence/absence coordinates
- `output/pl0/`: EU domain, EPM mask, EU terrain, land-use screen, and the two training blocks
- `output/pl2/`: pooled modelling frame, partitions, CV results, production model, projections
- `output/pl2/di_ruler_production.rds`: the frozen DI ruler (weights + scaling + normalisation constant)
- `output/*_parts/`: cached per-country / per-tile warps — expensive to rebuild, do not delete

### Coordinate Reference Systems
- Norwegian data: EPSG:25833 (UTM 33N)
- Modelling space: EPSG:3035
- Global data: EPSG:4326 (WGS84) — transformed as needed

### Spatial Data Patterns
- Extensive use of `st_union()`, `st_intersection()`, `st_transform()`
- Polygon validation with `st_make_valid()` / `st_is_valid()`
- Area calculations with `st_area()` and explicit unit conversions
- Attribute filtering (e.g. `M_TYPE` containing "A" for raised bogs)

## Machine Learning Modeling Framework

**One pooled model, not a nested cascade.** European and Norwegian data are row-bound into
a single training frame with a shared 3-class response (`nonpeat` / `otherpeat` / `bog`),
and a `dataset` column that is bookkeeping only and must never enter a model as a
predictor. The hierarchical `rf_global` covariate was retired — see
`plan_EUintegration.md` section 1.1.

1. **Learner**: one balanced 3-class `rfsrc()` (`case.wt = make.wt(y)`,
   `sampsize = make.size(y)`), argmax labels, no threshold to optimise. RFQ is two-class
   only and is not used.
2. **Data partitioning**: feature-space partitioning for honest cross-validation, cut along
   the top *materially shifting* predictor
3. **Projection**: Norway only, current and future, as a 3-layer probability simplex per
   scenario
4. **Reliability / extrapolation analysis**: DI-based, against a ruler frozen once from
   the production fit

### Key Algorithms and Parameters
- **Learner**: `randomForestSRC::rfsrc()`, balanced 3-class via `make.wt` / `make.size`
- **Distance metric**: VI-weighted Euclidean DI over *projection-dynamic predictors only*
  — static axes (terrain, paleo) cannot register projected change and would only dilute it
- **Parallelization**: `furrr` (multi-core)

## Working with This Codebase

### Data File Locations
- Large datasets in `data/` (gitignored); outputs in `output/` (gitignored)
- Processed spatial data consolidated in `data/DMraisedbog.gpkg`

### Common Operations
- Reading shapefiles: `st_read()` with explicit CRS handling
- Processing multiple files: `map()` over file lists
- Spatial filtering: combine attribute filters with spatial operations
- Writing outputs: `st_write()` with append options for GeoPackage layers

### Data Quality Considerations
- Validate geometries before spatial operations
- Handle CRS mismatches explicitly
- Check for spatial overlaps and duplicates in polygon datasets
- Apply appropriate quality filters (e.g. `DATAQUALITY == "G"` for Natura 2000)
- **Absence sampling**: both blocks draw absences by the same climate-stratified rule
  (`draw_stratified_absences()` in `R/functions.R`) and pass the same land-use / water
  screen; presences are one row per 250 m cell. Spatial thinning was retired with the
  archived `pl0_spatiallyThinCells.R`.

## File Patterns to Recognize

- `pl0_`/`pl2_`/`pl3_` prefixes: pipeline stage (see `R/RUNALL.R`); no `pl1_` — retired
- `R/archive/`: retired or exploratory scripts, deliberately outside the pipeline
- Scripts prefixed `collate*`: data integration and standardization
- `.Rproj`: RStudio project configuration
- `.gpkg`: GeoPackage spatial database
- `.tif`: raster data (climate, terrain, predictions)
- `.dbf/.shp`: legacy shapefiles requiring CRS handling

## Code Style and Conventions

- **Native pipe** `|>` preferred.
- **Pipe usage**: use pipes inside function-call arguments only sparingly.
  - Prefer: `coords |> select(x, y) |> write_csv("output/coords.csv", append = FALSE)`
  - Over: `write_csv(coords |> select(x, y), "output/coords.csv", append = FALSE)`
- **Messages**: use `cat()` sparingly — for critical checkpoints only; otherwise rely on
  good workflow comments.
- **sessionInfo**: every script ends with a `# sessionInfo ####` top-level heading so it
  is documented when RUNALL.R knits to HTML.

## Deliverables

- **Do not create shareable artifacts unless explicitly asked.** Written output belongs in
  the repo as plain files — `notebook.md` for dated project state, a markdown file at the
  root for standing references (e.g. `plan_EUintegration.md`), knitted HTML via `RUNALL.R`
  for script output. Deliver in the terminal and in the repo; publishing to an external
  host is opt-in, on request only.

## Git Workflow for Agents

Several agents may work on this repository at once. `main` is an integration branch that
only JV touches; everything else happens on task branches.

- **Branch and worktree per task.** Work on a branch in a git worktree (`EnterWorktree`),
  never directly on `main` in the main checkout. One task per branch; name it after the
  task, not the date.
- **Commit before stopping.** A session can end at any moment, and an uncommitted change
  is invisible to everyone else. Small commits, each a complete step, with a message that
  says what changed and why. Do not leave work only in the working tree.
- **Never merge into or push `main`.** Report the branch name and the commits when done;
  JV reviews and merges. Pushing a task branch is fine when asked.
- **`output/` and `data/` are shared state, not branch state.** They are gitignored and
  every worktree reads the same folders. Do not run pipeline scripts while another run is
  in progress; check for running `Rscript` processes first, and record any rerun in
  `notebook.md` (what was rerun, from which script, and why).
- **`notebook.md` is the coordination log.** Newest entry on top. Add an entry for any
  decision, rerun, or finding another agent would need to know about, and read the top
  entries before starting.
- **Manuscript figures** are built by `ms/make_figures.R`, run from the repository root.
  It reads pipeline outputs only and is deliberately outside `R/RUNALL.R`.

## External Tools

- Spatial CLI tools (e.g. `ogrinfo`) are available in bash.
- When running R from the CLI, use `Rscript --vanilla` (startup otherwise auto-loads
  mcptools and halts).

## Known Issues

**randomForestSRC + Positron crash with `do.trace`**
- **Problem**: R crashes (exit code -1073740791) when `randomForestSRC::imbalanced()` is
  called with `do.trace = TRUE` (or any numeric value) inside functions wrapped with
  `purrr::map_dfr()` in the Positron IDE.
- **Symptom**: "R exited unexpectedly" / "Error refreshing variables: The language runtime
  exited before RPC completed".
- **Solution**: always use `do.trace = FALSE` when calling `imbalanced()` inside functions
  or loops.
- **Root cause**: console trace output conflicts with Positron's variable-refresh
  mechanism in nested function contexts.
