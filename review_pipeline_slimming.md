# Pipeline slimming review — 2026-09-14

Scope: every script listed in `R/RUNALL.R` (13 pl0 + 12 pl2), `R/functions.R`, and what
`ms/make_figures.R` and `ms/ms.qmd` actually consume. Nothing in `R/archive/` was reviewed.
Status: section 0 resolved (manuscript moved to the pairwise arm, notebook 2026-09-14);
A1 and A2 implemented on this branch. The rest are recommendations.

Two aims, as asked: (A) what to remove so the codebase is slimmer, and (B) what to add so a
fresh RUNALL is easier to review.

---

## 0. Resolve before the next RUNALL: the calibrating CV arm

`R/pl2_fitErrorProfiles.R` sets `CV_ARM <- "pairwise"` and its header argues at length that
pairwise is the arm that calibrates (with LOPO kept as `COMPARE_ARM`). That decision is
recorded nowhere else:

- the header cites "notebook 2026-09-14", but `notebook.md` has no entry after 2026-09-12,
  and that entry says the *LOPO* switch "is still right";
- the commit that made the change (`1bb69db`, "Read reliability on the signed bio10 offset
  axis") has a one-line message;
- `ms/ms.qmd` still describes the reliability curves, Fig 3, Table 4 and S1 Fig as
  "held-out cells of the leave-one-partition-out arm" (lines 391, 628, 650, 757).

So the code, the notebook and the manuscript currently say three different things about which
arm the map reads. Whichever is right, it needs a dated notebook entry and a one-line fix in
the other two places. This is the single most important thing for the review of the next
RUNALL, because every reliability number downstream depends on it.

---

## A. What I would remove

### A1. `R/functions.R`: 1,193 of 2,058 lines have no live caller

Verified by grepping every pl0/pl2 script and `ms/make_figures.R` for each function name.

| Function (line) | Era | Status |
|---|---|---|
| `create_cut_based_partitions` (372) | binary, W_CV tuning | dead |
| `create_spacer_based_partitions` (502) | binary, W_CV tuning | dead |
| `evaluate_spacer_partitioning` (689) | binary, W_CV tuning | dead |
| `evaluate_cut_partitioning` (805) | binary, W_CV tuning | dead |
| `calculate_prediction_distances` (1065) | W_CV tuning | dead |
| `calculate_cv_distances` (1142) | W_CV tuning | dead (only called by the two dead evaluators) |
| `calculate_gmean` (1213) | binary | dead (only called by dead `train_*`) |
| `get_metrics_at_threshold` (1361) | binary thresholds | dead |
| `calculate_classification_thresholds` (1458) | binary thresholds, ROCR | dead |
| `train_brf_model` (1545) | binary BRF tuning | dead; also references an undefined global `progress_file` |
| `train_rfq_model` (1662) | binary RFQ tuning | dead |
| `train_single_model` (1756) | binary | dead |

Delete all twelve. What survives is exactly the set of helpers the 3-class pipeline uses:
`landuse_screen_keep`, `identify_dynamic_predictors`, `compute_feature_scaling`,
`calculate_weighted_di`, `calculate_multiclass_metrics`, `calculate_ovr_auc`,
`partition_by_presence_sorting`, `allocate_across_strata`, `draw_stratified_absences`, plus
`featuredist` if A2 keeps its one consumer (see below).

The `library()` block at the top of functions.R (readr, dplyr, terra, cluster, tidyr,
ggplot2, sf, purrr, twosamples) can go with them: every live helper already namespaces its
calls (`dplyr::`, `purrr::`, `tibble::`, `terra::`, `FNN::`, `stats::`). `cluster` and `sf`
are not referenced by any live helper at all. Each pipeline script loads what it needs itself,
so nothing downstream depends on functions.R attaching packages — but this is the one removal
I would confirm with a full RUNALL rather than by inspection.

### A2. Two diagnostic scripts to archive

**`R/pl2_exploreFeatureSpaceDistances.R`** (W_sample per feature). Superseded: the quantity it
measures (how far each predictor shifts current → future) is now computed as `shift_sd` by
`identify_dynamic_predictors()` and written into `weights_feature_data_partitioning.csv`,
which is what `pl2_exploreOccupancy.R` and the partition script actually read. The
importance-vs-shift picture is drawn in `pl2_freezeDIRuler.R` (`ruler_importance_vs_shift.png`).
Its output `Wsample_by_feature.csv` has no consumer and W_sample is not mentioned in the
manuscript. Archiving it also retires `featuredist()` (the last live user of `twosamples`
and `furrr`) and the slowest loop in functions.R (one FNN call per point).

**`R/pl2_exploreWeightedPCA.R`**. Not consumed by anything and not in the manuscript. It
carries a vestigial EU-anchor branch (`eu_presence_predictors.csv`, written by the archived
`pl3_extractEUpredictors.R`) that can never fire. The question it answers — how far the
future domain sits from training — is answered on the frozen ruler by the DI coverage check in
`pl2_evaluate.R` (`cv_novelty_coverage.png`), which is the version the argument rests on.
If you want to keep a PCA picture, strip the EU-anchor block and the `weighted_pca_*.csv`
outputs; otherwise archive.

### A3. Sections to delete inside live scripts

**`pl2_weightFeaturesDataPartitioning.R`** — the most expensive removal candidate. It fits four
rulers; every consumer filters `method == "Balanced Random Forest"`:

- *Norway-only rfsrc* (5 reps). Its stated purpose is "so the pl3 skill-vs-novelty curve can be
  recomputed under weights that never saw an EU row". pl3 is gone; no consumer. Delete.
- *randomForest* (3 reps, single-threaded permutation VI) and *glmnet* (3 × cv.glmnet
  multinomial over 85k × 43). Used only for the in-script Spearman rank table and one facet
  of one plot. Neither appears in the manuscript. Delete, or gate behind a `CROSS_CHECK <-
  FALSE` flag if you want them available for a reviewer response.
- The FOP / `pairs()` block (`MIAmaxent::plotFOP` on a binarised bog-vs-rest response). A
  binary diagnostic bolted onto a 3-class problem; nothing reads it. Delete.
- Dependencies that then disappear: `randomForest`, `glmnet`, `MIAmaxent`, `tictoc`.
- Header comment references `pl2_dissimilarityindex.R` and "pl3" — both retired.

**`pl0_rasterizePresenceAbsence.R`** — tables 3 and 4 (presences on the 5 km global grid,
Norway and EU) are relics of the retired global model. Only tables 1 and 2 are read (by
`pl0_labelNorwayBlock.R`). Drop tables 3–4, the `presences_eu` read, and the two
`presence_coords_global_*.csv` outputs; `global_raster` is then only needed for the CRS.

**`pl0_buildEPMmask.R`** — the "5 km version" section (`epm_map_cat_5km.tif`,
`epm_coverage_5km.tif`) is described as for "domain-level diagnostics and climate
stratification", but nothing reads either file. Delete the section and the two outputs.

**`pl0_buildEUdomain.R`** — `eu_domain_5km.tif` is written and never read (the domain is
consumed as the `EU_domain` gpkg layer). Keep the rasterisation if you want the 5 km cell
count in `eu_domain_summary.csv`; drop the `writeRaster`. The "warm-flank cells above the
legacy 17.0 threshold" row is a comparison against a threshold no longer used anywhere.

**`pl0_collatePredictors.R`** —
- `memory_usage()` is defined and never called.
- `ar50_250m_cover_EPSG3035.tif` is written and never read (artype_60 enters the regional
  stack directly). Drop the write.
- The past/future "scaling check" (two `plot()` calls over `checkvars`) is a one-time visual
  check whose conclusion is a comment. Replace with an assertion (`stopifnot(length(checkvars)
  == 0)`), or drop.
- The NA-to-zero block for threshold variables appears three times (global, regional current,
  regional future) as ~35-line copies. One helper `fill_threshold_na(stack, vars)` would make
  the three calls reviewable at a glance and remove the asymmetry that `gdd5` is in the
  regional list but not the global one — which I assume is deliberate but is currently only
  visible by diffing the two vectors.

**`pl0_prepareCHELSATrace.R`** — `dir.create("output/pl0/chelsa_trace")` creates a directory
nothing uses. The four `calc_*` helpers are two functions each written twice (mean/sum,
consecutive/all); one `calc_stat(x, fun, consecutive)` would do. `n_consecutive_icefree` is
written and then explicitly excluded downstream as collinear — fine to keep on disk, but say
so at the write rather than at the read.

**`pl0_collateEuropeanRaisedBog.R`** — `type_inventory <- "output/pl3/epm_peatl_type_total.csv"`
points at a directory that no longer exists, so the "which peatl_type values matched" report
silently never prints. Either regenerate that inventory inside this script (it is a useful
provenance table) or delete the block. The comment pointing at `scratch_exploreEPM2025.R`
should say `R/archive/`.

**`pl2_evaluate.R`** —
- The "Skill vs novelty" section (fold-level `Gmean_macro ~ DI_mean` plot, per-row DI
  ventiles, `accuracy ~ DI_mid` plot) is the pre-`pl2_fitErrorProfiles.R` version of the
  reliability curve and is superseded by it. Keep the DI *coverage* check (load-bearing) and
  the "where the CV's most novel bog cells sit" table (it motivates the offset axis); drop the
  rest.
- `future_rows` and `train_rows` re-read the partitioned CSV twice although `mf` is already
  in memory.
- `train_avg_dist` is now a frozen constant, so the `train_avg_dist` column in the metrics
  table and `mean_train_avg_dist` in `DI_summary` carry no information.

**`pl2_partitionDataByTopFeature.R`** — `partition_topfeature.tif` is written and never read;
keep the on-screen map, drop the file. The "observations dropped / drop rate" reporting is
always zero since the complete tiling, so it can shrink to the single assertion
`stopifnot(partitioning_result$n_dropped == 0)`.

**`pl0_collateLyngstadExtent.R`** — two commented-out `plot` attempts and the pasted warning
text can go; `%>%` → `|>` for consistency.

### A4. Orphan outputs (written by a live script, read by nothing)

`output/presence_coords_global_no.csv`, `output/presence_coords_global_eu.csv`,
`output/ar50_250m_cover_EPSG3035.tif`, `output/pl0/eu_domain_5km.tif`,
`output/pl0/epm_map_cat_5km.tif`, `output/pl0/epm_coverage_5km.tif`,
`output/pl2/Wsample_by_feature.csv`, `output/pl2/weighted_pca_*.{csv,png}`,
`output/pl2/partition_topfeature.tif`, `output/pl2/weights_feature_by_class.csv` (printed in
the knit; the per-class VI the manuscript uses comes from `variable_importance_brf.csv`).

Small summary CSVs that are written for the record and read only by eye
(`di_ruler_summary.csv`, `ruler_weight_mass.csv`, `eu_domain_area_by_country.csv`,
`landuse_parts_index.csv`, `occupancy_*` beyond `occupancy_grid.csv`) are fine to keep — they
are the review trail — but see B1 for making that trail explicit.

### A5. Documentation that has drifted

- `CLAUDE.md` "Key Dependencies" lists `probably`, `yardstick`, `rsample`, `pROC`, `ROCR`,
  `cluster`, `dbscan`, `RANN`, `CAST`, `GeoThinneR` and `randomForest`. After A1–A3, none of
  these is loaded by a live script (`CAST` survives only as a comment; `FNN`, `elevatr`,
  `ggrepel`, `sessioninfo` are the actual extras). The "Spatial thinning with GeoThinneR"
  bullet under Data Quality describes the archived `pl0_spatiallyThinCells.R`.
- `R/RUNALL.R` stage comment still says reliability mapping is "(to come)".
- `pl0_buildEUterrain.R` header explains itself by contrast with `pl3_extractEUpredictors.R`
  (archived); fine as history, but say "archived".
- The notebook's "Pipeline is 24 scripts" (2026-08-20) is 25 in RUNALL; after A2 it would be
  23.

### A6. Rough effect

| | before | after A1–A3 |
|---|---|---|
| `R/functions.R` | 2,058 lines | ~860 |
| live pl2 scripts | 12 | 10 |
| model fits in `pl2_weightFeatures…` | 16 (5+3+3+5) | 5 |
| extra packages | randomForest, glmnet, MIAmaxent, tictoc, furrr, twosamples, cluster | none of these |

The weights script alone should drop from "the expensive step" to a quarter of its runtime.

---

## B. Diagnostics to add for reviewing the next RUNALL

Ordered by how much review effort each saves.

### B1. A run manifest written by `RUNALL.R` (no per-script edits)

Record, per script: start, end, duration, and every file under `output/` whose mtime falls
inside that window. Write `output/run_manifest.csv`. This makes three things automatic that
today need grepping: which script produced which file, which files no script touched (orphans,
i.e. A4 next time), and whether a script ran for a plausible time. Cost: ~15 lines in the
render loop.

### B2. Input-freshness guards in the consumers that have been burned

On 2026-08-24 `pl2_fitErrorProfiles.R` and `pl2_mapReliability.R` ran on stale predictions
and produced a full set of wrong numbers that read as results. A two-line guard prevents the
recurrence:

```r
stopifnot(file.mtime("output/pl2/predictions_cv_topfeature.csv") >
          file.mtime("output/pl2/modeling_frame_regional_partitioned_topfeature.csv"))
```

and the analogue in `pl2_mapReliability.R` (`error_profiles.rds` newer than predictions;
`di_ruler_production.rds` newer than the weights file). Alternatively do it once in RUNALL
from a declared dependency list, but the explicit form in the two reliability scripts is
what a reader will look for.

### B3. Per-predictor ranges by block in `pl2_createModelingFrame.R`

Write `modeling_frame_predictor_ranges.csv`: min / q01 / median / q99 / max of every predictor,
split by `dataset` and by `scenario`. The gsp sentinel (4.29e8 mm in 8,987 rows, found
2026-09-12) would have been visible in this table on the first run. Add a warning when a
predictor's EU and NO ranges do not overlap at all, which is the one pattern that is never
legitimate for climate variables on a shared source.

### B4. Diff against the previous run

Before rendering, RUNALL copies the small summary CSVs (`modeling_frame_summary`,
`di_ruler_summary`, `metrics_cv_topfeature`, `cv_novelty_coverage`, `error_profiles_by_bin`,
`reliability_summary`, `reliability_at_lyngstad`, `occupancy_grid`) to `output/_previous/`;
after rendering it prints a numeric diff of each. Unchanged code should produce an empty diff
(see B5); a changed number then points at exactly the script to re-read.

### B5. Seed the fold fits in `pl2_evaluate.R`

`rfsrc()` is called 25 times with no seed and no `set.seed()` (the production fit in
`pl2_predict.R` uses `seed = -42`; the weights script seeds each rep). Without a per-fold
seed, two RUNALLs of identical code give different CV numbers, so B4's diff is always
non-empty and the review cannot separate code changes from forest noise. One line per fold:
`seed = -(1000 + i)`.

### B6. Row-accounting assertions along the chain

Cheap `stopifnot`s that each consumer should carry, so a silent row loss stops the run
instead of surfacing three scripts later:

- `pl2_createModelingFrame.R`: rows in `no_block_coords.csv` = NO rows in the frame + the
  count reported as dropped; same for the EU block.
- `pl2_partitionDataByTopFeature.R`: partitioned rows = current-scenario rows (already true
  by construction; assert it).
- `pl2_evaluate.R`: LOPO predictions = n_train rows, each (x, y, dataset) scored exactly
  once; pairwise = each row scored k−1 times.
- `pl2_fitErrorProfiles.R`: already asserts the join is 1:1 — good; add
  `stopifnot(all(fold_ref$ref > 0))` so an empty training-bog set in a fold cannot yield an
  NA reference silently.

### B7. Save what is currently only printed

Three tables the manuscript or a reviewer will want that exist only in knitted HTML:

- `pl2_predict.R`: OOB confusion and per-class metrics → `oob_metrics_production.csv`.
- `pl2_evaluate.R`: per-fold composition (n by class × dataset, train and test) →
  `cv_fold_composition.csv`. Table 3's "median n_train 5,042 vs 77,317" argument rests on it.
- `pl2_interpret.R`: the transition table → `lyngstad_transitions.csv`
  (`make_figures.R` currently recomputes it from the gpkg).

### B8. One place for the frozen constants

The design knobs are scattered: `ABSENCES_PER_BOG_PER_CLASS` and the land-use thresholds in
functions.R; `k_partitions`, `NTREE`, `N_BINS`, `N_BOOT`, `CV_ARM`, `NORM_SAMPLE` in five
scripts; the scenario and GCM as strings in two pl0 scripts (`ssp370`, `GFDL-ESM4`, which
must agree or the future stack silently goes empty). A `R/config.R` sourced by every script,
or at minimum a `pipeline_settings.csv` that each script appends its constants to, lets the
reviewer confirm the configuration of a run from one file. The scenario/GCM pair is the one
that most deserves this: `pl0_downloadCHELSA.R` and `pl0_collatePredictors.R` each define it
independently.

### B9. Knit HTML into one folder

`render()` writes each `.html` next to its script in `R/`. `output_dir = "output/html"` (with
a timestamped subfolder per RUNALL) keeps the review artefacts of one run together and out
of the source tree.

---

## Suggested order

1. Section 0 (arm decision, notebook entry, manuscript wording) — before anything else.
2. A1, A2, A3 (weights script first: biggest runtime and dependency win), A5.
3. B5 and B2 (small, and they protect the run you are about to do).
4. B1, B3, B6, B7 (the diagnostics that make the run reviewable).
5. RUNALL.
6. B4, B8, B9 as time allows.
