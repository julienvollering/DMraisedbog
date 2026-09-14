# Manuscript plan and checklist

Working plan for `ms/ms.qmd`. Started 2026-09-12. Everything here describes the pipeline as
listed in `R/RUNALL.R` on that date (24 scripts, pl0 + pl2, frozen DI ruler, reliability read
on the signed bio10-offset axis). Anything superseded before that date is deliberately absent.

## Audience and stance

Ecologists and ecological modellers. Lead with the model-free evidence (where the future
climate of today's raised bogs sits in climate space, and what is observed there today), then
the model, then the cross-validation that says how far the model can be trusted. Methods are
written so a reader can follow the *reasoning* for each design choice without reading code.

## Materials and methods -- section outline

1. Study system and question
2. Raised-bog occurrence data
   2.1 Norway: aerial-photo survey (Lyngstad), survey footprint, 250 m cells
   2.2 Europe: Natura 2000 + EUNIS presences, mire region IV domain
3. A three-class response: raised bog between two opposite non-bog classes
   3.1 Norway labels (AR50 peatland fraction)
   3.2 Europe labels (European Peatland Map 2025)
   3.3 One land-use / water screen for both blocks
   3.4 Climate-stratified absence draw
4. Predictors (climate, terrain, deglaciation history; current and 2071-2100 SSP3-7.0)
5. Model-free description of the projection problem (occupancy table)
6. The model: one balanced three-class random forest
7. Evaluating transfer into a warmer climate
   7.1 Partitions cut along summer temperature; pairwise and leave-one-partition-out arms
   7.2 Skill as a function of warming beyond the training bogs (offset axis, bins, cluster bootstrap)
   7.3 Multivariate novelty and area of applicability (frozen DI ruler)
8. Projection and reliability mapping
9. Software

## Results -- section outline

1. Training data (Table 1)
2. Where the future takes today's raised bogs, and what is observed there (Fig 2, Table 2)
3. Model fit and cross-validated skill (Table 3)
4. Skill against warming beyond the training bogs (Fig 3, Table 4)
5. Projected change at the mapped raised bogs (Fig 4, Fig 5, Table 5)
6. Reliability of the projection where it is read (Fig 6)

## Figures and tables (built by `ms/make_figures.R`, run from the repo root)

- [x] Fig 1  study area: EU domain + presences, Norway survey footprint + raised bogs
- [x] Fig 2  climate space (bio10 x gsp) with classes and future bog cells; occupancy grid
- [x] Fig 3  cross-validated skill vs signed bio10 offset from training bogs
- [x] Fig 4  P(bog) and majority class, current and future, over Norway
- [x] Fig 5  change at the 564 Lyngstad raised bogs
- [x] Fig 6  reliability layers over the future domain (offset, DI, expected recall)
- [x] S1 Fig skill vs DI (the retained multivariate axis)
- [x] Table 1 training frame composition
- [x] Table 2 occupancy grid
- [x] Table 3 cross-validated skill by arm and data block
- [x] Table 4 skill by offset bin
- [x] Table 5 Lyngstad transitions
- [x] S1 Table predictors

## Writing checklist

- [x] Delete template material from ms.qmd; point bibliography at Repeat-raisedbogDM.bib
- [x] Materials and methods drafted with subheaders
- [x] Results drafted with subheaders
- [x] TODO markers wherever input from JV is needed (search for `TODO`)
- [x] Render check: HTML variant renders clean (one missing citation key, Tegetmeyer2025); PDF blocked by local TinyTeX 2025 vs CTAN 2026 (`changepage.sty`), run `update-tlmgr-latest` or reinstall TinyTeX
- [ ] JV review of TODOs
- [ ] Introduction, Discussion (not started; out of scope for this pass)

## Known stale points to watch

- `pl2_mapReliability.R` was rerun in the main checkout on 2026-09-12 (finished 15:09);
  the manuscript numbers and figures use that run (expected-recall interval at the future
  offset now correctly NA).
- The knitted HTML for `pl2_fitErrorProfiles.R` and `pl2_mapReliability.R` predates the
  offset axis; the CSV/RDS outputs dated 2026-09-12 are current and are what the manuscript
  uses.
- DATA QUALITY: growing-season precipitation (gsp) carried CHELSA's no-data code in 8,987
  training rows (cold, high-elevation non-peat): the gsp GeoTIFFs are unsigned 32-bit with
  no NoData tag, so the sentinel 4294967295 read as 4.29e8 and bilinear resampling blended
  it into neighbours. FIXED 2026-09-14 in `R/pl0_collatePredictors.R` (sentinel masked on
  the cropped native grid before `project()`, gsp added to both NA-to-zero lists, an
  `assert_plausible()` guard on every CHELSA crop). PIPELINE NOT YET RERUN: run
  `R/RUNALL.R` from `pl0_collatePredictors.R` onward (the block draws are seeded and read
  only bio10, so training cells stay the same; the DI ruler, AOA threshold and the
  cold-wet occupancy cell will move), then `ms/make_figures.R`, then clear every
  `TODO(gsp)` in `ms.qmd`.
