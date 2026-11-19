# Plan: how generous is the bog map, and what counts as "no longer suitable"?

Status: implemented and run 2026-10-07 as `R/pl2_calibrationCheck.R` (results below, under
"Results 2026-10-07"; to be copied into `notebook.md` as entry 2026-10-07 (c): this session could
not write to the shared checkout). Background: `notebook.md` 2026-10-07 (b).

## Why

The production forest is trained at equal class priors, so argmax over-predicts bog against
real prevalence (OOB: 7,597 predicted vs 1,615 true, precision 0.19), and balanced P(bog)
is a likelihood ratio, not an occupancy probability. The manuscript currently calls a cell
"bog" when argmax says so and reports loss as flips out of that class. Two questions are
open, and the Lyngstad survey footprint (complete mapping, so prevalence is known there)
lets us answer both with existing rasters.

1. **Commission.** Inside the footprint, how many of the present-day argmax-bog cells are
   mapped bog? This is a real precision, replacing the 0.19 that comes from the designed
   training prevalence.
2. **Loss criterion.** Which future cells fall below the conditions where bogs occur today,
   as opposed to merely flipping class at equal priors?

## Inputs (all exist; no refit, no new prediction)

- `output/pl2/rf_local_pred_brf.tif` (`current_*`, `future_*`) and
  `rf_local_pred_brf_ssp126.tif` (`future_ssp126_*`)
- footprint: `nib-lyngstad-footprint` layer of `data/DMraisedbog.gpkg`;
  presences/absences by cell in `output/presence_coords_regional.csv` and
  `output/absence_coords_regional.csv` (`pl0_rasterizePresenceAbsence.R`)
- `output/pl2/reliability_*.tif` for the AOA mask

## Step A: commission inside the footprint

1. Rasterise the footprint on the 250 m prediction grid; keep cells with a prediction.
2. Mapped bog = presence cells; prevalence pi = n_presence / n_footprint_cells
   (state the denominator: cells with partial polygon coverage are included, see the
   comment in `pl0_rasterizePresenceAbsence.R`).
3. Argmax current class per cell. Report a 2 x 2: predicted bog vs mapped bog, giving
   recall, precision, commission (predicted-bog cells that are not mapped bog), and the
   ratio of predicted to mapped bog area. Also split by `otherpeat`, because a bog call
   on mapped other peatland is a different error from one on non-peat.
4. Threshold sweep for context: precision/recall of argmax against the same ratio at the
   balanced-P(bog) cut that matches the mapped area (the "area-matched" cut).

Output: `output/pl2/calibration_footprint.csv`.

## Step B: calibrated P(bog) by prior shift

1. Prior shift on the three balanced probabilities: p_cal(c) proportional to
   p_bal(c) x prior(c) / (1/3), renormalised, with priors from the footprint's mapped
   bog, other-peat and non-peat shares. Document which map defines `otherpeat` there.
2. Report calibrated P(bog) at present-day bog cells (distribution, not only the mean) and
   the share of mapped bog cells with calibrated P above 0.5.
3. Calibration check: bin cells by balanced P(bog); compare mean calibrated P to the
   observed mapped-bog share per bin inside the footprint (reliability diagram).

Caveats to state with the result: the shift assumes the prior is the same across the
domain; the model sees only climate/terrain/paleo, so calibrated P also absorbs the
non-climatic controls it cannot see; the EU rows still shape the fit.

Output: `output/pl2/calibration_prior_shift.csv` and a reliability-diagram PNG.

## Step C: percentile loss criterion at present bogs

1. Reference distribution: balanced P(bog) at present-day mapped bog cells in the footprint.
2. For each scenario, the share of present-bog cells whose future balanced P(bog) falls
   below the reference 10th percentile and below the reference minimum
   (minimum-training-presence style), and the area of cells that stay at or above it.
3. Same for calibrated P(bog) relative to the reference median: the share of present bog
   cells retaining at least half of the reference odds.
4. Compare with the argmax persistence already in
   `scenario_comparison_persistence_cells.csv`; report where the three criteria disagree.

Output: `output/pl2/scenario_comparison_criteria.csv`.

## Deliverable and placement

One new script, `R/pl2_calibrationCheck.R`, after `pl2_compareScenarios.R` in RUNALL
(read-only on existing rasters; minutes, not hours). Results would feed one sentence plus
one supplementary table in the manuscript, and decide whether Results keeps quoting the
absolute argmax km2 or switches to "bog-like climate" and relative change.

## Open decisions for JV

- Which criterion is the headline for "stand a chance": class flip (current), percentile at
  present bogs, or calibrated P?
- Is the footprint prevalence to be treated as the Norwegian prior, given the survey's
  partly-covered cells and the otherpeat definition (AR50 vs EPM)?

## Results 2026-10-07

Calibration check run: the argmax bog map is ~137x too generous, and class flip is the most
lenient loss criterion. Run once with `Rscript --vanilla` (read-only on the projections, ~5 min).
Outputs: `output/pl2/calibration_footprint.csv`, `calibration_prior_shift.csv`,
`calibration_reliability_bins.csv`, `calibration_reliability.png`,
`scenario_comparison_criteria.csv`, `scenario_comparison_criteria_agreement.csv`.

**Departure from the plan.** Every Norwegian bog cell (and 56k footprint absence cells) is a
production training row, so the projected current P there is in-sample. Current-climate values
at those cells are replaced by the production fit's OOB probabilities; the plain projection is
reported alongside. It matters: at the area-matched cut, precision/recall is 0.41/0.43 OOB vs
0.77/0.80 in-sample. OOB is still optimistic, because neighbouring bog cells are training rows
too (a quarter of mapped bog cells have OOB balanced P(bog) = 1). Two extra outputs beyond the
plan: the reliability bins as their own CSV, and a criteria-agreement table.

**A. Commission** (footprint, 1,563,694 cells with a prediction; 1,137 mapped bog, prevalence
0.073%). Argmax calls 155,995 cells bog: 137x the mapped bog area, precision 0.7%, recall 0.93
(OOB). Of the calls, 136k fall on artype_60 nonpeat, 2.7k on other peat, the rest on mixed
cells. The land-use screen barely changes this (108x). Area-matched cut on balanced P(bog): 0.979.

**B. Prior shift** (priors nonpeat 0.955, otherpeat 0.045, bog 0.00073; the non-bog remainder is
split by the artype_60-labelled counts, since mixed/unlabelled cells have no model class).
Calibrated P(bog) at mapped bogs: median 0.115, 41% above 0.5. At argmax-bog cells: median
0.0026. Reliability: the prior-shift form holds within ~1.5x for balanced P > 0.2 and overstates
the observed share 2-7x below that. Mean calibrated P over the footprint is 0.00104 against an
observed 0.00073.

**C. Loss criteria.** Reference = OOB balanced P(bog) at mapped bogs: P10 = 0.556, minimum = 0,
**so the minimum-training-presence criterion is vacuous (every cell passes) and should be
dropped.** Share of present bog still suitable, ssp126 / ssp370:

| criterion | mapped bog (footprint) | argmax bog today (Norway) |
|---|---|---|
| class flip | 13.6% / 0 | 12.8% (1,530 km2) / 0 |
| balanced P >= ref P10 | 0.09% (1 cell) / 0 | 0.005% (10 cells) / 0 |
| calibrated odds >= half ref median | 0 / 0 | 0 / 0 |

The class-flip rows reproduce `pl2_compareScenarios.R` exactly. Every cell the P10 criterion
keeps, class flip also keeps: **class flip is the most generous criterion**, because argmax can
call bog at P(bog) ~0.34-0.55, below where most mapped bogs sit today. Under ssp370 every
criterion loses everything.

**Expected calibrated bog area (Norway):** 121 km2 now, 75 (ssp370), 47 (ssp126). The ordering
is inverted (the milder scenario loses more), because the sum is dominated by low-P mass spread
over millions of cells, exactly where the reliability diagram shows the prior shift overstating.
Do not quote it.

**Bearing on the open decisions.** (1) Headline criterion: the ssp370 conclusion is the same under
all three; for ssp126 the "stand a chance" number depends entirely on the criterion (12.8% vs
~0%). (2) The footprint prior looks usable as the Norwegian prior for balanced P > 0.2, not for
domain-wide sums.

## Decision 2026-10-07 (JV)

Headline loss criterion = class flip (#1), mainly for consistency with the rest of the
pipeline (CV skill, reliability curves and transitions all use argmax). Check and probably
report the P10 criterion (#2) as the strict bound, so ssp126 persistence is a bracket between
the two. Calibrated odds (#3) is not a headline; it only explains why absolute argmax area is
not bog area. Reasons for and against each are in `notebook.md` 2026-10-07 (c).

## Note 2026-10-07: equal priors hold, but bog OOB rests on ~15 trees

The training frame is 25:25:1, but `make.wt` (equal total weight per class) plus `make.size`
(3 x n_bog per tree, `swor`) give near-equal in-bag shares: 33.7 / 33.5 / 32.8% (measured with
rfsrc's sampler on the real class counts). So the prior shift's 1/3 holds. Each tree takes ~98%
of the bog rows, though, so a bog row is OOB in ~1.5% of trees (~15 of 1,000; estimate). The
OOB reference at mapped bogs is coarse (quantiles 3/7, 5/9) and the P10 cut is noisy.
`samptype = "swr"` or a smaller per-class sample would fix it; needs a refit, not done.

## Note 2026-11-19: does the noisy P10 matter, and what a refit would cost

**Sensitivity (read-only on the projections, no refit; throwaway script, not in the
pipeline).** Share of today's argmax-bog cells over Norway (191,888) whose future balanced
P(bog) is at or above a cut:

| cut | 0.30 | 0.35 | 0.40 | 0.45 | 0.50 | 0.556 (P10) | 0.60 |
|---|---|---|---|---|---|---|---|
| SSP1-2.6 (%) | 34.4 | 18.9 | 10.6 | 3.9 | 0.44 | 0.005 | 0 |
| SSP3-7.0 (%) | 0.58 | 0.03 | 0.002 | 0 | 0 | 0 | 0 |

Under SSP1-2.6 the largest future P(bog) at any of these cells is 0.577 (99th percentile 0.488).
So SSP3-7.0 is ~0 at every cut, and the SSP1-2.6 strict bound is ~0 for any reference cut above
0.5. It moves only if the stable P10 lands in 0.40-0.50 (it would be 10.6% at 0.40, about the
class-flip 12.8%). The current reference quantiles are q5 = 0.429 (3/7) and q10 = 0.556 (5/9),
so the 10th percentile is computed from about 9 OOB trees, and that range is plausible.

**Three ways to stabilise the reference, and what each touches.**

1. *Sidecar fit for the reference only (recommended).* Same call as `pl2_predict.R` (same
   `make.wt` / `make.size` / `swor`, same frame) with ntree ~13,000 (bog rows then OOB in ~200
   trees), `importance = FALSE`, `forest = FALSE` (only `predicted.oob` is needed, which keeps
   memory down), inside `pl2_calibrationCheck.R`. It is the same estimator as the production
   forest, so its OOB P at bog rows is what the 1,000-tree forest would give with enough trees.
   It touches only `pl2_calibrationCheck.R` and its outputs, plus the P10 rows in the ms. The
   production model, maps, transitions, CV, ruler and AOA are unchanged. Cost: one fit, est.
   tens of minutes.
2. *Refit production with `samptype = "swr"`* (or a per-class sample of ~0.63 n_bog, `swor`).
   Bog rows OOB in ~37% of trees. The in-bag prior stays equal (1/3 holds), but the learner
   changes: each tree sees ~63% of the distinct bog rows, so the projection changes too.
   Ripple: `pl2_predict.R` (~1.5 h) -> `pl2_interpret.R`, `pl2_compareScenarios.R`,
   `pl2_calibrationCheck.R`, `ms/make_figures.R` (projection and scenario figures, domain
   shares, VI figure) and every projection number in the ms (transitions, 104 polygons,
   areas, S6 Table). For the CV skill and the ruler to describe the same learner, the recipe
   would also have to change in `pl2_weightFeaturesDataPartitioning.R` (-> `pl2_freezeDIRuler.R`
   -> DI / AOA) and `pl2_evaluate.R` (-> `pl2_fitErrorProfiles.R` -> `pl2_mapReliability.R`).
   In effect all of pl2 except the modelling frame, the occupancy grid and the partitions
   (pinned to bio10). Not worth it for one reference number.
3. *Spatially held-out reference* (e.g. leave-one-Lyngstad-region-out with the production
   recipe). Also removes the neighbour optimism that neither 1 nor 2 touches (a quarter of
   mapped bogs have OOB P = 1 because neighbouring bog cells are in-bag). But the result is the
   score of models that never saw the region, which is lower and partly an extrapolation score.
   It would be a different quantity from the projection it thresholds. Possible as a
   sensitivity row, not as the reference.

Not affected by any option: pl0, the modelling frame, the partitions, the model-free occupancy
grid.

## Results 2026-11-19: sidecar reference fit (option 1, JV)

`pl2_calibrationCheck.R` now fits a sidecar forest: the production recipe with 13,000 trees,
`forest = FALSE`, seed -4242. Its OOB values replace the production OOB at the Norwegian
training cells. Run through a one-script copy of RUNALL (manifest: ok, 6.9 min; the 13k-tree
fit itself takes a few minutes). All upstream summary files are unchanged. New outputs:
`calibration_oob_trees.csv`, `calibration_reference.csv`,
`scenario_comparison_criteria_sweep.csv`.

- **OOB rate measured** (200-tree in-bag record): a bog row is OOB in 1.55% of trees, i.e.
  15.5 of the production forest's 1,000 and ~201 of the sidecar's 13,000. Absence rows: 96%.
- **Reference at the 1,137 mapped bogs:** P10 = **0.568** (was 0.556). The coarse value was
  close. Share at OOB P = 1: 9.3% (was 41%, which was a few-tree artefact).
- **Strict criterion unchanged in substance.** SSP1-2.6: 1 of 191,888 argmax-bog cells
  (0.0005%), 0 of the mapped bogs; SSP3-7.0: 0. Class flip: 12.8% (1,530 km2) and 13.6% of
  mapped bogs under SSP1-2.6, 0 under SSP3-7.0 (reproduces `pl2_compareScenarios.R`).
- **Cut sweep, SSP1-2.6** (argmax bog today / mapped bog): 0.40 -> 10.6% / 13.5%;
  0.45 -> 3.9% / 8.5%; 0.50 -> 0.4% / 2.1%; >= 0.55 -> ~0. Max future P(bog) 0.577 / 0.559.
  SSP3-7.0 max 0.41: ~0 at every cut. The stable P10 is well above the steep 0.40-0.50 range,
  so the strict bound of ~0 is robust.
- **What moved, because the current-climate OOB values are now less noisy:**
  - Commission: argmax 156,025 cells, 137x mapped area, precision 0.70%, recall 0.956.
  - Area-matched cut 0.976: precision/recall 0.29 (was 0.41/0.43).
  - Calibrated P(bog) at mapped bogs: median 0.030 (was 0.115); share > 0.5 17% (was 41%).
  - Reliability: the prior shift now tracks the observed share within ~1.5x for balanced P
    > 0.3 (top bin 0.080 vs 0.087 observed). Below 0.3 it still overstates, by 5-10x.
  - Calibrated expected area (still inverted, still not to be quoted): 103 / 75 / 47 km2.
- Minimum-presence cut still vacuous (reference minimum 0.007, every cell passes): drop it.

## Manuscript placement (prepared 2026-11-19, not yet written)

The branch is merged with main (098c629), so `ms/ms.qmd` is current. Numbers to be read in the
setup chunk from `output/pl2/scenario_comparison_criteria.csv` and
`output/pl2/calibration_footprint.csv` (as `sc_cells` already is), never typed.

1. **Methods, "The model"** (the paragraph on the rebalanced prior): add one sentence with the
   footprint check. Inside the survey footprint, where raised bog is mapped completely
   (prevalence 0.07%), the cells predicted raised bog cover ~137x the mapped area (precision
   0.7%, recall 0.93, OOB at training cells). Hence argmax area is an area of bog-like
   climate, not bog extent.
2. **Methods, "Projection and reliability mapping"**: after "To compare the scenarios we
   followed ...", define the strict criterion. A cell stays suitable if its future raised-bog
   score is at least the 10th percentile of the OOB score at the mapped raised bogs today (a
   percentile training-presence threshold), alongside the class-flip criterion. Cite Pearson
   et al. 2007 / Radosavljevic & Anderson 2014 (verify; not yet in the .bib).
3. **Results, "A milder scenario"**: after the polygon sentence, the bracket. Of today's
   raised-bog cells, 12.8% stay raised bog by class under SSP1-2.6, and ~0% stay above the
   strict cut. Under SSP3-7.0 none do by either criterion. Add one clause on cut sensitivity
   (the SSP1-2.6 strict bound is ~0 for any cut above 0.5 and reaches the class-flip value only
   near 0.40).
4. **S6 Table caption**: replace "over-calls raised bog relative to its prevalence" with the
   footprint ratio. Add a row "Of today's raised-bog cells, score >= P10 of mapped bogs (km²)".
5. **New S7 Table**: criteria x scenario for the two populations (mapped bog in the footprint,
   argmax bog over Norway). Rows: class flip and P10, plus the cut sweep. The minimum-presence
   row and the expected-area rows are dropped (vacuous / inverted, see Results 2026-10-07).
   Built by `ms/make_figures.R` into `ms/tables/`. The cut sweep moves from the throwaway
   script into `pl2_calibrationCheck.R`.
6. **Setup-chunk guards**: the class-flip rows equal `scenario_comparison_persistence_cells.csv`;
   P10 keeps no cell that class flip drops; no criterion keeps a cell under SSP3-7.0 (the
   prose says so).
7. **Discussion TODO (v)** and the abstract skeleton: "still stand a chance" becomes a bracket
   (12.8% by class, ~0% by the strict cut, of today's bog-like cells; 19% of polygons by class).
8. **Before rendering**: `pl2_calibrationCheck.R` has no entry in `output/run_manifest.csv`
   (it was run by hand), so rerun it through a sliced RUNALL after the edits above. That
   records it, refreshes its outputs, and lets the setup-chunk freshness guard cover it. Then
   rerun `ms/make_figures.R`.

Depends on the refit decision: only the P10 value and the P10 rows (items 2-5, 7). Nothing
else in the list changes.
