# Plan: how generous is the bog map, and what counts as "no longer suitable"?

Status: plan only, nothing run (2026-10-07). Background: `notebook.md` 2026-10-07 (b).

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
