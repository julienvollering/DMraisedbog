# Plan: the climate history of today's raised bogs against the future climate

Status: plan only, nothing run. Written 2026-10-06 on branch `worktree-ms-refresh-from-rerun`.

## 1. The question, and what would answer it

> Is the expected future climate of today's Norwegian raised bogs **outside or inside the
> climate space those bogs have experienced since their fen-to-bog transition?**

`fig_occupancy` answers a *spatial* version of this: it compares the future with climates
bogs occupy **today** (1981-2010), across space. Raised bogs are millennial structures, so
the relevant reference is the climate each bog has **persisted through**. A future that is
novel against today's niche may be well inside the Holocene thermal maximum (HTM, ca. 8-6
ka) envelope, and then "novelty" overstates the risk of loss.

The analysis has to give a verdict that is explicit about its one unknown. That unknown is
**t0, the age at which a bog became a bog**, so the headline result is a function of t0,
not a number:

- **Coverage curve**: share of future bog cells that fall inside the experienced climate
  space, for every assumed transition age t0 from 0 (today's niche, the current figure) to
  ~9 ka.
- **t\***: the smallest t0 at which coverage reaches a stated threshold (e.g. 95%), or
  "never". Read t\* against the independent evidence on when bogs ombrotrophied.

If coverage is already high at t0 = 2-3 ka, the conclusion is robust to our ignorance. If it
only reaches the threshold when the HTM is included, the conclusion rests on bogs being
>= 6-8 ka old, and the manuscript should say exactly that.

## 2. What we already have, and what is missing

Already in the repo:
- `R/pl0_downloadCHELSA.R`: TraCE21k at 500-yr spacing, -140..20 (hundred-year units; 20 =
  2000 CE, so BP = (20 - t) x 100; the HTM, 8-6 ka BP, is 5 slices), but **only bio01,
  bio12, glz, orog**.
- `R/pl0_prepareCHELSATrace.R`: reduces the series to 8 time-**averaged** `paleo_*`
  predictors plus years since deglaciation. These are the model's static paleo axes; the
  time dimension is collapsed there and is exactly what we now need back.
- `R/pl2_exploreOccupancy.R` + `ms/make_figures.R` (Fig 2): bio10 x gsp occupancy.
- `future_bog` construction: Norwegian bog cells carried to the future scenario
  (persistence comparison, same cells).

Gaps, found while planning:
1. **bio01 is the wrong temperature axis.** The occupancy figure and the partition axis are
   bio10 (warmest quarter). TraCE's HTM signal at high latitude is mainly a *summer*
   (orbital) warming; annual means are muted by cooler winters. Using bio01 would
   understate the early-Holocene warmth and bias the verdict toward "future is novel".
   Verified 2026-10-06 by HTTP HEAD: `bio10`, `bio08`, `bio18`, `bio04`-`bio06`, `scd` exist on
   the TraCE server at the same file pattern; `gdd*`, `gsl`, `gst`, `npp` do not.
2. **gsp does not exist in TraCE.** The moisture axis must be bio12 (or bio18, warmest-quarter
   precipitation). Fig 2 already has a bio12 robustness grid. Re-draw the *present-day*
   occupancy on the same axes so the two figures are comparable.
3. Neither `ecospat`, `hypervolume` nor `ENMTools` is installed (checked).

## 3. Design

### 3.1 Axes
Primary: **bio10 x bio12** (2-D, readable, both available in TraCE and CHELSA v2/future).
Sensitivity: add bio06 (coldest-month; winter severity) as a third axis, and swap bio12 for
bio18. Keep the axis count at 2-3: the occupancy figure shows that the shifting predictors
are all collinear temperature measures, so more axes add dimensionality without adding
information, and hull/density methods degrade quickly.

### 3.2 Unit of analysis: TraCE pixels, not 250 m cells
TraCE is ~1 km. 1,116 Norwegian bog cells at 250 m collapse to a few dozen to a few hundred
distinct ~1 km pixels. Extract climate per pixel and weight or resample by pixel, so n is
honest. (The same autocorrelation caution `pl2_exploreOccupancy.R` already makes for the warm
evidence.) Bootstrap by pixel, or by Lyngstad polygon where available.

### 3.3 Experienced climate, with a TraCE-to-CHELSA bridge
For every bog pixel i and slice s: c_i(s) = (bio10, bio12). Two ways to place TraCE on the
CHELSA v2 scale, run both:

- **Delta (primary)**: c_i(s) = CHELSA_1981-2010(i) + [TraCE(i, s) - TraCE(i, 0 ka)]. Removes
  the TraCE-vs-CHELSA offset at the present, so past and future are on one scale (CHELSA
  v2 current vs GFDL-ESM4 SSP3-7.0, `FUTURE_SCENARIO`/`FUTURE_GCM` in `R/config.R`).
- **Raw TraCE** as sensitivity.
First report the present-day offset (TraCE 0 ka vs CHELSA) per axis: if it is large, that
alone is a result about how far the HTM comparison can be trusted.

### 3.4 Per-pixel window, capped by deglaciation
For a given t0, pixel i's history is slices with age <= min(t0, `paleo_years_icefreeland`_i).
A bog cannot be older than its ground, so this cap is automatic and removes the
glaciated or inundated part of the series. Caveat: TraCE `orog`/`glz` carry no isostatic
uplift, so coastal sites that emerged from the sea late are treated as land too early; the
cap is a *maximum* age. Flag these pixels (low elevation, near coast) rather than hide them.

### 3.5 "Inside the experienced space": three definitions, reported together
For each t0, with E(t0) = the set of (bio10, bio12) points experienced by *any* bog pixel
in its own window, and F = the future points of the same pixels:

1. **Per-pixel (local history)**: is F_i inside the range of *pixel i's own* trajectory?
   Reported as **headroom** in degrees C: future bio10 minus the historical maximum at that
   pixel. The most direct reading of "this bog has been this warm before" and it needs no
   multivariate machinery. Map it.
2. **Pooled hull**: convex hull of E(t0) in bio10 x bio12 (alpha-hull as variant), share of F
   inside. Answers "has *a* bog been through this climate", the niche-style question.
3. **Pooled density quantile**: kernel density of E(t0), F inside the 95% (and 99%)
   contour. Robust to a few extreme slices that a hull would lock in.

Definitions 1 and 2 differ on purpose. If F is in the pooled space but outside its own
history, bogs are surviving a future that only *other* bogs have seen; say which one the
conclusion rests on. Same cells, never a resample, as with `future_bog` today.

### 3.6 Where ecospat / hypervolume / ENMTools fit
Use each for what it does well, not as the engine:

- **ecospat (COUE) as the headline niche figure and table.** Niche A = experienced space E(t0)
  (or one time slice s), niche B = future bog climate. `ecospat.grid.clim.dyn()` then
  `ecospat.niche.dyn.index()` return exactly the quantities in the question:
  **expansion** (future climate outside the experienced niche = the novelty we are asking
  about), **stability** (inside), **unfilling** (experienced but not in the future), plus
  Schoener's D. Pass the *available* climate (all Norwegian land cells at that slice) as the
  background so density is corrected for availability; this matters because the HTM
  landscape was differently distributed in climate space. Run on a few t0 (0, 2, 4, 6, 8 ka)
  and per slice (e.g. 0, 2, 4, 6, 8 ka), not on every t0.
- **hypervolume** (Blonder): sensitivity for the 3-axis version and for set operations
  (`hypervolume_set`: unique/intersection fractions). Cost and bandwidth sensitivity make
  it a check, not the base.
- **ENMTools**: identity/background tests answer "is the niche *different* from chance",
  which is not our question (we want containment of a future in a history, not a test of
  equivalence). Leave out unless a reviewer asks; do not install for the base analysis.
- The transparent **headroom + hull + density-contour** set (3.5) carries the verdict;
  ecospat and hypervolume confirm it. Agreement across the methods is the robustness claim.

### 3.7 Bog-age information to narrow t0 (this is the tricky part)
Do not pick one t0. Provide:

- **A sweep** t0 = 0, 0.5, 1, 2, ..., 9 ka (slice spacing is 0.5 ka).
- **A weak prior from the literature**, drawn as a shaded band on the coverage curve, never
  as a single line. Edvardsson (2013) is the anchor: bog-pine chronologies in *Sphagnum peat
  over fen peat* at three bogs in southern Sweden span 5284-3728 BC and 2668-1108 BC
  (~7.2-5.7 ka and ~4.6-3.1 ka BP), with bog-tree establishment during the relatively warm,
  dry HTM and lateral expansion as it turned humid. That says some southern Fennoscandian
  bogs were already ombrotrophic well before 6 ka; it does **not** date Norwegian bogs, and
  tree establishment is a lower bound on the transition, not the transition.
  Needs a literature pass before the band is drawn: a Norwegian/west-Nordic basal and
  fen-to-bog transition compilation (e.g. peatland basal-date compilations such as Loisel
  et al. 2014 and successors, plus regional stratigraphy). *Basal peat age is not the
  fen-to-bog age*: it is a bound (the bog is younger than its basal peat), so treat it as
  an upper limit on t0 alongside the deglaciation cap.
- **Within-bog structure**: large domes are old at the centre and young at the margins
  (lateral expansion). The Lyngstad polygons give size; compute distance to polygon edge.
  Run the whole analysis on "core" cells (large polygons, interior) as a sensitivity: these
  are the cells for which a long history is most plausible, and the ones the question is
  about.

### 3.8 Known limits of the climate record (state them in the Methods)
- TraCE21k is one transient simulation, downscaled; GCM-scale (~3.75 degree) forcing means local
  detail is not resolved, and the HTM amplitude is model-dependent. Validate bio10 at a
  handful of sites against Fennoscandian proxy summer-temperature reconstructions, or at
  least report the simulated HTM-minus-present anomaly next to published ones.
- **Smoothing mismatch**: TraCE slices are centennial means; the future layer is a 30-yr
  mean. A smoother past makes the future look *more* extreme than the past (the past's
  warm years are averaged away). Bogs survived the warm *years* of the HTM, which TraCE
  cannot show. Mitigation: compare with the modern interannual SD added (envelope widened by
  the 1981-2010 interannual spread), and state the direction of the bias.
- Absence from the experienced space is a statement about *this* climate record, not about
  bog physiology; hydrology (bio12 vs evaporative demand) is only partly captured by
  precipitation alone. A warm-dry past is not a warm-wet future; the moisture axis is
  the reason the analysis is 2-D, and an effective-moisture proxy (e.g. bio12 against
  temperature-based demand) is a follow-up, not a base requirement.
- Present-day-fixed locations: bogs do not move in this analysis; fine for persistence, but
  it cannot say anything about bogs that existed in warmer-climate places and are now gone.

## 4. Pipeline changes

New and edited scripts (names follow the pl0/pl2 convention; this is **downstream and read-only
to the model**: nothing here touches the frozen DI ruler, the production fit or
the `paleo_*` predictors, so no model rerun is implied):

| Step | Script | What |
|---|---|---|
| 1 | `R/pl0_downloadCHELSA.R` (edit) | add `bio10`, `bio18` (and `bio06` for the 3-axis check) to `paleo_bioclim_vars`; the loop already skips existing files. Check file integrity as with the other downloads. ~33 slices x 2-3 vars. |
| 2 | `R/pl0_prepareCHELSATraceSeries.R` (new) | crop the new variables to the Norway window, **keep the time dimension**: one multi-layer GeoTIFF per variable (`data/CHELSA/paleo_derived/trace_series_<var>_Norway.tif`). Do **not** edit `pl0_prepareCHELSATrace.R`: its `paleo_*` outputs feed the model and must not change. |
| 3 | `R/pl2_exploreClimateHistory.R` (new), listed in `RUNALL.R` right after `pl2_exploreOccupancy.R` | Sections: bog pixels and deglaciation cap; TraCE-vs-CHELSA offset check; delta-corrected series; E(t0) and the three containment definitions; coverage curve; ecospat COUE at selected t0; hypervolume sensitivity. Writes `output/pl2/climatehistory_*.csv` and PNGs; calls `record_settings()` with t0 grid, axes, thresholds. |
| 4 | `ms/make_figures.R` (edit) | new figure from the CSV outputs only. |
| 5 | `plan` / `notebook.md` | record the rerun and the finding. |

Packages: add `ecospat` (and `hypervolume` for the sensitivity) to CLAUDE.md "Key
Dependencies" once installed; check that the `ecospat` install does not drag in retired
spatial packages.

Reuse: `future_bog` construction and `record_settings()`; the pixel-level `extract` pattern in
`extract_native_predictors()` (`R/functions.R`); Lyngstad polygons for polygon size and
edge distance.

## 5. Figures and tables

1. **Space-time figure (the "third dimension")**, 7.5 in wide, shown like `fig_occupancy`:
   - **A.** bio10 x bio12 small multiples, one panel per ~1 ka slice (0, 2, 4, 6, 8 ka): bog
     pixels' climate for that slice (points + hull/95% contour), future bog climate in red
     in every panel. The reader sees the red cloud move from outside the contour to inside
     as the panels go back to the HTM, or not.
   - **B.** Same data as a *time axis*: bio10 against age (ka BP), per-pixel trajectories
     as light lines, median and 5-95% band, with the future range drawn as a bar at the
     right edge and the Edvardsson-derived t0 band shaded. Where the future bar overlaps
     the history, it reads at once.
2. **Coverage curve**: % future bog pixels inside E(t0) against t0, for definitions
   1-3, delta and raw, with the t0 prior band. This is the headline.
3. **Headroom map** (definition 1), Norway, one value per pixel at a chosen t0.
4. **Table**: COUE (expansion / stability / unfilling, Schoener's D) at t0 = 0, 2, 4, 6, 8 ka,
   plus per-definition coverage at the same t0.

## 6. Decisions and checks before running

- **Run first, decide after** (cheap): the TraCE-vs-CHELSA offset at the present for bio10
  and bio12, and whether 0 ka in TraCE is 2000 CE (script treats `time = 20` as present; the
  slice list also contains `0` twice because the two `seq()` calls overlap at 0, harmless but
  worth confirming).
- **Which future?** Primary scenario is the one in `R/config.R` (SSP3-7.0, GFDL-ESM4). The
  branch `worktree-ssp126-scenario` adds an SSP1-2.6 leg and `FUTURE_SCENARIOS`; when merged,
  loop over both. The milder scenario is where t* is most likely to be small, so it matters for the
  message.
- **Coverage threshold** (95%? any pixel outside?) and the **t0 prior band**: JV's call, and
  needs the Norwegian literature pass in 3.7.
- **Scope of bogs**: Norwegian presence cells only (the same cells as `future_bog`). The
  European bogs would answer a different question (history of the warm-edge analogues) and
  are out of scope for the base analysis, but they are a natural follow-up since the
  warm-evidence concern in `pl2_exploreOccupancy.R` is about exactly those sites.
- Shared-state rules: check for running `Rscript` before the download/crop steps; the cropped
  TraCE cache folder is recreated by `pl0_prepareCHELSATrace.R` and deleted at its end
  (`unlink(temp_dir)`), so the new script must not rely on it.

## 7. Suggested order of work

1. Download bio10/bio18; crop to a time-series stack; check offset to CHELSA and slice timing.
2. Headroom and coverage curve for the delta version, bio10 only (fastest read of the answer).
3. Add bio12, the hull and density definitions, pixel-level bootstrap.
4. ecospat COUE; hypervolume sensitivity.
5. Literature pass for the t0 band; core-cell sensitivity.
6. Figures into `ms/make_figures.R`; text for Methods, Results, and the Discussion caveats
   in 3.8.
