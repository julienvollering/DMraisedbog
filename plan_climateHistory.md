# Plan: the climate history of today's raised bogs against the future climate

Status: plan written 2026-10-06 on branch `climate-history-plan`, revised after review; steps 1-3 of
section 7 implemented and run the same day (first results in section 0 below).

Departures from the plan in the implementation:
- **No edit to `pl0_downloadCHELSA.R`.** `R/pl0_prepareCHELSATraceSeries.R` reads the Norway window of
  each centennial slice over `/vsicurl/` (~4 s per slice) instead of downloading ~56 GB of global
  files. It must not write into `data/CHELSA/paleo/`: `pl0_prepareCHELSATrace.R` globs every bio01/bio12
  file there and would either stop on mismatched slices or silently change `paleo_*`.
- **Future climate is read from native CHELSA at each pixel**, for both scenarios, rather than from the
  modelling frame, so SSP1-2.6 works before its projection leg has run.
- **Definition 1b is one-sided on moisture** (not drier than the pixel's driest slice): the future is
  wetter than every slice of every pixel's past, the safe direction for an ombrotrophic bog.
- **The GCM-cell bootstrap is weak**: 614 pixels fall in 6 approximate CCSM3 cells, 510 of them in one.
  Bands are reported but are not a confidence statement.
- **Not done yet**: ecospat/hypervolume (step 4; `ecospat` pulls in retired `biomod2` and `cluster` -
  COUE could instead be computed directly on the KDE grid), the multidecadal-variability widening (3.8),
  core-cell sensitivity and the literature t0 band (step 5), `ms/make_figures.R` (step 6).

## 0. First results (2026-10-06, delta, bio10 x bio12)

Outputs: `output/pl2/climatehistory_*.csv|png`. 1,116 NO bog cells -> 614 TraCE pixels in 6 approximate
CCSM3 cells (510 in one), so effective n is ~1-6.

| | t0 = 0 | 6 ka | 10 ka |
|---|---|---|---|
| SSP3-7.0, pooled hull | 0.2% | 3.9% | 3.9% |
| SSP3-7.0, per-pixel headroom (median) | +3.5 C | +2.6 C | +2.6 C |
| SSP1-2.6, pooled hull | 37% | 99% (t\* = 6 ka) | 99% |
| SSP1-2.6, per-pixel coverage / median headroom | 0% / +1.2 C | 11% / +0.3 C | 11% / +0.3 C |

Reference, today's pooled NO+EU bog climate (hull): 21% (SSP3-7.0), 99% (SSP1-2.6).

- TraCE's simulated HTM is weak here: max 8-6 ka bio10 only +0.7 C over the last slice (mean +0.1 C).
  The SSP1-2.6 verdict (+0.3 C headroom) sits inside any plausible proxy-vs-TraCE bias; the SSP3-7.0
  one (+2.6 C) does not.
- TraCE last slice vs CHELSA 1981-2010: bio10 +0.3 C, bio12 -18%, bio18 -22%. The "near zero by
  construction" expectation (3.3) holds for temperature only.
- Moisture never binds (future wetter than every past slice), so temperature carries the verdict.

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

- **Coverage curve**: share of future bog pixels that fall inside the experienced climate
  space, for every assumed transition age t0 from 0 (today's climate of the Norwegian bogs)
  to ~9 ka. Note t0 = 0 is *not* `fig_occupancy`'s reference: that figure compares against
  all training bogs, Norwegian and European. Draw the pooled NO+EU present as a reference
  line so the two figures can be read together.
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
3. **500-yr spacing samples one centennial mean in five**, so a per-pixel maximum can miss
   the warmest HTM century and bias the verdict toward "novel". Use every centennial slice
   for 10-0 ka (~100 files per variable); read only the Norway window over `/vsicurl/`
   rather than storing global files.
4. Neither `ecospat`, `hypervolume` nor `ENMTools` is installed (checked).

## 3. Design

### 3.1 Axes
Primary: **bio10 x bio12** (2-D, readable, both available in TraCE and CHELSA v2/future).
Sensitivity: add bio06 (coldest-month; winter severity) as a third axis, and swap bio12 for
bio18. Keep the axis count at 2-3: the occupancy figure shows that the shifting predictors
are all collinear temperature measures, so more axes add dimensionality without adding
information, and hull/density methods degrade quickly.

### 3.2 Unit of analysis: TraCE pixels, resampled by GCM cell
TraCE is ~1 km. 1,116 Norwegian bog cells at 250 m collapse to a few dozen to a few hundred
distinct ~1 km pixels; extract climate and report coverage per pixel. But pixels are not
independent: after the delta correction (3.3) each pixel's history is its local baseline
plus a CCSM3 anomaly from a ~3.75 degree grid, and the GFDL future anomaly is similarly
smooth. All of Norway sits in a handful of GCM cells, so the effective n is that handful.
Bootstrap by GCM cell; a pixel or polygon bootstrap would give far too narrow intervals.

### 3.3 Experienced climate, with a TraCE-to-CHELSA bridge
For every bog pixel i and slice s: c_i(s) = (bio10, bio12). Two ways to place TraCE on the
CHELSA v2 scale, run both:

- **Delta (primary)**: bio10 additive, CHELSA_1981-2010(i) + [TraCE(i, s) - TraCE(i, 0 ka)];
  precipitation as a **ratio**, CHELSA_1981-2010(i) x TraCE(i, s) / TraCE(i, 0 ka) (an
  additive precipitation delta can go negative at dry sites). Puts past and future on one
  scale (CHELSA v2 current vs `FUTURE_GCM` under `FUTURE_SCENARIOS`, `R/config.R`).
- **Raw TraCE** as sensitivity.

Baseline mismatch: TraCE ends in 1990, so its last centennial slice is mostly 20th-century
climate, cooler than CHELSA 1981-2010. Anchoring the anomaly to 1981-2010 then makes the
past too warm by that difference (bias toward "inside"). Confirm what period `time = 20`
covers, quantify the difference, and state it. The present-day TraCE-vs-CHELSA offset is
near zero by construction (CHELSA-TraCE is itself delta-downscaled onto a CHELSA
climatology), so it is a plumbing check only and says nothing about HTM fidelity; that
needs the proxy comparison in 3.8.

### 3.4 Per-pixel window, capped by deglaciation
For a given t0, pixel i's history is slices with age <= min(t0, `paleo_years_icefreeland`_i).
A bog cannot be older than its ground, so this cap is automatic and removes the
glaciated or inundated part of the series. Caveat: TraCE `orog` palaeotopography (check whether it
includes ICE-5G/6G isostasy and sea level) is too coarse to resolve shoreline displacement
along fjords and coast, so sites that emerged from the sea late may be treated as land too
early; the cap is a *maximum* age. Flag these pixels (low elevation, near coast) rather than hide them.

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
   contour. Robust to a few extreme slices that a hull would lock in, but that is also its
   weakness here: the question is about extremes, and a contour trims rare warm periods
   like the HTM by design. It is not monotone in t0 either (adding cool early slices can
   push HTM points outside), so t\* is defined on definitions 1-2 only; 3 is a sensitivity.

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
  Schoener's D. Two pitfalls: (a) `ecospat.niche.dyn.index()` can restrict the indices to
  climates present in both backgrounds, which would drop exactly the novel future climate
  we ask about; set `intersection = NA` and say so. (b) Its availability correction
  (background = all Norwegian land at that slice) answers "do bogs prefer this climate",
  not containment for fixed sites, and D measures overlap; so COUE expansion will not
  equal the coverage curve, and is reported alongside it, not as a check of it. Run on a
  few t0 (0, 2, 4, 6, 8 ka) and per slice, not on every t0.
- **hypervolume** (Blonder): sensitivity for the 3-axis version and for set operations
  (`hypervolume_set`: unique/intersection fractions). Cost and bandwidth sensitivity make
  it a check, not the base.
- **ENMTools**: identity/background tests answer "is the niche *different* from chance",
  which is not our question (we want containment of a future in a history, not a test of
  equivalence). Leave out unless a reviewer asks; do not install for the base analysis.
- The transparent **headroom + hull + density-contour** set (3.5) carries the verdict;
  hypervolume checks it, ecospat puts it in the niche-dynamics vocabulary reviewers expect.

### 3.7 Bog-age information to narrow t0 (this is the tricky part)
Do not pick one t0. Provide:

- **A sweep** t0 = 0, 0.5, 1, 2, ..., 9 ka (slice spacing is 0.5 ka).
- **A weak prior from the literature**, drawn as a shaded band on the coverage curve, never
  as a single line. Edvardsson et al. (check year: 2012/2016) is the anchor: bog-pine chronologies in *Sphagnum peat
  over fen peat* at three bogs in southern Sweden span 5284-3728 BC and 2668-1108 BC
  (~7.2-5.7 ka and ~4.6-3.1 ka BP, BP = 1950; TraCE ages here are b2k, a 50-yr
  difference to state once), with bog-tree establishment during the relatively warm,
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
  detail is not resolved, and the HTM amplitude is model-dependent; CCSM3/TraCE likely
  shows weaker Fennoscandian HTM summer warmth than proxies (bias toward "novel"). Validate bio10 at a
  handful of sites against Fennoscandian proxy summer-temperature reconstructions, or at
  least report the simulated HTM-minus-present anomaly next to published ones.
- **Smoothing mismatch**: TraCE slices are centennial means; the future layer is a 30-yr
  mean. A smoother past makes the future look *more* extreme than the past (the past's
  warm years are averaged away). Bogs survived the warm *years* of the HTM, which TraCE
  cannot show. Mitigation: widen the envelope by the SD of 30-yr means about a centennial
  mean (multidecadal variability, e.g. from TraCE's own variability or a long control
  run), **not** the interannual SD, which is far larger and would overcorrect toward
  "inside"; state the direction of the bias.
- **Not the same climate, even inside the envelope**: future change is ~100x faster than
  Holocene change; CO2 is ~260 ppm then vs >500 ppm under SSP3-7.0 (plant water use); HTM
  warmth came with colder winters, so the bio06 sensitivity will flag "novel" from winter
  warming almost by construction (say so before reading it); and surviving the HTM is not
  surviving unchanged (the Edvardsson bog pines record dry surfaces).
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
| 1 | `R/pl0_downloadCHELSA.R` (edit) | add `bio10`, `bio18` (and `bio06` for the 3-axis check), plus centennial slices for 10-0 ka (gap 3; keep the existing 500-yr set for `paleo_*`); the loop already skips existing files. Check integrity as with the other downloads. Changes the `record_settings()` row; note it in `notebook.md`. |
| 2 | `R/pl0_prepareCHELSATraceSeries.R` (new), listed in `RUNALL.R` after `pl0_prepareCHELSATrace.R` | crop the new variables to the Norway window, **keep the time dimension**: one multi-layer GeoTIFF per variable (`data/CHELSA/paleo_derived/trace_series_<var>_Norway.tif`). Do **not** edit `pl0_prepareCHELSATrace.R`: its `paleo_*` outputs feed the model and must not change. |
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
   - **A.** bio10 x bio12 small multiples, one panel per ~2 ka slice (0, 2, 4, 6, 8 ka): bog
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

- **Run first, decide after** (cheap): what period the `time = 20` slice averages (the
  scripts treat it as present; it sets the baseline mismatch in 3.3), and the present-day
  offset as a plumbing check. The slice list contains `0` twice (the two `seq()` calls
  overlap); harmless, since the download skips existing files.
- **Which future?** `FUTURE_SCENARIOS` in `R/config.R` (SSP3-7.0 primary, SSP1-2.6 mild,
  GFDL-ESM4): loop over both. The milder scenario is where t* is most likely to be small, so it matters for the
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

1. Download bio10/bio18 (centennial, 10-0 ka); crop to a time-series stack; check the
   `time = 20` period and the baseline mismatch.
2. Headroom and coverage curve for the delta version, bio10 only (fastest read of the answer).
3. Add bio12, the hull and density definitions, pixel-level bootstrap.
4. ecospat COUE; hypervolume sensitivity.
5. Literature pass for the t0 band; core-cell sensitivity.
6. Figures into `ms/make_figures.R`; text for Methods, Results, and the Discussion caveats
   in 3.8.
