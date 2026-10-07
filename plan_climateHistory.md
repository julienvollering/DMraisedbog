# Plan: the climate history of today's raised bogs against the future climate

Status (2026-10-07): **analysis done, figure approved by JV, manuscript revision in progress** (section 8).
- Done: pipeline (3 new pl0 scripts, `R/pl2_exploreClimateHistory.R`, all in RUNALL); two records (TraCE21k,
  EPOCH-2) on the 1981-2010 baseline with Berkeley Earth offsets; literature t0 band (3.7a, 2-6 ka core,
  1-8 ka outer); `fig_climatehistory` and `tbl_climatehistory` in `ms/make_figures.R`; current numbers in 0.
- In progress: `ms/ms.qmd` (plan in section 8). The branch's `ms.qmd` was rebased by hand on `main`'s
  (which added SSP1-2.6, f2541e2), because merging `main` into the branch was blocked in the agent
  session. **JV: merge `main` into this branch before rendering**; `ms/make_figures.R` will conflict
  (main's Fig 7 / S6 Table block and this branch's climate-history block are both new sections; keep both),
  then rerun `ms/make_figures.R`.
- Not done: core-cell sensitivity (step 5; proposed to drop unless a reviewer asks, JV to decide: the t0
  sweep already varies the window that a core/margin split would vary); `notebook.md` entry (text handed to JV, the agent
  cannot write the main checkout's notebook); Väliranta 2017 and MacDonald 2015 unverified, not cited.
- Dropped (JV): ecospat / hypervolume / COUE (3.6); correcting for multidecadal variability (3.8, caveat only);
  any blending of TraCE and EPOCH-2 (3.9); the two supplementary figures (slices, exceedance map).

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
- **Decided 2026-10-07 (JV)**: no ecospat / hypervolume / COUE (step 4 dropped; 3.6 is superseded).
  The multidecadal-variability mismatch (3.8) is stated as a Methods caveat, not corrected for.
- **Literature t0 band**: 3.7a, revised after JV's check to 2-6 ka (outer 1-8 ka).
- **Figures/tables (step 6)**: `ms/make_figures.R`, section after Fig 2: `fig_climatehistory`,
  `tbl_climatehistory`, `ch_*` study numbers (the supplementary figures were dropped).
- **TraCE's weak HTM**: 3.9, EPOCH-2 run in parallel (JV 2026-10-07); baseline offsets from Berkeley Earth (decision 1, JV), implemented; results in 0.
- **18 of 614 pixels have no ice-free date** (`icefree_ka` = 0; coastal Helgeland, GCM cell 4_17, TraCE
  palaeogeography): their window is the youngest slice only, at every t0. Stated in S2 Text, not changed.

## 0. Current results (2026-10-07: two records, both on the 1981-2010 baseline)

Outputs: `output/pl2/climatehistory_*.csv`; figure `ms/figures/fig_climatehistory.png`, table
`ms/tables/tbl_climatehistory.csv`. 614 bog pixels in 6 approximate CCSM3 cells (510 in one).
Records: TraCE21k (delta) and EPOCH-2 (@mauriClimateEuropeHolocene2015), each placed on CHELSA
1981-2010 with the observed JJA warming from its own present (Berkeley Earth; median at the
pixels D_trace = +0.41 C for 1900-1990, D_epoch2 = +0.52 C for 1850-1900).

Share of bog pixels facing an **unprecedented** climate (outside the bio10 x bio12 hull of what any
bog lived through since t0), TraCE / EPOCH-2:

| | t0 = 2 ka | 4 ka | 6 ka |
|---|---|---|---|
| SSP3-7.0 | 99.7 / 99.5% | 99.3 / 99.2% | 98.9 / 99.0% |
| SSP1-2.6 | 56 / 82% | 40 / 58% | 9 / 11% |

- The two records agree on the pooled hull at 6 ka and on SSP3-7.0 everywhere; they differ for
  SSP1-2.6 at 2-4 ka (EPOCH-2's 1000-yr slices step later than TraCE's centennial ones).
- EPOCH-2's warmer HTM does not rescue SSP1-2.6 at the pixel level: its 6 ka peak is +1.37 C against
  ~1850 but only ~+0.8 C against 1981-2010, below the SSP1-2.6 future (+1.2 C). Under both records
  >= 95% of pixels are warmer in summer than any slice of their own past, by a median +0.4 C
  (EPOCH-2) to +0.7 C (TraCE) at t0 = 6 ka; SSP3-7.0 by +2.7 to +3.0 C.
- Applying D_trace lowered TraCE's pooled-hull coverage (SSP1-2.6 at 4 ka: 88% -> 60% inside), which is
  why the 2026-10-06 numbers below are superseded.

## 0a. First results (2026-10-06, TraCE only, no baseline offset; superseded by 0)

Outputs: `output/pl2/climatehistory_*.csv|png`. 1,116 NO bog cells -> 614 TraCE pixels in 6 approximate
CCSM3 cells (510 in one), so effective n is ~1-6.

| | t0 = 0 | 6 ka | 10 ka |
|---|---|---|---|
| SSP3-7.0, pooled hull | 0.2% | 3.9% | 3.9% |
| SSP3-7.0, per-pixel exceedance (median) | +3.5 C | +2.6 C | +2.6 C |
| SSP1-2.6, pooled hull | 37% | 99% (t\* = 6 ka) | 99% |
| SSP1-2.6, per-pixel coverage / median exceedance | 0% / +1.2 C | 11% / +0.3 C | 11% / +0.3 C |

Reference, today's pooled NO+EU bog climate (hull): 21% (SSP3-7.0), 99% (SSP1-2.6).

- TraCE's simulated HTM is weak here: max 8-6 ka bio10 only +0.7 C over the last slice (mean +0.1 C).
  The SSP1-2.6 verdict (+0.3 C exceedance) sits inside any plausible proxy-vs-TraCE bias; the SSP3-7.0
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
   Reported as **exceedance** in degrees C: future bio10 minus the historical maximum at that
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
- The transparent **exceedance + hull + density-contour** set (3.5) carries the verdict;
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

### 3.7a Literature pass, first stab (2026-10-07)

Sources from `ms/ZoteroAnnotations` unless marked *web* (found 2026-10-07; abstract only, not yet in Zotero
or the `.bib`, so verify before citing). Ages are cal BP as published; "14C BP" where the source is uncalibrated.

**Fen-to-bog (ombrotrophication) ages**

| Site / compilation | Region | Fen-to-bog age | Kind of evidence | Source | JV check |
|---|---|---|---|---|---|
| Viss mosse | S Sweden, 55.9 N | ~7.7 ka | onset of Sphagnum peat over fen peat, one core | `edvardssonHolocenePeatlandDevelopment2012` (also thesis 2013) | confirmed |
| Viss/Hällaryd/Åbua bog pines | S Sweden | pines in Sphagnum peat 7.2-5.7 ka and 4.6-3.1 ka; Åbuamossen expands laterally 5-3 ka | lower bound on transition (trees grew on a bog that already existed) | `edvardssonHoloceneClimateChange2013` | rejected; Åbua RW analysis are much younger (p.16) |
| Bergslagen, 3 bogs | central Sweden, raised-bog N limit | Sphagnum peat 5-4 14C ka (~5.7-4.5 cal ka); fens 8.6-3.5 14C ka | basal Sphagnum peat dates | `almquist-jacobsonIntegratedModelRaisedBog1995` | confirmed |
| Søndre Akershus | SE Norway | floating mats from 8.7 ka, tarn overgrown ~6.7 ka, ombrotrophic hummocks after | one basin, pollen/macrofossil | `sorensenHolocenVegetasjonshistorieOg2015` | p185 dates fen-bog transistion to "phase" 7600-2500 yr bp |
| Sellevollmyra, Andøya | N Norway, northernmost concentric raised bog | 440 cm core spans ~7 ka | *web*; whether all 7 ka is ombrotrophic is not yet confirmed | Vorren et al. 2007, Boreas | revised to 6410-5930 cal. yr BP |
| Lofoten, two bogs | N Norway | humification record over 7.4 ka | *web*; ombrotrophic peat implied, not confirmed | Vorren et al. 2012, Boreas |  "ombrotrophic since c. 3700 cal. a BP...and (ii) a mire that transformed into an ombrotrophic bog after c. 2000 cal. a BP"
| Finnish fen-bog ecotone | Finland | mostly late Holocene; final ombrotrophication 3-2 ka | *web*; macrofossils + 14C | Väliranta et al. 2017, The Holocene |
| Northern peatlands, 7,571 cores (329 transition dates) | circumboreal | peak fen-to-bog 5-3 ka | *web*; AGU abstract only | MacDonald, Holmquist, Kremenetski & Loisel 2015 |

Norwegian context: Lyngstad (2023, p. 24) describes raised bog as the end point of "many thousand years" of
undisturbed peat accumulation, starting "soon after the last glaciation"; mean mire depth in lowland S Norway
is ~2 m (p. 11). That is about peat initiation, not ombrotrophication.

**Revision after JV's check (2026-10-07).** The Edvardsson bog-pine row is dropped (Åbuamossen's
chronologies are much younger than the band they were used for). The Norwegian dates now carry the band:

| Norwegian site | Ombrotrophic since | Note |
|---|---|---|
| Sellevollmyra, Andøya (`vorrenHighresolutionStratigraphyNorthernmost2007`) | 6.41-5.93 cal ka | central core of the northernmost concentric raised bog; ombrogenic peat 415 cm there, 250 cm towards the edge |
| Rystad 1, Lofoten (`vorrenClimateChangesLast2012`) | ~3.7 cal ka | deep watershed bog |
| Rystad 2, Lofoten (same) | after ~2.0 cal ka | marginal mire turned bog |
| Rustadmåsan, S Akershus (`sorensenHolocenVegetasjonshistorieOg2015`) | gradual, within 7.0-2.5 ka (pollen zone 5) | basin overgrown 7.6-7.0 ka, initially minerotrophic; the light Sphagnum peat is from the last ~2 ka |

The two Swedish anchors stay as the regional bounds: Viss mosse 7.7 ka (southernmost, warmest; a pine
forest on the bog by ~7.2 ka, `edvardsson2010development`), and Bergslagen 5.7-4.5 cal ka (the raised-bog
northern limit, most like Norway's climate).

**How to read it.** Two things work against a single early date:
- *Basal is not ombrotrophic.* Basins in Norway were overgrown early (Rustadmåsan 7.6-7.0 ka) and peat
  initiation is often early Holocene, but ombrotrophy came later and gradually
  (`almquist-jacobsonIntegratedModelRaisedBog1995`: "basal" dates refer to fen peat).
- *Lateral expansion.* A dome is old at its initiation locus and young at its margins: Sellevollmyra's
  ombrogenic peat thins from 415 to 250 cm towards the edge, and Rustadmåsan was overgrown from the margin
  inwards. Lyngstad cells cover whole polygons, so most cells are younger than the oldest date for their
  bog. Granlund et al. 2022 show fen-to-bog fronts still advancing.

**Revised band: t0 = 2-6 ka (core), 1-8 ka (outer), centre ~4 ka.** Upper core edge: the only dated
Norwegian central-dome onset (Sellevollmyra, 6.4-5.9 ka). Lower core edge: the late Norwegian transitions
(Rystad 2 after 2.0 ka; Rustadmåsan's surface Sphagnum peat) and the late-Holocene Finnish and circumboreal
peaks (*web*, still to verify). Outer range: 8 ka admits the earliest onsets anywhere in the region (Viss
mosse 7.7, Rustadmåsan's possible early start); 1 ka admits the youngest margins. The band is for the
**typical bog cell**, not the oldest point in a bog. The previous 3-7 ka was pulled up by the rejected
bog-pine row and the Akershus misreading.

**Consequence for the verdict** (numbers from the TraCE-only run without baseline offset; superseded by
section 0, kept for the record). SSP3-7.0 is novel across the whole band (pooled hull <= 4%, per-pixel 0%).
SSP1-2.6: pooled-hull coverage rises from 76% at 2 ka to 88% at 4 ka and reaches 95% only at ~6 ka, the
top of the core band; per-pixel coverage stays at 5-11%. So "inside" under SSP1-2.6 needs the oldest
dated Norwegian onset for every cell, and even then only in the pooled sense.

**Proxy HTM against TraCE** (plan 3.8). Northernmost Fennoscandia: HTM maximum 1.5 +/- 0.5 C above present,
stable 8-6 ka (`sjogrenOverviewHoloceneClimate2021`). N European pollen stack: HTM 8.0-4.8 ka (Seppä et al.
2009, Clim. Past; *web*, amplitude not yet read). TraCE at our pixels: +0.7 C max, +0.1 C mean (8-6 ka minus
last slice). So TraCE likely understates HTM summer warmth by ~1 C. Read the exceedance against that: SSP1-2.6
(+0.3 C) is inside the gap, SSP3-7.0 (+2.6 C) is beyond the proxy maximum even at its upper end (2.0 C).
Karger et al. 2023 also note temperature-coupling errors between 10 and 6 ka
(`kargerCHELSATraCE21kHighresolution1thinspkm2023`).

**Still open:** Väliranta 2017 and MacDonald et al. 2015 (*web*, not in Zotero) back the lower edge
only alongside the Norwegian dates and can be dropped if unverified; a Vestlandet/Trøndelag onset date
would test the 6 ka upper edge in the bogs' core region; the Seppä 2009 HTM amplitude.

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

### 3.9 A second, independent temperature record: EPOCH-2 in parallel with TraCE (revised 2026-10-07)

**Problem.** At the 614 bog pixels TraCE gives a mid-Holocene (8-6 ka) bio10 only +0.1 C above its last
slice (warmest single slice +0.7 C). Pollen stacks put July temperature at 8-6 ka about **+1.5 C** above
present (`seppaLastNinethousandYears2009`; its sites include several in the bogs' region), and Karger et
al. (2023, Fig. 5) show TraCE21k and CHELSA-TraCE21k underestimating the Greenland ice-core HTM. That
error is **not systematic across time** (Karger Fig. 5), so TraCE is not corrected or blended.

**Decision (JV, 2026-10-07).** Use @mauriClimateEuropeHolocene2015 (EPOCH-2) **as it is, in parallel
with TraCE**: two independent records of the bogs' past, each run through the same containment analysis,
both shown in the verdict figure. No hybrid, no bias correction of one by the other. EPOCH-2 is at least as
trustworthy as TraCE for a summer-temperature question.

**Data.** NOAA, https://www.ncei.noaa.gov/pub/data/paleo/pollen/epoch2/ — `epoch-2.nc`,
`epoch-2_uncertainties.nc`, downloaded to `data/EPOCH2/` on 2026-10-07. 1 degree grid (63 x 44 cells,
28.5-71.5 N, 11.5 W-50.5 E); 13 slices at 12000, 11000, ..., 1000, 100 cal BP; anomalies in C (temperature)
and mm/day (precipitation). Used: `tanom_jja`; for the pooled hull, `panom_ann`; uncertainties
`tanom_jja_se_*`.

**First check (2026-10-07).** Bilinear `tanom_jja` at the 614 pixels, none NA. Median by slice: 8 ka -0.34,
7 ka +1.24, 6 ka +1.35, 5 ka +1.09, 4 ka +0.81, 3 ka +0.31, 2 ka +0.48, 1 ka +0.42, 100 BP 0.00. The 100 BP
slice is exactly 0 everywhere over Norway: **the anomalies are relative to ~100 cal BP (~1850 CE)**. The
6-7 ka values agree with Seppä's +1.5 C; the 8 ka dip is the 8.2 ka event.

**Construction (EPOCH-2 track).** For pixel i and slice k: bio10_E(i, k) = CHELSA_1981-2010(i) -
D_E(i) + tanom_jja(i, k). Only the 13 slices themselves enter the window (ages <= t0, capped by
deglaciation as for TraCE); nothing is interpolated between them. For the pooled hull, bio12_E(i, k) =
CHELSA_1981-2010 bio12(i) + 365 x panom_ann(i, k), floored at a small positive value; the preindustrial
precipitation offset is ignored (moisture never binds, see 0). JJA is taken as bio10's season; check in
CHELSA which quarter is warmest at the pixels and note any mismatch. Uncertainty band: repeat with
tanom_jja +/- 1 SE.

**Decision 1: putting each record on the 1981-2010 baseline (open).** The future is CHELSA 2071-2100,
i.e. an absolute climate. Each past record arrives as an anomaly against its own "present": TraCE against
its last slice (1900-1990 mean), EPOCH-2 against ~1850. To compare with the future, each anomaly must be
added to the CHELSA 1981-2010 climate *minus the warming between its own present and 1981-2010*:
- D_E = warming from ~1850 to 1981-2010 (EPOCH-2);
- D_T = warming from 1900-1990 to 1981-2010 (TraCE). The current analysis silently sets D_T = 0, which
  makes TraCE's past slightly too warm.

D is not a detail: it moves every past slice by the same amount. With the median numbers above, SSP1-2.6
is +1.2 C above 1981-2010 and EPOCH-2's 6 ka peak is +1.35 C above ~1850. If D_E = 0, the future is inside
the HTM; if D_E = 1 C, the HTM is only +0.35 C above 1981-2010 and SSP1-2.6 is ~0.85 C beyond it. Norway's
observed summer warming since the late 19th century is of this order, so the SSP1-2.6 answer depends on D_E
about as much as on the record.

Recommended: take both offsets from **one observational dataset** that reaches back to 1850, so the two
records are treated identically: Berkeley Earth gridded land temperature (1 degree, monthly, from 1850),
JJA mean at each pixel. D_E = JJA(1981-2010) - JJA(1850-1900) (1850-1900 is the IPCC AR6 preindustrial
reference); D_T = JJA(1981-2010) - JJA(1900-1990). Report both values. Alternatives: a published
Norwegian value (simpler to cite, but the station series start ~1900, so ~1850 must be extrapolated); or
GFDL-ESM4's own historical warming (consistent with the future scenario, but a model, not observations).

**Presentation without making everything conditional.**
- Panel A: both records' bio10 trajectories (TraCE centennial line, EPOCH-2 millennial points with SE),
  future ranges at the right.
- Panel B: the coverage points for both records (two point styles, or two small facets).
- Text: state what holds under both records once (SSP3-7.0 is outside under both, to be confirmed), then
  where they disagree once (expected: SSP1-2.6 own-past coverage), with the reason (TraCE's weak HTM).
  No sentence carries more than one "if".

**Pipeline.**
- New `R/pl0_prepareEPOCH2.R` after `pl0_prepareCHELSATraceSeries.R` in RUNALL: downloads the two `.nc`
  files if absent, records the URL and slices with record_settings().
- Baseline offsets: a small `R/pl0_prepareBaselineOffsets.R` (Berkeley Earth download, JJA means for the
  two periods, D at the pixels), if decision 1 goes that way.
- `R/pl2_exploreClimateHistory.R`: a third `history` version `epoch2` on the 13-slice age grid (the window
  and coverage functions take the age vector as an argument instead of the TraCE one). TraCE delta gains
  the D_T offset. Curves, t*, exceedance and bootstrap extend to the new version.
- `ms/make_figures.R`: panels A and B as above; `tbl_climatehistory` gains EPOCH-2 rows.

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
3. **Exceedance map** (definition 1), Norway, one value per pixel at a chosen t0.
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
2. Exceedance and coverage curve for the delta version, bio10 only (fastest read of the answer).
3. Add bio12, the hull and density definitions, pixel-level bootstrap.
4. ecospat COUE; hypervolume sensitivity.
5. Literature pass for the t0 band; core-cell sensitivity.
6. Figures into `ms/make_figures.R`; text for Methods, Results, and the Discussion caveats
   in 3.8.

## 8. Folding the analysis into the manuscript (2026-10-07)

### 8.1 The argument: three lines of evidence, each covering another's blind spot

| Leg | Question it answers | Reference for "has this climate been seen before" | Its blind spot | Covered by |
|---|---|---|---|---|
| 1. Observed occupancy (space, no model) | How common are raised bogs today in the climates the future brings? | all training cells, NO + EU, today | occurrence mixes formation and persistence; the warm analogues are few and clustered; no site-level answer | 3 (persistence at fixed sites); 2 (site level) |
| 2. Three-class forest + transfer CV (space, model) | Which mapped bogs keep a raised-bog climate, and what do the others become? | the same training data, generalised | extrapolation; learns occurrence, not persistence | 1 (same evidence, no learner); 3 |
| 3. Climate history (time, no model of bogs) | Have Norway's bogs lived through this climate since they became bogs? | the bogs' own past, at their own sites | t0 unknown; palaeo records coarse and disagree on HTM warmth; centennial vs 30-yr means | 1-2 (no t0, no palaeo record) |

Legs 1-2 share a weakness that only leg 3 addresses: a raised bog is a landform that can persist in a
climate where it could not form (now Discussion point iii, which concedes "loss of suitability, not a
disappearance date"). Leg 3 asks the persistence question directly, at the bogs' own sites. Its own
weaknesses (t0, the record) do not touch legs 1-2. So the three are complementary, not redundant, and the
paper can say that the conclusion stands on three legs that fail in different ways.

What the legs share, and must not be sold as independent: **the future climate** (CHELSA GFDL-ESM4, both
scenarios) and the CHELSA 1981-2010 baseline. The independence is in the *reference* (what counts as
precedent), not in the future. Say this once.

**What they say** (current numbers):

- **SSP3-7.0, all three agree on loss.** Model: 0 of 544 polygons still raised bog (96% non-peat).
  Occupancy: 75% of the future bog cells land in the > 18 C, wet cell of the grid, where observed
  prevalence is 10-fold lower than at 16-18 C (0.36 vs 3.63%). History:
  98.9-99.0% of bog pixels face a climate no bog has lived through in the last 8 ka (both records), and
  only 21% of the future even lies inside today's pooled NO + EU raised-bog climate (reference hull).
- **SSP1-2.6, all three put the bogs at the edge, not beyond it.** Model: 104 of 544 polygons (19%) stay
  raised bog, none lies beyond the warmest training bog, the score still falls at 98%. History: the
  unprecedented share depends on t0, 56-82% at 2 ka, 40-58% at 4 ka, 26-36% at 5 ka, 9-11% at 6 ka:
  bogs at least 5-6 ka old have, as a group, lived through comparable summers in the HTM. Occupancy:
  98.7% of the SSP1-2.6 future lies inside today's pooled raised-bog climate.
  This is not a contradiction between the legs but the persistence reading: under SSP1-2.6 most mapped
  bogs lose the climate in which raised bogs are common today (legs 1-2), yet a climate like it has been
  lived through by bogs during the HTM (leg 3). SSP1-2.6 is where "loss of suitability" and "loss of the
  bog" can come apart; SSP3-7.0 is where they cannot be told apart on any evidence we have.

**Writing rules.** State each scenario's verdict once, with the legs together, in Results and Discussion;
no sentence carries more than one "if" (the t0 dependence is carried by the figure and the band, not by
nested conditionals). No leg is presented as a "validation" of another. The two palaeo records are
reported as a range ("26-36%"), never as a primary and a check.

### 8.2 Where it goes (to-the-point: about +450 words main-text Methods, +300 Results)

| Place in `ms.qmd` | Change |
|---|---|
| YAML abstract skeleton | add the history leg as point (5); the CV point moves to (6) |
| Study system, question and design | one short paragraph after the five design decisions: three lines of evidence; roadmap sentence updated |
| Predictors | fill "[TODO cite]" for CHELSA-TraCE21k |
| New Methods subsection after "Model-free description of the projection problem": **"The climate history of today's raised bogs"** | 4 short paragraphs: question and unit (pixels, bio10 x bio12); two records and the baseline; window, t0 and the measure (pooled hull, "unprecedented"); the t0 band from the literature. Caveats in one sentence each, pointing to S2 Text |
| New Results subsection after the occupancy subsection: **"The future against the bogs' own past"** | 2 paragraphs + the figure (it becomes Fig 3; Quarto renumbers); sensitivities in one sentence pointing to S7 Table |
| "A milder scenario" | one sentence tying the 19% kept polygons to the history's SSP1-2.6 reading |
| Discussion TODO list | (i) rewritten as the three legs; new (viii) persistence vs suitability, sharpened by the history (SSP1-2.6 vs SSP3-7.0); new (ix) the palaeo caveats (TraCE's weak HTM, why two records, centennial smoothing biases toward "unprecedented", rate and CO2, t0 is an assumption) |
| Supplement | **S2 Text** (records, baseline offsets, window and deglaciation cap incl. the 18 pixels without an ice-free date, t0 literature table, measures and sensitivities, caveats); **S7 Table** (`tbl_climatehistory`) |
| `setup` chunk | read the `ch_*` study numbers and `tables/tbl_climatehistory.csv`; GUARDS for the qualitative claims (SSP3-7.0 >= 95% unprecedented across the outer band under both records; SSP1-2.6 unprecedented share falls below 50% inside the core band under both records) |

Deliberately left out of the main text: the per-pixel "own past" measure (JV: less relevant to persistence,
bio10 only; in S7 Table), the density-quantile definition, raw TraCE and bio18 (S7 Table), the GCM-cell
bootstrap (weak, 6 cells; one sentence in S2 Text), the core-cell sensitivity (not run).

### 8.3 Open for JV

1. Three legs or two? Recommended **three**, because occupancy and model are already presented as
   separate lines (Discussion point i, "the model-free occupancy table carries the conclusion and the
   model adds the map"). The history is then the third and the only one in time.
2. The persistence reading of SSP1-2.6 (8.1) is the main new claim for the Discussion; check that you
   want to make it.
3. Title and abstract: the working title is model-centred; a history leg might justify "three lines of
   evidence" in the title. Not changed.
4. Figure caption: drafted in the `ms.qmd` (based on the rough draft in chat).
