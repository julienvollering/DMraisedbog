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
- **Decided 2026-10-07 (JV)**: no ecospat / hypervolume / COUE (step 4 dropped; 3.6 is superseded).
  The multidecadal-variability mismatch (3.8) is stated as a Methods caveat, not corrected for.
- **Literature t0 band**: 3.7a, revised after JV's check to 2-6 ka (outer 1-8 ka).
- **Figures/tables (step 6)**: `ms/make_figures.R`, section after Fig 2: `fig_climatehistory`,
  `figS_climatehistory_slices`, `figS_climatehistory_exceedance`, `tbl_climatehistory`, `ch_*` study numbers.
  Awaiting JV's approval before the manuscript text.
- **TraCE's weak HTM**: planned fix in 3.9 (proxy-anchored EPOCH-2 track), awaiting JV.
- **Not done yet**: core-cell sensitivity (step 5).

## 0. First results (2026-10-06, delta, bio10 x bio12)

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

**Consequence for the verdict.** SSP3-7.0 is novel across the whole band (pooled hull <= 4%, per-pixel 0%).
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

### 3.9 TraCE's mid-Holocene warmth: a proxy-anchored track from EPOCH-2 (option 1, planned 2026-10-07)

**Problem.** At the 614 bog pixels TraCE gives a mid-Holocene (8-6 ka) bio10 only +0.1 C above its last
slice (warmest single slice +0.7 C). Pollen stacks put July temperature at 8-6 ka about **+1.5 C** above
present, and annual temperature at +2.0-2.5 C (`seppaLastNinethousandYears2009`; its 36 sites include
Haugtjern, Svartvatnet, Dalane, Grostjørna and others in the bogs' region). Karger et al. (2023, Fig. 5)
show both TraCE21k and CHELSA-TraCE21k underestimating the Greenland ice-core HTM, and temperature-coupling
errors at 10-6 ka. The gap (~1-1.5 C) is larger than the SSP1-2.6 exceedance (+0.3 C at t0 >= 6 ka) and
smaller than the SSP3-7.0 one (+2.6 C). So the SSP1-2.6 verdict currently depends on TraCE's error, while
the SSP3-7.0 verdict does not.

**Approach.** Run a second temperature history next to TraCE and report the verdict for both: "model
(TraCE)" and "proxy-anchored (EPOCH-2)". Neither replaces the other: TraCE is internally consistent and
centennial; the pollen reconstruction has the HTM amplitude but is millennial and smoothed. If the two
bracket the answer, the manuscript says so.

**Data.** @mauriClimateEuropeHolocene2015. NOAA:
https://www.ncei.noaa.gov/pub/data/paleo/pollen/epoch2/ — `epoch-2.nc`, `epoch-2_uncertainties.nc` (also
GeoTIFFs under `tiff_anomalies/`, `tiff_uncertainties/`). 1 degree grid, 28.5-71.5 N, 11.5 W-50.5 E (covers
all bog pixels); 13 slices at 100, 1000, 2000, ..., 12000 cal BP; anomalies against "preindustrial".
Variables used: `tanom_jja` (summer temperature, C); kept for sensitivity: `panom_ann`, `anom_p_e`
(moisture balance). Uncertainties combine transfer-function and interpolation error.

**Construction (per pixel i, centennial slice s at age a_s).**
1. *EPOCH-2 at the pixels.* Bilinear from the 1 degree grid to pixel centres (nearest valid cell where the
   coastal cells are NA; count them). Interpolate linearly in age from the 1000-yr steps to every
   centennial age: a_E(i, s). Ages are cal BP (b1950); TraCE ages are b2k. Shift by 50 years once.
2. *TraCE anomaly and its millennial mean.* a_T(i, s) = TraCE bio10(i, s) - TraCE bio10(i, t20). Smooth
   with a centred 1000-yr running mean (11 slices) to m_T(i, s), the TraCE equivalent of a 1000-yr step.
3. *Hybrid anomaly.* a_H(i, s) = a_E(i, s) + [a_T(i, s) - m_T(i, s)]: EPOCH-2 sets the millennial level,
   TraCE adds the centennial wiggles around it. This keeps the warmest-century information that a pure
   1000-yr series would lose (gap 3 in section 2), without TraCE's level.
4. *Reference period.* EPOCH-2 anomalies are relative to preindustrial; CHELSA current is 1981-2010.
   bio10_H(i, s) = CHELSA_1981-2010(i) - D_PI(i) + a_H(i, s), where D_PI is the summer warming from
   preindustrial to 1981-2010 at the pixel. Source for D_PI (decision below): CRU TS JJA 1901-1930 vs
   1981-2010 at the pixels, which treats the early 20th century as preindustrial and so understates D_PI
   slightly (bias toward "inside"). A fixed regional value is the fallback. D_PI enters every past slice
   equally, so it shifts the history as a whole; report it next to the result.
5. *Precipitation.* Keep TraCE's ratio-delta bio12 for the primary hybrid track: moisture never binds
   (the future is wetter than every past slice), and pollen precipitation is the weaker reconstruction.
   Sensitivity: replace it with EPOCH-2 `panom_ann` (as a ratio on CHELSA) to check that a drier HTM does
   not change the pooled-hull result.
6. *Uncertainty.* Repeat with a_E +/- 1 SE from `epoch-2_uncertainties.nc` (`tanom_jja`) to give a band
   on the hybrid track. This is separate from the GCM-cell bootstrap, which stays weak (6 cells).

**Checks before reading results.**
- EPOCH-2 at the pixels, 8-6 ka, should come out near Seppä's +1.5 C (Tjul) for the same region; if it is
  far off, find out why before using it.
- bio10 is the warmest quarter, EPOCH-2 is JJA. Check in CHELSA 1981-2010 which quarter is warmest at the
  pixels (JJA vs JAS on the coast); if JAS is common, note the mismatch rather than correct it.
- The proxy HTM is a smoothed, multi-centennial signal; the 30-yr future is not. This is the same
  smoothing caveat as 3.8, now in the proxy track too (stated in the Methods, not corrected).
- Model-proxy disagreement on HTM warmth: the argument that proxies are seasonally biased concerns
  annual means. For a summer variable, the orbital forcing and summer-calibrated proxies agree on a warm
  HTM summer, which is why the proxy track is credible for bio10. Say this in the Discussion.

**Expected outcome (to test, not assume).** SSP3-7.0 stays novel: +3.5 C above 1981-2010 against an HTM
of ~+1.5 C above preindustrial, i.e. about 2 C beyond even after D_PI. SSP1-2.6 (+1.2 C above 1981-2010)
moves to roughly the proxy HTM level, so per-pixel coverage for t0 >= ~5-6 ka could rise from ~11% to a
large share. If so, the SSP1-2.6 sentence becomes conditional on which record is right and on bog age.

**Pipeline.**
- New `R/pl0_prepareEPOCH2.R` (after `pl0_prepareCHELSATraceSeries.R` in RUNALL): download the two `.nc`
  files to `data/EPOCH2/`, crop to the Norway window, write `tanom_jja`, `panom_ann` and their uncertainties
  as multi-layer GeoTIFFs; record_settings() with source URL and slices.
- `R/pl2_exploreClimateHistory.R`: add `history$hybrid` (and `hybrid_lo` / `hybrid_hi` for the SE band)
  built as above; D_PI computed or read in. `coverage()` already loops over versions, so the coverage
  curves, t*, exceedance per pixel and the bootstrap extend without new machinery. Write `a_E`, `m_T`
  and D_PI per pixel to `climatehistory_hybrid_inputs.csv` for checking.
- `ms/make_figures.R`: panel A adds the hybrid median trajectory (and the EPOCH-2 1000-yr points with
  their SE) over the TraCE one, so the gap is visible; panel B shows both tracks (design after JV's
  current figure review); `tbl_climatehistory` gains hybrid rows.

**Decisions for JV.** (a) D_PI source: CRU TS 1901-1930 vs 1981-2010 (recommended), or EPOCH-2's own 100
BP slice, or a fixed regional value. (b) Whether the hybrid (EPOCH-2 level + TraCE wiggles) or the plain
EPOCH-2 series (interpolated, no wiggles) is the proxy track; the hybrid is recommended, the plain series
is a check. (c) Whether the verdict figure shows both tracks or the proxy track goes to the supplement.

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
