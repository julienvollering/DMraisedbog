# Implementation plan — EU data integration

Revised 2026-08-18 (implementation started same day). Supersedes the hierarchical
(two-scale nested) architecture.
Project state log stays in `notebook.md`; this file is the standing implementation
reference for the EU-data work.

---

## 1. Architecture decisions

**1.1 Drop `rf_global` as a local covariate; retire the hierarchical cascade.**
`pl0_modelGlobalScale.R` trains on the same 41 CHELSA+paleo stack the local model already
holds, at 5 km. The covariate therefore adds no predictor axis — only EU *label*
information, through a channel redundant with `bio10`/`gdd*` and carrying the global
model's background-design defects invisibly. VI rank 11/44 (0.0035 vs `bio10` 0.0129) is
consistent with this but is not the reason. The decisive reason: it puts the local model's
DI on two training sets at once — the same objection that killed pairwise RFQ
(notebook 2026-06-05, pt 3).

**1.2 Fix the pl3 anchor leakage.** `pl3_extractEUpredictors.R:3-5` claims the anchors are
leakage-free, but `pl0_modelGlobalScale.R:10` trains the global model on the thinned EU
presences and `pl3_extractEUpredictors.R:64` reads their `rf_global` from that model's
*in-sample* prediction surface. Deleting the `rf_global` block (lines 62–66) fixes it;
`feat` becomes 43.

**1.3 EU-in-training vs EU-in-evaluation is a dial, and the sweep is the experiment.**
The cut must be in **feature space** (withhold EU beyond a warm-flank cut), never a random
row fraction — a random holdout sits at low DI relative to what was kept and measures
nothing.

| Cut | Train | Test | Test DI |
|---|---|---|---|
| Norway's warm edge | Norway only | all EU | highest |
| mid | Norway + cool EU | warm EU | medium |
| warm | Norway + most EU | warmest EU | lowest |
| past all EU | Norway + all EU | — | (production) |

**DI is per-row (per-pixel) throughout; the cuts shift its *distribution*, they are not
themselves the unit of measurement.** The "Test DI" column names where each cut's DI
distribution sits, not a scalar the method returns. Each cut therefore yields a whole
(novelty, skill) *curve segment*, not a single point: bin the test rows by DI within a cut
and read skill per bin. The frozen ruler (§1.4) makes those bins commensurable across cuts,
so the segments can be overlaid — and if they superimpose, novelty is doing the explaining;
if they do not, *which* data was withheld matters beyond its novelty, and that is a result.
Do not naively pool all cuts' test rows into one binned curve: the warmest EU rows sit in
every cut's test set, so pooling double-counts them and correlates the errors. Carry
per-bin counts (`output/pl0/eu_absence_allocation.csv`) beside the curve — the warmest bins
are thin on the other-peat side (§6).

What the cuts buy is *range*: a span of novelty no single fit can produce, obtained by
design rather than by binning whatever novelty the folds happened to produce.
Operationalises "DI as forecast lead time".

- **Evaluation** = the sweep.
- **Production** = fit on everything; future Norway's DI read against the pooled training
  set, skill read off the curve.
- Assumption to state: the curve extrapolates to the production fit's own novelty level.
- Limits to state: the test set shrinks as the cut moves warmer (noisiest at the low-DI
  end); the curve stops where the warmest EU data stops — beyond that only an extrapolated
  trend exists, and saying so beats quoting a single accuracy.

**1.4 Freeze the DI ruler across the sweep.** DI is VI-weighted, so per-fit VI would change
the ruler between measurements. Take VI **and** the normalisation constant once from the
**production** fit (§5 step 8, Norway + all EU) and reuse for every cut — there is no
separate pilot fit, and the ordering works because the production fit precedes the sweep
(§5 step 11). The leak is a low-dimensional axis scaling, not
per-point information — state it in one sentence. Sensitivity: recompute the curve with
Norway-only-fit weights (no refits needed, distances only). Fallback if unstable:
unweighted / PCA-whitened distance. Report the systematic VI shift along the sweep
(warm-flank variables should gain importance) as a result; do not let it move the ruler.

**1.5 Encode differential data trust in the inclusion rule, not the loss.** No `case.wt`
trust weights. (i) `case.wt` is a *bootstrap sampling* weight — it shifts effective
prevalence, not label noise; a rarely-drawn noisy label still votes at full strength.
(ii) It conflates with class balancing if they share one vector. (iii) A scalar is global
over feature space while trust is not, and global EU downweighting attenuates exactly the
warm-flank contrast EU was recruited for. Instead: stated evidence bars for inclusion
(§3), thinning for the representativeness half, and dataset-blocked CV reporting
Norway-only vs pooled skill side by side. Any weight that survives must be **swept, not
tuned** — insensitive ⇒ drop it; sensitive ⇒ that is the finding.

**1.6 Integrated SDM — deferred, not rejected.** The right tool for genuinely disjoint
PO/PA data, and it would dissolve the pseudo-absence question inside a joint likelihood.
Three standing objections: (a) **identification** — separating true intensity from
observation bias works by borrowing strength across data types *over shared support*, and
Norway PA / EU PO support is disjoint, so the bias field is extrapolated precisely into
warm continental Europe, failing in the direction that makes bogs look *safer*; (b) **no
temporal transfer** — an SPDE spatial random field absorbs residual paleo/hydrological
structure and has no principled projection under climate change, which is disqualifying for
a projection RQ; (c) **wrong deliverable** — DI/AOA reliability mapping is defined against a
training feature set and does not port to a latent-field model. Proportionate future use:
a coarse-resolution ISDM as a **cross-check on warm-flank intensity only**, not as the
production architecture. Full reviewer-facing form in §7.3.

**1.7 The EU block is presence–absence by construction; "different data types" is not an
argument this project may use.** §3 builds a designed PA dataset on the EU side — presences
under stated evidence bars, absences from `domain \ M` split by EPM `map_cat`, inside an
explicitly bounded domain. The presence-only character of the raw Natura 2000 / EUNIS
material is upstream of the block, not a property of it. Consequences, both load-bearing:

- The PO/PA distinction may **not** be invoked to dismiss integrated SDMs (§1.6 does the
  work instead). It is self-defeating: any argument that PO and PA cannot be combined is
  equally an argument against pooling, which is the production architecture.
- What actually differs between the Norway and EU blocks is **absence provenance and label
  reliability**, not data type. That is handled by the inclusion rules (§3), thinning, the
  `dataset` column, and dataset-blocked CV (§1.5, §7.1) — measured, not assumed.

---

## 2. What each EU source is for

| Source | Role | Why |
|---|---|---|
| Natura 2000 7110/7120 | **presences** + exclusion mask | Polygon-based, hydrogeomorphological, grain-appropriate |
| Tanneberger mire region IV | **domain** | Covers Norway's future envelope (§4.1) |
| EPM2025 raster `map_cat` | **absence label split** | 0 / 1 / 2 categorical; no fractions needed |
| EUNIS `Wetland_Plot` Q11 | **exclusion mask only** | ~1 m² relevés cannot establish cell-level absence; and Q-codes are phytosociological, not hydrogeomorphological |
| EUNIS `Wetland.gpkg`, `Wetlands_Q_prob_100m.gpkg` | **never** | Modelled suitability/probability surfaces — training on them is circular |
| EPM2025 vector `peatl_type` | mask assist only | 60+ unharmonised raw values; disagrees with known raised bogs more often than it agrees (§4.4) |

**Rejected: EUNIS Q-codes as absence labels.** Poor fen and quaking mire are *constituents*
of raised bog complexes (lagg, hollows), so a Q22 relevé inside a bog complex becomes a
spurious absence — systematically, not randomly. And the response is active raised-bog
*function* (notebook 2026-06-05), a landform+hydrology concept; vegetation-type labels
would silently redefine it.

---

## 3. Rules for building the EU block

**Domain.** Tanneberger mire region IV **∩ EU-27** ∩ EPM2025 country coverage ∩ global
stack extent, **minus Norway**.

The EU-27 clip is not optional: Natura 2000 exists only in member states, so unclipped
Region IV would treat Russia, Belarus and Norway — 64.9% of its area — as absence domain
where raised bogs exist but are unmapped, producing systematic false absences concentrated
in the continental east. Norway is excluded separately because it is the Norway block.
The clip is nearly free in climate terms (§4.1). Do **not** extend the global stack east or
north for this reason: added extent would be non-Natura territory.

Verify the CRS transform explicitly — `mire_region_general.shp` has a `.prj` but no EPSG
code (`st_crs()$epsg` is `NA`). Use `ADM0_A3` rather than `ISO_A2` on the Natural Earth
countries layer: `ISO_A2` is `-99` for France and Norway, and France holds more raised-bog
polygons (3,207) than any other country.

**Response — 3 classes, matching the Norway block** (`non-peat` / `other-peat` / `bog`),
per notebook 2026-06-05. Norway gets it from `artype` + Lyngstad; EU from EPM + the
raised-bog mask. Same response definition on both blocks is what makes pooling coherent.

**Exclusion mask `M`** — maximally inclusive; deliberately *broader* than the presence set,
since a broader mask removes false absences at no cost to the presence definition:

1. Natura 2000 `HABITATCODE` 7110, 7120 **and 7130** (blanket bog) — **all**
   `DATAQUALITY` values, **all** `REPRESENTATIVITY` grades, all sizes.
2. EUNIS Q11 **and Q12** (blanket bog) plot locations, buffered.
3. EPM vector polygons with bog-like `peatl_type` (`bog`, `raised bog peat`,
   `intermediate bog peat`, `blanket bog & peatlands`, `quaking bog`, `string_bog`,
   `Mountain blanket bogs`, `Lowland Atlantic blanket bogs`, …) where available.
4. Buffer the union by ≥1 cell.

**Presences `P`** — the existing route, plus one filter:

5. Natura 2000 7110/7120 (**not** 7130 — blanket bog is mask-only), `DATAQUALITY == "G"`,
   within the domain.
6. Small polygons (<1 km²) taken directly; Q11 points inside large polygons as now.
7. **Exclude `REPRESENTATIVITY == "D"`** ("non-significant presence", 145 rows).
   **Retain** the 47 rows where the field is NA: excluding on missing metadata would induce
   bias correlated with national reporting completeness — the same error pattern avoided
   with `PRECISION` — and 47 rows is immaterial either way.

**Site attributes are not used.** `COVER_HA`, `PERCENTAGE_COVER`, `CONSERVATION` and
`RELSURFACE` are not trusted and play no part in any rule. The cover-fraction recovery arm
considered in §4.3 is **dropped**.

**Absences `A`** = domain \ `M`, split categorically by EPM `map_cat`:

| `map_cat` | meaning | label |
|---|---|---|
| 0 | no peatland | `non-peat` |
| 1 | peat dominated | `other-peat` |
| 2 | peat in soil mosaic | **excluded** — this class *is* the mixed-cell problem, already labelled by the data producers |

8. Majority-resample `map_cat` to the target grid. **No fractional cover layers**, here or
   in presences.
9. NA handling: inside a covered country NA means 0; outside coverage it means unknown.
   These must not be conflated when the raster is built. Uncovered cells are excluded, not
   defaulted to non-peat.

**Discarded presences are excluded, never relabelled.**

10. Large Natura polygons that fail the presence rules remain in mask `M`, so they are
    excluded from both absence classes. They become unlabelled — never `other-peat` or
    `non-peat`. This is structural under rule 1, not a filter that can silently fail.

**Sampling and bookkeeping.**

11. `non-peat` is effectively unlimited, so stratify the draw over climate space — never
    uniform-random (that failure is already logged for `pl0_modelGlobalScale.R`).
    Explicitly guarantee representation in the 1,738 warm-flank cells.
12. Retain a `dataset` column (EU / NO) throughout; report class counts per dataset; set
    balancing so the EU block does not silently reset the effective prior.

---

## 4. Evidence base

### 4.1 The domain covers Norway's future envelope, and the EU-27 clip is nearly free

Two load-bearing checks. Region IV is partly *defined* by mire occurrence, and the EU-27
clip removes 64.9% of its area — either could have cut off the warm flank. Neither does.

| Domain | cells @5 km | `bio10` min/med/q95/q99/max | warm-flank cells |
|---|---|---|---|
| Region IV, full (1,412,689 km²) | 23,365 | 6.1 / 15.8 / 17.2 / 19.9 / 21.0 | 1,743 |
| **Region IV ∩ EU-27 (495,223 km²)** | 21,218 | 6.1 / 15.9 / 17.3 / 20.1 / 21.0 | **1,738** |
| Norway current | — | 6.2 / 11.8 / 15.0 / 15.8 / 17.0 | — |
| Norway future | — | 9.5 / 14.9 / 17.1 / 18.8 / 20.7 | — |

Warm-flank cells are those above Norway's current maximum (17.0). The clip costs **5 of
1,743** of them; the clipped domain still contains **100%** of Norway's future cells with
**0%** above its q99. Area falls by two-thirds while cell count falls only 9%, because the
5 km stack extent already excluded most of Russia.

**Built and verified 2026-08-18** by `pl0_buildEUdomain.R`. Area 495,223 km² and the
per-country breakdown reproduce exactly; the bio10 envelope reproduces exactly under
`rasterize(touches = TRUE)`, which is now the stated cell rule (a 5 km cell partly inside
the domain still supplies EU information, and the centre rule drops border cells that
carry much of the warm flank). Cell count 21,226 against 21,218 above. The Norway
subtraction removed **0 km²**, confirming it is a no-op under the EU-27 clip. Cyprus and
Malta are the only EU-27 states EPM does not cover, and neither touches Region IV.

**The warm-flank threshold survives the reframing — corrected.** An intermediate reading
took Norway's envelope from the full 250 m raster (max 17.9) and concluded that warm-flank
supply collapsed from 1,975 cells to 702. That was wrong: the 250 m stacks cover terrain
Lyngstad never surveyed, so their maximum overstates what Norway's *training data* reach.
Measured over the actual Norway training population (`pl0_labelNorwayBlock.R`), current
bio10 runs **1.9 / 11.5 / 15.9 / 16.5 / 17.1** — a maximum of **17.1**, essentially the
17.0 this section already used. The warm flank is intact and rule 11 stands.

The threshold is therefore **published by `pl0_labelNorwayBlock.R`, not by the domain
script**, and `pl0_sampleEUabsences.R` consumes it from there. It is a property of the
training population, not of the raster extent, and every cut in the §1.3 sweep moves with
it. Norway's future cells remain 100% inside the domain's bio10 range, with 0.4% above its
q99 (the §4.1 "0%" was a masked-population number).

Region IV area by country (>50 km², `ADM0_A3`): Russia 813,511 · Finland 123,044 ·
Sweden 117,372 · Latvia 64,316 · Belarus 61,287 · Lithuania 61,187 · Norway 42,668 ·
Estonia 41,978 · Poland 28,104 · Austria 25,525 · France 22,480 · Germany 9,253 ·
Czechia 1,964. Not covered by Natura 2000: Russia, Belarus, Norway.

### 4.2 Presence loss from large polygons is severe but climate-neutral

Region IV holds 2,387 raised-bog polygons across 793 sites (7110: 1,643; 7120: 744);
1,319 polygons are `DATAQUALITY == "G"`.

| | polygons | area |
|---|---|---|
| Small (<1 km², used directly) | 613 | 160 km² |
| Large (≥1 km²) with a Q11 plot → retained | 51 | ~2,006 km² |
| Large without a Q11 plot → **discarded** | 655 | **15,729 km²** |

**99% of G-quality raised-bog area sits in large polygons; only 7% of those contain a Q11
plot.** Median site cover fraction is 0.21 — the bog is typically about a fifth of the
designated complex.

But the loss is approximately unbiased in climate. `bio10` (min/q25/med/q75/max):

| | |
|---|---|
| small, retained | 9.6 / 14.9 / **15.6** / 16.5 / 20.1 |
| large with Q11, retained | 13.4 / 14.7 / **16.2** / 16.4 / 17.8 |
| large, discarded | 8.7 / 15.1 / **15.9** / 16.6 / 18.6 |

Medians within 0.6 °C, ranges heavily overlapping, and the retained small polygons span the
widest envelope. **Loss of power, not bias** — report this table.

The mask is cheap: all 2,387 polygons cover 43,933 km² = **3.1% of Region IV** before
buffering, so generosity in `M` costs almost nothing in absence supply.

### 4.3 Cover-fraction recovery — considered and rejected

`PERCENTAGE_COVER` is entirely NA (0/2,387); `COVER_HA` is populated (2,386/2,387), so site
cover fraction was computable. Recovering discarded large-G sites at `cover_frac ≥ 0.75` and
representativity A/B would have added 79 sites / 581 km² — a 4× increase on the 160 km²
retained via small polygons.

**Rejected** (2026-08-17, user): the site attributes are not trusted, and the gain does not
justify the moving parts. The evidence agrees — the recovered envelope is *narrower* than
what is already held (`bio10` 14.4 / 15.5 / 16.9 against 9.6 / 15.6 / 20.1 for the retained
small polygons, nothing above 16.9). Recovery adds density in the core and nothing at the
flanks, so it was a power gain, not a coverage gain, and the warm flank — the reason the EU
block exists — is unaffected by dropping it.

### 4.4 Source characterisations

- **EUNIS `Wetland_Plot.gpkg`** — 92,468 plots, EPSG:3035, Level 3, 18 codes (Q11 5,357;
  Q12 1,170; Q2x 24,806; Q31 322; Q4x 14,514; Q5x 46,299); 40,218 unique locations, 5,688
  carrying >1 code. `PRECISION` is **100% NA for Q11** (0/5,357) while Q12 is 91%
  populated — so any precision filter would delete the presence class outright or induce
  class-dependent positional error. Moot now that Q11 is mask-only.
- **EPM2025** — a **compilation of national peatland maps**, not a climate-driven model,
  so using it as an absence source is not circular with respect to the CHELSA predictors
  (§7.4). Raster is 1 arcsec, EPSG:4326, per-country zips, values {1, 2} + NA.
  **Neither the NoData value nor the band type is harmonised across countries**
  (found 2026-08-18): DEU Byte/NoData 0, CZE+POL+AUT+FRA Byte/NoData 3, LVA+EST
  Byte/NoData 15, SWE+LTU Float32/NoData 99, FIN Int8/NoData -128. "No peatland" is
  encoded as whatever that country's NoData happens to be, so warping with the nodata
  mask on would resample the peat pixels only and lose class 0 outright. `map_cat` is
  therefore built by disabling the mask, taking the majority over raw codes, and
  remapping each country's own NoData code to 0 — a relabelling, not a recount — with an
  assertion on the resulting class set. Two further traps: the per-country rasters are
  **bounding boxes, not outlines**, so without a cutline a neighbour's background
  clobbers real data across the border; and reprojecting a lon/lat bbox into LAEA leaves
  uncovered corners that must stay NA rather than become a manufactured "no peatland". The
  vector `peatl_type` carries 60+ unharmonised values across ~55M polygons and at existing
  EU presences types them **fen 10,623 vs bog 4,468** — unusable as a bog/fen
  discriminator at point level.
- **The `peatl_type` bog-like mask (rule 3) reaches only 8 of the 10 domain countries**
  (measured 2026-08-18). Bog-like polygons: Finland 15,407,229 · Latvia 292,122 ·
  Germany 36,207 · Lithuania 22,405 · Estonia 13,396 · France 820 · Czechia 93 ·
  Poland 11 · **Sweden 0** · **Austria — no vector layer shipped at all**. Sweden's EPM
  vector simply does not use any bog-like type, so across 125,903 km² — the second-largest
  country in the domain — rule 3 contributes nothing and mask `M` rests on Natura 2000 and
  EUNIS alone. A Swedish raised bog mapped only as generic peat therefore survives into
  the absence pool labelled `other-peat`. This is a **country-level** gap, not random
  noise, so it cannot be waved away as attenuation: it belongs in the limitations and is
  a natural thing for the dataset-blocked reporting of §1.5 to pick up.
- **Global stack** — `predictors_global_5km_EUNorway_EPSG3035.tif`, extent 2.0–7.0 Mm ×
  1.0–5.5 Mm, 5 km, 41 layers. Those 41 plus `elevation` and `slope` are exactly the 43
  predictors the Norway 250 m stack carries once `artype_60` is removed, so the two blocks
  align by name with nothing left over.
- **EU terrain (`pl0_buildEUterrain.R`)** — DTM50 is Norway-only, so EU `elevation` and
  `slope` come from AWS Terrain Tiles via {elevatr} at z=9 (~87 m at these latitudes),
  aggregated to 250 m by mean and then differentiated to slope — the DTM50 order, not the
  reverse, which would read systematically steeper. Fetched in 111 cached 100 km chunks
  with a 3 km halo so slope is valid at the seams. Two source quirks worth knowing: the
  tiles carry **bathymetry** (the Baltic reads to −261 m, correct and simply not land) and
  occasional **decode artifacts** (6 cells above 5,000 m in the Ötztal Alps, where the true
  maximum is under 4,000). Values outside −500…5,000 m are set NA so a future re-draw that
  lands on one is dropped rather than trained on. At the drawn EU rows: 0 NA, elevation
  −11…2,985 m, slope median 0.97°.

---

## 5. Implementation sequence

Steps 1 and 2 are **swapped from the order they are numbered in**: the domain is what
bounds the EPM warp, and restricting the warp to the ten contributing countries is the
difference between minutes and hours over the 1-arcsec source. There is no circularity —
EPM *coverage* is a country-level fact read from the raster filenames, not from the
warped rasters.

1. **`pl0_buildEUdomain.R`** — **done 2026-08-18.** Region IV ∩ EU-27 ∩ EPM coverage ∩
   stack extent, minus Norway; explicit CRS check on the Mollweide shapefile. Writes
   `EU_domain` and `EPM_coverage` to the gpkg plus `output/pl0/eu_domain_5km.tif` and the
   §4.1 verification tables.
2. **`pl0_buildEPMmask.R`** — **done 2026-08-18.** Per-country warp under a cutline of
   country ∩ buffered domain → majority-resample `map_cat` to 250 m → merge → 5 km by
   `resample(method = "mode")` off the 250 m product (a second source pass would double
   the read cost for a diagnostic layer; the approximation labels nothing). Emits a
   separate **coverage** mask so in-country NA (= 0) and out-of-coverage NA (= unknown)
   stay distinguishable. Norway is excluded even inside the buffer: it is the other block.
   Result over the domain at 250 m — non-peat 90.2%, other-peat 9.3%, soil mosaic 0.5%.
3. **`pl0_buildEUterrain.R`** (new) — **done 2026-08-18.** EU `elevation` and `slope` at
   250 m (§4.4). Needed because `pl3_extractEUpredictors.R`'s per-point {elevatr} fetch is
   one network round-trip per row: fine for ~700 presences, hopeless for tens of thousands
   of absences.
4. **`pl0_collateEuropeanRaisedBog.R`** (revised) — **done 2026-08-18.** Presences per
   rules 5–7, exclusion mask `M` per rules 1–4. Result: **701 presences** (445 small
   Natura polygons + 256 EUNIS Q11 points inside 35 large polygons) and an `M` covering
   **13.2% of the domain** (69,289 km²), leaving 87% as absence supply. The §4.2 audit is
   recomputed in-script and still holds — medians within 0.5 °C across retained-small,
   retained-large and discarded, with the retained small polygons spanning the widest
   envelope: loss of power, not bias. Counts sit below the §4.2 table because that table
   predates both the representativity-D filter and the EU-27 clip.
   Buffers fixed at **one 250 m cell** (EUNIS plot buffer, and the final dilation of `M`).
5. **`pl0_labelNorwayBlock.R`** (new) — **done 2026-08-18.** Closes the §6 blocker.
   `artype_60` stops being a **mask** and becomes a **label**: presence -> `bog`,
   `>= 0.5` -> `otherpeat`, `< 0.1` -> `nonpeat`, `0.1-0.5` -> **excluded**, NA -> excluded.
   The excluded middle band mirrors EPM `map_cat` 2 on the EU side, so "other-peat" means
   the same thing on both blocks — without it Norway's mixed cells would land in
   `nonpeat` while Europe's equivalents are dropped. Banding of the 1,759,832-cell
   footprint: non-peat 1,322,932 · mixed 177,860 · other-peat 61,850 · `artype_60` NA
   197,190. The `>= 0.5` cut reproduces the current frame's 61,850 absences exactly, which
   is the check that the relabelling is the same operation seen differently. Result:
   **57,987 rows** (1,137 bog + 28,425 + 28,425). The survey footprint is *not* relaxed —
   absences are only meaningful where Lyngstad looked.
6. **`pl0_sampleEUabsences.R`** — **done 2026-08-18.** Same shared rule, applied to
   domain \ `M` split by `map_cat` per rules 9–12, from a candidate pool of 6,842,060
   cells. Result: **26,200 absences** (13,100 per class) + **524 presence cells**.
   Presences are counted as unique 250 m **cells**, not geometries: several EUNIS Q11
   relevés can fall in one cell of one Natura polygon and collapse to a single training
   row, so 701 geometries are 524 cells. Budgeting off the geometry count would make "per
   bog" mean different things on the two blocks and break the symmetry the shared rule
   exists to provide.

   **One sampling rule, stated once** (`draw_stratified_absences()` in `R/functions.R`),
   called by both block scripts with the same fixed 1-degree `bio10` bins, so a stratum
   index means the same thing on both sides and the dataset-blocked comparison of §1.5
   compares like with like. Budget = `ABSENCES_PER_BOG_PER_CLASS` (25) x that block's bog
   count, equal per absence class. Equal-per-class deliberately departs from
   area-proportional shares: bog is a climatic *intermediate* flanked by two opposite
   classes, so both flanks need real sample, and proportional draws left EU other-peat too
   thin to carry the contrast. Stated, symmetric, and **swept, not tuned** (§1.5). Pooled
   frame ~ **93,738 rows**, Norway:EU ~ 1.6:1 — neither block swamps the other.
7. **`pl2_createModelingFrame.R`** — **done 2026-08-18.** `rf_global` dropped (43
   predictors); the `artype_60` mask dropped and `artype_60` itself removed from the
   predictors so the Norway label cannot leak into the features; `dataset` retained; a
   guard fails loudly if the two blocks do not carry an identical predictor set. Pooled
   **current** frame:

   | dataset | non-peat | other-peat | bog | total |
   |---|---|---|---|---|
   | NO | 28,425 | 28,425 | 1,137 | 57,987 |
   | EU | 13,100 | 13,100 | 524 | 26,724 |

   Plus a **future** block (Norway projection domain): the block's own coordinates — so
   `pl2_exploreFeatureSpaceDistances.R`'s `(x, y)` join onto presence locations stays
   exact — union a 200,000-cell random sample of the projection domain, so the DI
   histogram still describes where Norway is going rather than only where the block drew.
   `scenario_current.tif` / `scenario_future.tif` rewritten at 43 layers, unmasked.
8. **pl2 re-run** — partitioning weights, tuning, CV. **Not a re-run:** `pl2_predict.R`,
   `pl2_evaluate.R` and `pl2_exploreFeatureSpaceDistances.R` all tested `response == 1` /
   `response == 0` and would have failed against a 3-class factor. That was deliberate —
   a compatibility column would have let them keep running on a response that no longer
   means what they assume. **Converted 2026-08-19:**

   - `pl2_weightFeaturesDataPartitioning.R` — RFQ dropped (two-class only); the ruler is
     now the production learner, one balanced 3-class `rfsrc()` fitted on the *whole*
     pooled frame with no absence subsampling, so the weights are literally the
     production fit's. `randomForest` BRF and a multinomial elastic net stay as
     cross-checks, and Spearman rank agreement against the production ruler is reported.
     The Norway-only ruler §1.4 asks for is emitted in the same file under its own
     `method` label, so the sensitivity curve needs no refits.
   - `pl2_exploreFeatureSpaceDistances.R` — prediction locations restricted to
     `dataset == "NO"` bog cells. **This was load-bearing, not cosmetic:** the future
     block covers the Norway projection domain only, so the unrestricted join returned an
     all-NA row for each of the 524 EU bog cells — 32% of the prediction set, silently.
   - `pl2_partitionDataByTopFeature.R` / `partition_by_presence_sorting()` — the anchor
     class is named (`presence_level = "bog"`) instead of implied by `== 1`; all other
     classes are assigned by range overlap exactly as absences were. Reports the full
     partition × class table and warns on an empty cell, since a partition missing a
     class silently turns that fold into a two-class problem.
   - `pl2_evaluate.R` — BRF-vs-RFQ comparison gone (one learner now); no threshold to
     optimise, since with three classes and a balanced forest the label is the argmax of
     the simplex. Reports per-class recall/precision/specificity/one-vs-rest AUC, macro
     G-mean, the full 3×3 confusion matrix, **and every metric split by the test row's
     `dataset`** — that split is §7.1 arm 1. Skill is binned over per-row DI as well as
     summarised per fold (§1.3).
   - `pl2_predict.R` — one production model; predictions are a 3-layer simplex per
     scenario. Also writes `output/pl2/di_ruler_production.rds` (weights + feature
     scaling + normalisation constant), which is §1.4's frozen ruler made concrete.
   - `pl2_interpret.R` — transitions are now class transitions at Lyngstad polygons.
     `bog -> otherpeat` is the reading the binary version could not produce at all, both
     classes having been the same `0`.
   - `calculate_weighted_di()` — rewritten. Exact nearest neighbour via kd-tree, and the
     normalisation constant from a subsample of training rows against all of them via
     chunked BLAS instead of an O(n²) R loop. Verified identical to the old
     implementation (max |ΔDI| = 1.1e-16 at equal settings). Gains `scaling` and
     `train_avg_dist` arguments so pl3 can freeze the ruler per §1.4.

   Not converted: `pl2_dissimilarityindex.R`. It is a superseded exploration of the
   cut-based and spacer-based partitioning schemes (replaced by presence-sorted
   partitioning), is in no pipeline in `RUNALL.R`, and still carries an interactive
   `?trainDI` call. It has the same binary assumptions and the same EU-presence join
   defect; it should be deleted or rewritten against the frozen ruler, not patched.
9. **Paired arms on Norway folds** — Norway-only fit and pooled fit scored on the *same*
   Norway test folds, so the with/without-EU comparison is made where the deliverable
   lives (§7.1, arm 1). Keep the pre-EU pl2 results alive as a reported arm; do not
   overwrite them with the pooled re-run.
10. **`pl3_extractEUpredictors.R`** — **done 2026-08-18.** `rf_global` block deleted and
   the stale "leakage-free anchors / not training data" header corrected; `feat` is 43,
   so the script's own guard now fails until step 5 lands, which is the intended signal.
11. **pl3 sweep** — cut positions along the warm gradient; one fit per cut; frozen DI ruler
   per §1.4.
12. **Projection agreement** — future-Norway projection under the Norway-only fit and under
    the production fit; correlation, area above the operating threshold, and the
    **disagreement map** read against DI (§7.1, arm 2). Plus the DI shift of future Norway
    against Norway-only vs pooled training (§7.1, arm 3).
13. **pl3 error profiles and reliability map.**

`pl0_modelGlobalScale.R` stops being a pipeline dependency. Keep it as a reported
side-analysis and as the source of EU presence coordinates.

---

## 6. Open items

- ~~**Norway-side 3-class response**~~ — **built 2026-08-18** (`pl0_labelNorwayBlock.R`,
  §5 step 4). The pipeline is unblocked; the next step is the pooled modelling frame.
- **Buffer distance for mask `M`** — "≥1 cell" is stated but not fixed; it differs between
  the 5 km and 250 m grids and should be set once, explicitly. `pl0_buildEPMmask.R`
  provisionally uses 5 km for the *domain* buffer it warps within; the mask buffer is
  still open.
- **`ABSENCES_PER_BOG_PER_CLASS` = 25 is a placeholder to be swept, not a finding.** It
  sets both blocks' size and therefore their relative weight in the pooled fit. §1.5's
  rule applies: insensitive => drop it from the discussion; sensitive => that is the
  result. Two earlier budget rules were tried and rejected — pinning to the *raw*
  Lyngstad footprint ratio (1,548 absences per presence) gives 1.08M EU rows against ~63k
  Norway rows, the exact prior-reset rule 12 forbids; pinning to area-proportional class
  shares leaves EU other-peat at 2,255 rows, too thin to carry the three-class contrast.
- **The warmest `bio10` bins are thin on the other-peat side** — 372 cells in the 18-19
  bin and 155 in 20-21, against 13,309 in 17-18. Non-peat is comfortable throughout. The
  very top of the §1.3 sweep therefore rests on a small other-peat sample; report the
  per-stratum counts (`output/pl0/eu_absence_allocation.csv`) alongside the curve rather
  than quoting skill there without them.
- ~~**`nonpeat` conflates climatic unsuitability with land-use conversion.**~~ **RESOLVED
  2026-08-19** by a land-use / water screen applied to BOTH blocks and ALL THREE classes
  (`pl0_buildLandUseScreen.R`, `landuse_screen_keep()` in `R/functions.R`). ESA WorldCover
  10 m, four fractional bands at 250 m; a cell is dropped if cropland+built-up **> 50%**
  or water+no-data **> 50%**. One threshold for both quantities, chosen for simplicity;
  the removal sweep over 25/33/50/67% is recorded in `notebook.md` and is **not** rerun as
  a sensitivity arm. CORINE was rejected on arithmetic — its 25 ha minimum mapping unit
  exceeds a 250 m cell's 6.25 ha, so it cannot express a within-cell fraction at all.
  Effect at the training cells: NO bog 98.2% kept, NO non-peat 92.2%, EU bog 95.2%, EU
  non-peat 82.2%. The blocks fail differently and the screen catches both — Norway drops
  mostly on **water** (1,456 cells vs 762 on human cover), Europe mostly on **human
  cover** (2,058 vs 276) — which is the concrete case for why screening Norway alone would
  have manufactured a warm-flank asymmetry rather than removing one.
- **Superseded diagnosis, kept because it explains the fix.** Found
  2026-08-19 when the first full projection put all 564 Lyngstad polygons into `nonpeat`.
  For 663 of 1,137 future Norwegian bog cells the nearest training analogue is a Norwegian
  non-peat cell, and those come from just **45 distinct cells** in the Oslofjord lowlands
  (median 9.26 E / 59.33 N) — flat (slope 0.73) and near sea level (5 m), so not non-peat
  for terrain reasons. `artype_60 < 0.1` cannot separate "climate cannot support peat"
  from "drained and ploughed", and the resulting error runs in the unsafe direction: it
  charges anthropogenic drainage to warming. Note this was *invisible* under the old
  `artype_60` mask, which deleted this population outright. Refined the same day: the anchors are
  **coastal**, averaging **48.5% water** against 5.1% for ordinary Norwegian non-peat (58%
  of them >25% water, vs 7%); human cover is enriched 3.3x but is the second effect, and
  forest is *not* implicated (0.239 vs 0.368 — less forested than ordinary non-peat). A
  built-up/cultivated screen alone catches only 12 of the 45. Norway's warmest maritime
  climate is on the coast, so coastal cells are the nearest analogue for warmed inland
  bogs while being non-peat largely because half the cell is sea. Note also that bog
  presences carry *more* human cover than random non-peat (0.092 vs 0.073), so any screen
  must run over all three classes with its effect on each reported, or it reshapes the
  presence set silently. Candidate fixes: a water-fraction screen (largest single effect),
  a built-up/cultivated screen per the existing artype_60 banding, or leave converted land
  **unlabelled** — the same
  treatment rule 10 gives discarded Natura polygons — rather than counting it as evidence
  of climatic unsuitability. Until resolved, the projection's headline number should not
  be quoted without it.
- **Does the thinner warm flank change the sweep design?** 702 cells is the whole
  high-novelty end of the §1.3 curve. Either accept a noisier top end and say so, or
  reconsider whether the cut axis should be `bio10` alone.

### Resolved 2026-08-17

- **Do not extend the global stack.** Added extent would be non-Natura territory, so the
  27% of Q11 records outside the current extent stay outside. Domain is clipped to EU-27
  instead (§3, §4.1).
- **EPM coverage vs Region IV countries** — subsumed by the EU-27 clip.
- **Blanket bog** → into mask `M` (Natura 7130, EUNIS Q12), never a presence.
- **Site attributes** (`COVER_HA`, `PERCENTAGE_COVER`, `CONSERVATION`, `RELSURFACE`) → not
  trusted, not used. The 8 sites with `COVER_HA` exceeding their own area are moot.
- **`REPRESENTATIVITY`** → exclude D from presences; retain the 47 NA rows.
- **`cover_frac` recovery arm** → dropped (§4.3).

- **The two blocks' climate is not at the same resolution.** Norway's predictors come
  from the 250 m regional stack, Europe's from the global 5 km stack, because no 250 m EU
  climate stack exists. CHELSA is 30 arcsec native so both are resamples of one source and
  neither is "the" native grid — but the EU rows are the coarser, and any feature living on
  fine topographic gradients is smoother on the EU side. Terrain is unaffected
  (`pl0_buildEUterrain.R` builds EU elevation and slope at 250 m). Decide whether to build
  a 250 m EU climate stack — the notebook's retired argument (b) says this is compute, not
  validity — or to state the asymmetry as a limitation. It is a candidate explanation for
  any systematic EU-vs-NO difference the §1.5 dataset-blocked comparison turns up, so it
  should be settled before that comparison is interpreted.
- **Sweden's missing bog-like `peatl_type`** (§4.4) — accept the gap and report it, or
  source a Swedish raised-bog layer separately. Not blocking, but settle it before the
  limitations are written.

### Resolved 2026-08-18

- **EPM circularity** → no issue. EPM2025 is a compilation of national maps, not a
  climate-driven model, so an EPM-derived absence label is not a re-import of a climate
  suitability surface (§4.4, §7.4).
- **"Different data types" as an argument** → retired (§1.7). The EU block is PA by
  construction; the argument would cut against pooling.

---

## 7. Anticipated objections and the parallel-arm design

Three objections are foreseeable and all three are answerable from work already planned,
provided the arms below are run and reported rather than assumed.

### 7.1 "The EU data are of doubtful quality — would you get the same answer without them?"

The §1.3 sweep does **not** answer this on its own. The sweep varies EU-in-training but
ties evaluation to **EU test points**; the objection is about **Norway conclusions**. Those
are two axes and only one is currently swept:

| | evaluated on Norway | evaluated on EU |
|---|---|---|
| **trained Norway-only** | pl2 as it exists (pre-EU arm) | sweep endpoint (“B”) |
| **trained pooled** | **arm 1 — must be run** | sweep interior / endpoint (“D”) |

Note also that sweep endpoint B is *not* the same object as “the analysis without EU data”:
B still uses EU to evaluate. The untouched-by-EU analysis is the existing pl2, which is why
§5 step 7 says keep it alive as a reported arm.

Three arms close the gap. None needs new architecture.

1. **Dataset-blocked CV on Norway folds.** Norway-only fit vs pooled fit, scored on the
   *same* Norway test folds. Already named in §1.5; now a numbered step (§5 step 7).
2. **Projection agreement, not just skill agreement.** The deliverable is a map, so compare
   maps: future-Norway projection under the Norway-only fit vs the production fit —
   correlation, area above the operating threshold, and the **disagreement map**. Prediction
   to state in advance: disagreement concentrates on the high-DI warm flank and is near-zero
   elsewhere. This figure is the reviewer answer and it is stronger than a table of AUCs.
3. **The DI shift.** Future Norway's DI against Norway-only training vs against pooled
   training. This is the affirmative case for the EU block, in the units the paper already
   uses.

**Frame the claim conditionally.** “We get the same answer without the EU data” is unsafe in
both directions: if the two projections agree everywhere, the EU block was unnecessary; if
they disagree, the naive reading is that EU data drove the result. The defensible claim is

> conclusions are unchanged wherever Norway-only training had support, and the fits diverge
> only in the warm-flank region where the Norway-only model was extrapolating — where the
> alternative to the EU data is not a different answer but an unsupported one.

Arm 2 *demonstrates* that claim rather than merely surviving it. It is the same statement as
the §1.3 skill-vs-novelty curve, read spatially.

### 7.2 "Why not a nested (hierarchical) SDM?"

Because it was built, measured, and retired — §1.1, and it can be shown rather than argued:

- The global model trains on the **same 41-layer stack** at coarser grain, so `rf_global`
  adds no predictor axis — only EU *label* information, through a channel redundant with
  `bio10`/`gdd*`.
- The nested design's purpose is importing broad-scale information the local domain lacks.
  Pooling does that **more directly**: as labelled rows whose novelty is measurable, instead
  of as a fitted surface that launders the global model's background-design defects into
  something un-auditable.
- Decisive: it puts the local model's DI on two training sets at once, which makes the
  paper's actual deliverable — the reliability map — uninterpretable.

VI rank 11/44 (0.0035 vs `bio10` 0.0129) is corroboration, not the argument.

### 7.3 "Why not an integrated SDM?"

Not because the data are of different types — see §1.7; that argument is unavailable and
would cut against pooling. The answer is in three parts:

1. **The part of an ISDM that would have mattered here was done explicitly.** The ISDM
   selling point is dissolving the pseudo-absence problem inside a joint likelihood; §3
   dissolves it with an independent peatland map and an auditable exclusion mask, where the
   assumptions are inspectable rather than buried in an identification argument.
2. **Identification fails where it matters** (§1.6a). Disjoint PA / PO support means the
   observation-bias field is extrapolated into warm continental Europe — and it fails in the
   direction that makes bogs look safer. This is the *spatial-domain* argument, and note it
   runs **against** ISDM, not for it.
3. **The random field does not transfer, and the deliverable does not port** (§1.6b, c).

### 7.4 "Your EU absences come from a mapped product — isn't that circular?"

No. EPM2025 is a **compilation of national peatland maps**, not a climate-driven model, so
an EPM-derived absence label does not re-import a climate suitability surface into the
training data. This is exactly why §2 admits EPM `map_cat` while permanently excluding
EUNIS `Wetland.gpkg` / `Wetlands_Q_prob_100m.gpkg`, which *are* modelled surfaces — the
distinction is principled, not convenience. State it in one sentence in the methods.
