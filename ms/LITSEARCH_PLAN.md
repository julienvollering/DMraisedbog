# Literature search plan for `ms/ms.qmd`

Drafted 2026-11-17. A plan for finding the literature the manuscript needs and for showing
JV *which part* of each paper matters, so that reading time goes on the passages that will
actually be cited. Nothing here has been run yet.

## 0. What the search has to serve

The search is **claim-driven**: every paper is looked for because a specific sentence in the
manuscript needs support, context or a counter-argument. The search starts from a **claim
register** taken from `ms.qmd`. Each claim has an ID, and every paper found is tagged with
the claim IDs it bears on.

My reading of the aims, to be confirmed by JV (see questions at the end):

- **RQ1 (ecological).** Which of Norway's mapped typical raised bogs keep a climate that
  supports raised bogs by 2071–2100 (SSP3-7.0; SSP1-2.6 as the milder case)?
- **RQ2 (evidence without a model).** Is the same answer visible in the observed prevalence
  of raised bogs in the climates the future brings, without fitting any model?
- **RQ3 (methodological).** How far can a model trained on pooled data, with European
  warm-flank evidence, be trusted when it is moved into warmer climates? This is measured
  as skill against the signed bio10 offset, plus a DI applicability mask.

Draft claim register (to be completed from `ms.qmd` before searching):

| ID | Section | Claim / need |
|----|---------|--------------|
| I1 | Intro | Raised bogs are ombrotrophic landforms with a climatically bounded distribution |
| I2 | Intro | The typical raised bog region of Europe; Norway as its north-western outpost |
| I3 | Intro | Earlier climate-envelope / bioclimatic models for bogs and mires, and what they found |
| I4 | Intro | Observed and projected climate impacts on bog hydrology, *Sphagnum* and peat accumulation |
| I5 | Intro | Management context: restoration prioritisation, Natura 2000 habitat 7110/7120, climate-informed conservation |
| I6 | Intro | Truncated niches: a single-country model cannot see the warm edge; pooling data across regions as a fix |
| M1 | Methods | Three-class / multiclass response in distribution modelling; habitat or ecosystem (not species) distribution models |
| M2 | Methods | Balanced random forests for imbalanced classes (`make.wt`, `make.size`) |
| M3 | Methods | Feature-space / environmental blocking in CV; evaluating transferability along a gradient |
| M4 | Methods | Dissimilarity index / area of applicability; extrapolation metrics (MESS, ExDet, etc.) |
| M5 | Methods | Climate-stratified background or absence sampling |
| M6 | Methods | Data-source citations: every `[TODO cite: …]` in `ms.qmd` (Natura 2000 dataset, AR50, EPM2025, WorldCover, DTM50, elevatr, CHELSA-TraCE21k, R, sf, terra, randomForestSRC) |
| D1 | Disc (i) | Model-free / occupancy-based evidence alongside correlative models |
| D2 | Disc (ii) | Do correlative models just default to absence outside their range? Evidence on extrapolation behaviour of RF |
| D3 | Disc (iii) | Formation versus persistence: self-maintaining bog hydrology, lags, extinction debt, resilience of existing peatlands |
| D4 | Disc (iv) | Geographic clustering and representativeness of warm-edge evidence; space-for-time substitution caveats |
| D5 | Disc (v) | Uncertainty from a single GCM; scenario versus model spread in impact studies |
| D6 | Disc (vi) | Terrain-source confounds when pooling data from different sources |
| D7 | Disc (vii) | How the absence or background design shapes apparent prevalence |
| D8 | Disc | Comparable projections for other mire types or regions (blanket bog, aapa mires, Fennoscandia, Baltic) and how our numbers compare |

## 1. Access and tooling

**Network.** `api.openalex.org` is currently **blocked** by this cloud environment's
network policy (proxy 403). Before running anything, the host must be added to the
environment's allowed hosts (I can propose that change and JV approves it on a card).
Optional extras, each an exact host:

- `api.unpaywall.org`: open-access PDF locations, as a cross-check on OpenAlex's own OA links
- `www.ebi.ac.uk`: Europe PMC full-text XML for OA papers
- `api.crossref.org`: metadata fallback and reference lists missing from OpenAlex

**Publisher PDF hosts are too many to allowlist one by one.** Full text will therefore come
from (a) OA copies that OpenAlex points to on hosts we allow, and (b) **JV's own Zotero
PDFs**. The `.bib` `file =` fields show they sit in `OneDrive/ZoteroLibrary` on JV's
computer. That folder can be attached to this session through Remote Control, or JV can
drop the relevant PDFs into a folder I can read. This also covers paywalled papers JV
already has access to.

**API key / polite pool.** OpenAlex now asks heavy users for a free API key. If it is
needed I'll request it through the secure credential form. I will not put JV's email in
request headers without asking.

**Code.** There is no R in this container. The search runs from a small Python script in
`ms/litsearch/` (stdlib + `requests`), writes plain CSVs, and so can be re-read or re-run
in R (`openalexR`) on JV's machine. The script sits outside `R/RUNALL.R`, like
`ms/make_figures.R`. Every query is logged with its timestamp and hit count, so the search
can be reported (PRISMA-style counts in S1 Text if the journal or reviewers want it).

## 2. Finding papers

Four complementary routes. Results are pooled and deduplicated on OpenAlex ID / DOI.

### 2.1 Seed set: what is already in the library
Resolve all 168 entries of `Repeat-raisedbogDM.bib` (131 have DOIs) to OpenAlex IDs. This
(a) prevents re-recommending what JV already has, (b) gives the seed papers for citation
chaining, and (c) lets me tag existing entries with claim IDs, which is often the fastest
win: several Discussion points can be supported from papers already in the `.bib`.

### 2.2 Keyword searches per claim block
Use OpenAlex `search` (title + abstract + fulltext where indexed) and
`title_and_abstract.search` with Boolean strings, filtered to `type:article|review`, and with
year limits where the field moves fast (methods ≥ 2010; climate projections ≥ 2005).
Example strings (to be refined after the first pass):

- I1/I3: `("raised bog" OR "ombrotrophic" OR "Sphagnum bog") AND ("climate envelope" OR "bioclimatic" OR "distribution model" OR "climatic limit")`
- I4/D3: `(bog OR peatland OR mire) AND (warming OR drought OR "water table") AND (persistence OR resilience OR "ecohydrological feedback" OR tipping)`
- I6/M3/D2: `("species distribution" OR "habitat distribution" OR "ecosystem distribution") AND (transferability OR extrapolation OR "novel climate") AND ("cross-validation" OR "environmental blocking" OR "truncated niche")`
- M4: `("area of applicability" OR "dissimilarity index" OR "MESS" OR "ExDet" OR "environmental novelty") AND model`
- M1/M2: `("random forest") AND (multiclass OR "class imbalance" OR "balanced") AND (ecolog* OR habitat OR land cover)`
- I2/I5/D8: `(mire OR peatland OR bog) AND (Norway OR Fennoscandia OR Sweden OR Finland OR Baltic) AND (climate change OR projection OR Natura 2000)`

### 2.3 Citation chaining (the main source of the best hits)
For about 15–25 anchor papers (e.g. Gallego-Sala et al. 2010, Parviainen & Luoto 2007,
Clark et al. 2010, Essl et al. 2012, Heikkinen et al. 2022/2023, Meyer & Pebesma 2021,
Ludwig et al. 2023, Chevalier et al. 2022, Linnenbrink et al. 2024, Tanneberger et al.
2021, Swindles et al. 2019, Belyea & Baird 2006):

- **forward** (`filter=cites:W…`): newer work that builds on or argues with the anchor, which
  is how 2023–2026 papers turn up;
- **backward** (`referenced_works`): the foundations the anchor rests on;
- **related works** (OpenAlex `related_works`): a cheap third net.

A paper reached from two or more anchors is ranked up (co-citation signal).

### 2.4 Targeted look-ups
For M6 (data and software citations) and named items, I look each one up directly
(DOI / exact title, or the data provider's recommended citation) rather than searching.
Norwegian and Nordic grey literature (NINA, NIBIO and Miljødirektoratet reports) is poorly
indexed in OpenAlex. For that I'll list what is missing and ask JV, rather than guess.

**Stopping rule.** Iterate per claim block until a new round of chaining adds < 5% new
relevant records, or the block has 3–5 strong references.

## 3. Screening and ranking

1. **Automatic de-duplication and pre-filter**: drop retracted, non-English/non-Scandinavian,
   and out-of-scope types; flag preprints that have a published version.
2. **Title/abstract screen** (by me, from the OpenAlex abstract). Each record gets:
   - `relevance` 0–3 against the claim register,
   - `claims` (one or more IDs),
   - `role`: *support* / *contrast* (disagrees or complicates) / *method* / *data* / *comparison*,
   - a one-line reason.
   Records whose abstract is missing in OpenAlex (common for some publishers) are marked
   `abstract_missing` and screened on title + venue, then checked in full text if they pass.
3. **Ranking within each claim**: relevance, then co-citation from anchors, then how
   directly the paper's study system matches (bogs > peatlands > other ecosystems), then
   citation count normalised by age (so recent papers are not penalised).
4. **Contrasting evidence is kept on purpose.** For each Discussion point at least one paper
   that complicates our reading is kept (e.g. evidence that bogs persist or expand under
   warming, or that RF extrapolation is better than feared), so the Discussion can address it
   rather than be surprised by a reviewer.

Expected volume (rough): 1,500–3,000 unique records → 150–250 pass title/abstract →
40–70 read in full text → perhaps 30–50 new references.

## 4. Highlighting the relevant parts of each paper

This is the part that saves JV's time. For each paper that reaches full-text review:

### 4.1 Extraction
Get full text from the OA copy or JV's PDF; extract text **with page numbers**
(`pdftotext -layout`, page by page; section headings detected where possible).

### 4.2 Passage finding
Read the paper against its tagged claim IDs and pick the 2–6 passages that bear on them.
Typical passages: a result with its number, a stated limitation, a method detail we
mirror or depart from, or an explicit disagreement.

### 4.3 Verbatim and verified
Every highlighted passage is an **exact quote**. A script checks that each quote occurs in
the extracted text of the cited page (after normalising whitespace and hyphenation), and
drops any that do not. No paraphrase is presented as a quote. If only the abstract was
available, the card says so.

### 4.4 Deliverables, in the repo (`ms/litsearch/`)

- `claims.csv`: the claim register.
- `search_log.csv`: every query/chain step, date and hit count.
- `candidates.csv`: all screened records (OpenAlex ID, DOI, title, year, venue, relevance,
  claims, role, reason, in_bib).
- `evidence/<bibkey>.md`: **one evidence card per paper**:
  - full reference + DOI, *why it matters for this manuscript* (2–3 sentences);
  - highlighted passages: quote, page, section, claim ID, role, and a one-line note on how
    it could be used ("cite after '…' in Methods §Partitions");
  - a caveat line (study system, scale or design differences that limit the comparison).
- `claim_matrix.md`: the **claim → evidence map**. For each claim ID: the supporting,
  contrasting and method references with their single best quote, plus the gaps where
  nothing good was found. This is the document to write the Introduction and Discussion
  from.
- `new_refs_dois.txt`: DOIs for Zotero's *Add Item by Identifier*, so new references enter
  JV's Zotero library (and thus the Better BibTeX export) instead of being hand-edited into
  the `.bib`.
- **Optional: highlighted PDFs.** Copies of the PDFs with the selected passages highlighted
  (PyMuPDF), colour-coded by role (e.g. yellow = support, red = contrast, blue = method,
  green = data). Each highlight carries a comment with the claim ID. Zotero shows
  annotations embedded in a PDF and can import them into its own annotation layer, so the
  highlights land where JV already reads. These go to a folder JV chooses (not the repo:
  copyright).

The cards and matrix are plain Markdown, so they diff in git and can be pulled into a
Quarto document if wanted.

## 5. Order of work

1. JV answers the questions below; network access is approved.
2. Build the claim register from `ms.qmd`; resolve the `.bib` to OpenAlex; tag existing
   entries with claims. **Checkpoint:** show JV the register and the coverage that the
   existing library already gives.
3. Keyword searches + chaining per block; screening → `candidates.csv`.
   **Checkpoint:** JV skims the top ~30 per block and vetoes or adds items.
4. Full-text pass on the shortlist → evidence cards, claim matrix, DOI list,
   optional highlighted PDFs.
5. Fill M6 (data/software citations) directly and close those `[TODO cite]` markers on a
   branch.
6. Commit each step on the task branch; record the search in `notebook.md` if it changes
   the manuscript's framing.

The blocks are independent, so steps 3–4 could run as parallel sub-agents (one per block)
if JV wants it faster. By default I'll run them in sequence in one session.

## 6. Known limitations

- OpenAlex abstracts are missing for part of the publisher catalogue; title-only screening is
  weaker, and those papers rely more on chaining.
- Grey literature (Norwegian reports, EU habitat reports) is under-indexed; JV's knowledge
  fills that gap.
- Highlights are only as complete as the full text I can reach; paywalled papers not in
  JV's Zotero get abstract-level cards only.
- Relevance judgements are mine. The checkpoints in §5 are where JV corrects them.

## Questions for JV before starting

1. **Aims/RQs.** Are RQ1–RQ3 above right? And is the paper primarily ecological (fate of
   Norway's bogs) with the transfer evaluation as support, or are the two equal
   contributions? This sets how much of the search goes to blocks M3/M4/D2 versus
   I1–I5/D3/D8.
2. **Target journal.** GCB or PLOS ONE (or other)? This sets the reference budget (roughly
   60–90 vs no hard limit) and how much methods literature belongs in the main text.
3. **Scope.** Include Scandinavian-language literature and grey reports via your input?
   Any time window?
4. **Full text.** May I attach your Zotero PDF folder (via Remote Control), or will you
   drop PDFs in a folder? Do you want the highlighted-PDF output, or are the evidence cards
   enough?
5. **Network.** OK to add `api.openalex.org` (and optionally `api.unpaywall.org`,
   `www.ebi.ac.uk`, `api.crossref.org`) to your Default cloud environment? That environment
   is used by all your sessions.
