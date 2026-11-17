# Literature search plan for `ms/ms.qmd`

Drafted 2026-11-17 in a cloud session; revised the same day after checking access on JV's
machine (where the search will run) and after JV answered the questions at the end. A plan
for finding the literature the manuscript needs and for showing JV *which part* of each
paper matters, so that reading time goes on the passages that will actually be cited.
Nothing here has been run yet.

## 0. What the search has to serve

The search is **claim-driven**: every paper is looked for because a specific sentence in the
manuscript needs support, context or a counter-argument. The search starts from a **claim
register** taken from `ms.qmd`. Each claim has an ID, and every paper found is tagged with
the claim IDs it bears on.

The aims, as confirmed by JV (answer 1):

- **RQ1 (ecological; the paper's question).** Which of Norway's mapped typical raised bogs
  keep a climate that supports raised bogs by 2071–2100 (SSP3-7.0; SSP1-2.6 as the milder
  case)?
- **RQ2 (a second, independent leg for RQ1).** Is the same answer visible in the observed
  prevalence of raised bogs in the climates the future brings, without fitting any model?
- **Validation step (not a research question).** How far the pooled model can be trusted
  when moved into warmer climates: skill against the signed bio10 offset, plus a DI
  applicability mask.

**Weighting of the search.** Most effort goes to the ecological blocks (I1–I5, D1, D3, D8)
and to the evidence-without-a-model framing (D1). The methods blocks (I6, M1–M4, D2) get
only what a GCB reader needs to accept the validation step: a few well-chosen references
each, not a review of the transferability literature.

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

The search runs locally on JV's Windows machine. Access was checked on 2026-11-17.

**Network.** All of these respond from this machine; no allowlisting is needed:

- `api.openalex.org`: the search, chaining and metadata backbone;
- `api.crossref.org`: metadata fallback and reference lists missing from OpenAlex;
- `www.ebi.ac.uk` (Europe PMC): full-text XML for OA papers.

Unpaywall is **not used**: it requires an email address in every request, and JV declined
that (answer 5). OpenAlex's own OA links cover the same ground.

**OpenAlex API key.** JV's personal key is in the user `.Renviron` as `openalexR.apikey`
(tested 2026-11-17: status 200, no shared-budget notice). That file lives in
`OneDrive - Høgskulen på Vestlandet/Documents/.Renviron`, which is where Positron's R looks;
`Rscript --vanilla` reads no `.Renviron` at all, and command-line R's home on this machine
is `C:\Users\julienv\Home`. So every search script starts with an explicit
`readRenviron("<that path>")`. The file stays where it is (it also sets `HOME`, `PATH` and
`R_LIBS` for Positron). The key is never printed, logged or committed.

**Full text.** From **OA copies only**: those OpenAlex or Europe PMC points to, downloaded
to a scratch folder outside the repo. JV's Zotero PDF folder (`ZoteroLibrary/`) is
**not read** (answer 4); the `file =` paths in the `.bib` are ignored.

**JV's existing reading notes.** `ms/ZoteroAnnotations/` holds one Markdown export per
`.bib` entry (168 files, named by Better BibTeX citekey) with JV's own highlights and page
numbers: 133 files contain highlights, 2,773 quotes in total. These are the evidence
source for the existing library (§2.1) and, together with OA full text, the only text read
beyond abstracts. The folder is untracked in git and exists only in the main checkout; it
is read from there, not copied.

**Code.** R, to match the repo: R 4.6.2 with `openalexR` 3.1.2 (installed in the user
library, which `Rscript --vanilla` finds), and `httr2` for Crossref / Europe PMC. Text is
extracted with `pdftotext -layout` (ucrt64, on PATH) or `pdftools`. Scripts live in
`ms/litsearch/`, are run with `Rscript --vanilla` from the repository root, and sit outside
`R/RUNALL.R`, like `ms/make_figures.R`. Every query is logged with its timestamp and hit
count, so the search can be reported (PRISMA-style counts in S1 Text if the journal or
reviewers want it). No Python is needed.

## 2. Finding papers

Four complementary routes. Results are pooled and deduplicated on OpenAlex ID / DOI.

### 2.1 Seed set: what is already in the library
Resolve all 168 entries of `Repeat-raisedbogDM.bib` (131 have DOIs) to OpenAlex IDs. This
(a) prevents re-recommending what JV already has, (b) gives the seed papers for citation
chaining, and (c) lets me tag existing entries with claim IDs, which is often the fastest
win: several Discussion points can be supported from papers already in the `.bib`.

Tagging starts from JV's own highlights in `ms/ZoteroAnnotations/<citekey>.md`: each
highlight is mapped to the claim IDs it bears on, so the evidence cards for papers JV has
already read are built from passages JV chose, with their page numbers. The 35 entries
without highlights (and any claim a paper's highlights do not cover) get an OA full-text
pass where an OA copy exists, and an abstract-level card otherwise.

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
Scandinavian-language and grey literature are out of scope (answer 3); data-provider
documentation needed for M6 is the only exception.

**Reference budget.** GCB, 60–90 references in total (answer 2). `ms.qmd` currently cites
about 23 keys, so the search should yield roughly 40–65 more, drawn from the existing
`.bib` first and new papers second. Rough split: ecological blocks (I1–I5, D1, D3–D5, D8)
about two-thirds; methods blocks (I6, M1–M5, D2, D6, D7) about one-third; M6 data/software
citations on top as needed.

**Stopping rule.** Iterate per claim block until a new round of chaining adds < 5% new
relevant records, or the block has its share of the budget in strong references
(typically 2–4 per methods claim, 3–6 per ecological claim).

## 3. Screening and ranking

1. **Automatic de-duplication and pre-filter**: drop retracted, non-English and
   out-of-scope types (incl. grey reports); flag preprints that have a published version.
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
40–70 get a card (OA full text where available, abstract otherwise) → perhaps 20–40 new
references, the rest of the budget coming from the existing `.bib`.

## 4. Highlighting the relevant parts of each paper

This is the part that saves JV's time. For each paper that reaches full-text review:

### 4.1 Extraction
For papers already in the `.bib`, start from JV's highlights in `ms/ZoteroAnnotations/`.
Otherwise get the OA copy (PDF or Europe PMC XML); extract text **with page numbers**
(`pdftotext -layout`, page by page; section headings detected where possible). Downloads
and extracted text go to a scratch folder outside the repo. No OA copy → abstract-level
card.

### 4.2 Passage finding
Read the paper against its tagged claim IDs and pick the 2–6 passages that bear on them.
Typical passages: a result with its number, a stated limitation, a method detail we
mirror or depart from, or an explicit disagreement.

### 4.3 Verbatim and verified
Every highlighted passage is an **exact quote**. A script checks that each quote occurs in
the extracted text of the cited page (after normalising whitespace and hyphenation), and
drops any that do not. Quotes taken from JV's own highlights are marked as such and keep
JV's page number (they are not re-verified, since the PDF is not read). No paraphrase is
presented as a quote. If only the abstract was available, the card says so.

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

No highlighted-PDF output: the evidence cards are the deliverable (answer 4).

The cards and matrix are plain Markdown, so they diff in git and can be pulled into a
Quarto document if wanted.

## 5. Order of work

1. ~~JV answers the questions below; the OpenAlex key is in place; `openalexR` is
   installed.~~ Done 2026-11-17.
2. Build the claim register from `ms.qmd` (including its 25 `TODO` markers); resolve the
   `.bib` to OpenAlex; tag existing entries with claims, starting from
   `ms/ZoteroAnnotations/`. **Checkpoint:** show JV the register and the coverage that the
   existing library already gives.
3. Keyword searches + chaining per block; screening → `candidates.csv`.
   **Checkpoint:** JV skims the top ~30 per block and vetoes or adds items.
4. Full-text pass (OA only) on the shortlist → evidence cards, claim matrix, DOI list.
5. Fill M6 (data/software citations) directly and close those `[TODO cite]` markers on a
   branch.
6. Commit each step on the task branch; record the search in `notebook.md` if it changes
   the manuscript's framing.

The blocks are independent, so steps 3–4 could run as parallel sub-agents (one per block)
if JV wants it faster. By default I'll run them in sequence in one session.

## 6. Known limitations

- OpenAlex abstracts are missing for part of the publisher catalogue; title-only screening is
  weaker, and those papers rely more on chaining.
- Grey and Scandinavian-language literature is out of scope by choice.
- Full text is OA only. New paywalled papers get abstract-level cards; for papers already
  in the `.bib`, JV's highlights stand in. If a paywalled paper turns out to be essential,
  JV can highlight it in Zotero and re-export its annotations at the step-3 checkpoint.
- Quotes from JV's highlights are not re-verified against the PDF.
- Relevance judgements are mine. The checkpoints in §5 are where JV corrects them.

## Questions for JV before starting (answered 2026-11-17; folded into the plan above)

1. **Aims/RQs.** Are RQ1–RQ3 above right? And is the paper primarily ecological (fate of
   Norway's bogs) with the transfer evaluation as support, or are the two equal
   contributions? `ms/PLAN.md` ("Audience and stance") already leans ecological: lead with
   the model-free evidence, then the model, then how far it can be trusted. This sets how
   much of the search goes to blocks M3/M4/D2 versus I1–I5/D3/D8.

> Its really mainly RQ1, with RQ2 as an alternative leg to stand on, and RQ3 just as a validating step (not a RQ).

2. **Target journal.** GCB or PLOS ONE (or other)? `ms.qmd` still has this as a
   `TODO(JV)`, with GCB as the working guess. This sets the reference budget (roughly
   60–90 vs no hard limit) and how much methods literature belongs in the main text.

> GCB is the target journal, and the reference budget is 60–90. The methods literature should be in the main text, but not too much.

3. **Scope.** Include Scandinavian-language literature and grey reports via your input?
   Any time window?

> Not necessary.

4. **Full text and highlights.** I can read `ZoteroLibrary/` and `ms/ZoteroAnnotations/`
   directly. Is it fine to use both? Do you want the highlighted-PDF output (needs
   `pip install pymupdf`), or are the evidence cards enough?

> Use only ZoteroAnnotations, do not read ZoteroLibrary. I do not need the highlighted-PDF output, the evidence cards are enough.

5. **Email-bearing services.** Will you create the OpenAlex API key yourself, or may I
   request it with your email? And may your email go in Unpaywall requests (optional;
   skipped otherwise)?

> I have added the OpenAlex API key to my .Renviron file, and I do not want my email to go in Unpaywall requests.
