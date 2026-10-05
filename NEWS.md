# scStudio 0.3.1 (unreleased) — audit 2026-10

## Figures, colours and navigation (feature/plot-style)

- **One figure style.** Every figure is ggplot2 in one theme; Seurat's DimPlot /
  FeaturePlot / VlnPlot / DotPlot / DoHeatmap are replaced by `sc_dim_plot()`,
  `sc_feature_plot()`, `sc_violin_plot()`, `sc_dot_plot()` and `sc_heatmap()`
  (cell maps with corner axes and in-place labels; expression scaled per gene;
  dot and heatmap on a diverging z-score scale).
- **A level keeps its colour.** A colour registry (`group_colors()`, stored in
  `obj@misc`) gives every cluster, cell type, sample and condition the same
  colour in every step; cell types inherit their cluster's colour; fixed
  meanings (kept / removed, mutant / WT, low / high, cell-cycle phase) keep
  fixed colours.
- `compose_grid()` lays out multi-panel figures with aligned panels, so the KM
  risk table, the PCA / HVG pair and the annotation map + bar now draw on every
  installation (they needed patchwork >= 1.2 before).
- **Step navigation moved to the bottom of each step**: previous on the left,
  next on the right (highlighted once the step has run). The top bar keeps
  the progress chip only.
- **Animations fill the empty result area** instead of a 620 px strip.

## New analysis steps (feature/analysis-modules)

- **Pseudobulk DE** (new phase *Compare conditions*): counts of each cell type
  summed per sample, tested across samples with edgeR quasi-likelihood
  (default), limma-voom or DESeq2; optional pairing / batch covariate; FDR
  within each cell type; cell types with < 2 samples per group are skipped and
  listed. Overview, volcano, per-sample cell counts, full table.
- **Differential abundance**: cell-type proportions per sample, propeller method
  (logit + limma robust eBayes; t test for two conditions, F test for more; BH)
  written with limma, so speckle is not needed.
- Both steps share the sample / condition / cell-type selectors and accept a
  **sample sheet** (CSV/TSV, one row per sample) to add a condition the object
  lacks; a condition that varies within a sample is refused.
- Runnable log code for both, tested to reproduce the app's p-values;
  explainer animations; a "Condition comparisons" report section.


## Added
- **Step animations.** Every step's empty preview plays a short looping scene of
  what the step does (QC outliers being flagged, a PCA cloud flattening onto its
  components, cells pulling into a UMAP, a Kaplan-Meier split...). Plain Canvas,
  no libraries, works offline, follows EN / 中 and light / dark, pauses when off
  screen, honours reduced motion. Engine `inst/app/www/explain.js`, scenes in
  `explain-sc.js`; checked by `node tools/check_explain.js`.
- **Pipeline state** (`R/fct_state.R`): every registry entry declares the steps it
  consumes; re-running a step marks finished downstream steps *stale* (amber dot
  and a banner); importing a dataset starts a new epoch that clears the old
  results, statuses and log; QC and doublet removal re-run from their own input
  instead of re-filtering filtered cells.
- **Every single-cell figure is downloadable** (PNG / JPEG / PDF, DPI) — the same
  closure that draws it on screen (`render_step_plot()`).
- **Survival layer**: 95% bands and a number-at-risk table on every KM plot, the
  HR (95% CI) beside it, Cox HR per SD as the primary result for continuous
  variables, `cox.zph` check, multivariable Cox with an events-per-variable
  check, small-sample warnings, and a selection-adjusted p-value for the
  optimal cut-point (maxstat or the Lausen–Schumacher approximation).
- **Patient-level composition**: filter the cells (e.g. tumour, baseline), choose
  the denominator (fraction *of T cells*), drop patients with too few cells —
  then join, then cut.
- `AGENTS.md` / `CLAUDE.md`: the design contract every model follows;
  `tools/sync_mirror.sh` keeps the two repositories identical where they share
  code.

## Changed
- The app follows the OS light/dark preference; the language choice is remembered
  and a first visit follows the browser language.
- System font stack only (Google Fonts was render-blocking where unreachable).
- The "next step" link moved into the top bar (the floating chip covered the
  figure-export row); the always-"not set up" Python indicator is gone.
- The export menu offers what the active pipeline can export (WES: MAF, mutation
  table, sample / gene summaries; plus the clinical cohort).
- Reproducibility log: a re-run replaces its entry; every logged command is
  runnable base R / Seurat / scop / maftools code with the parameters used (no
  package internals), and replays to the app's results.
- scop wrappers rewritten against scop 0.9.2 (`fct_scop.R`); dead wrappers removed.

## Fixed
- **Doublets** are called per sample (`samples =`), with a seed.
- **SCT**: PCA now runs on the Pearson residuals (no re-scaling); switching back
  to LogNormalize no longer normalises the SCT assay.
- **Integration**: CCA/RPCA results are selectable downstream; "none" and a new
  PCA drop stale Harmony / integrated reductions; the preview compares the
  uncorrected and corrected spaces with a batch-mixing entropy.
- **Clustering**: Leiden backend detected (leidenbase / igraph), k is applied, the
  active resolution is explicit.
- **Markers** filtered on adjusted p before ranking; dot plot as described.
- **Import**: SCE / v4-integrated / multi-layer / multimodal h5 objects are
  standardised to an `RNA` assay.
- **QC**: per-sample MAD option, MAD = 0 guarded, haemoglobin pattern no longer
  matches HBEGF / HBS1L, Ensembl ids detected; violins with thresholds.
- **Normalisation preview** shows the depth effect before vs after.
- **Enrichment / GSEA, RNA velocity, trajectories (Slingshot / Monocle / PAGA /
  Palantir / WOT), CellChat, LIANA, copykat, signature scoring** all called with
  their current arguments; copykat low-confidence and undefined calls no longer
  become "normal"; mascarade outlines draw.
- **Survival**: tied tertiles fail loudly, a single-group split is refused,
  duplicated patients are refused, the stale event coding no longer zeroes all
  events, and the shared cohort is never overwritten by an analysis subset.

## Validation
- Tests: 608 passing (R 4.2.3 + Seurat 4.3); the scop, CellChat, LIANA,
  copykat, SingleR, Azimuth and scDblFinder paths are checked against source but
  were not executed. See `docs/AUDIT_2026-10.md`.

---

# scStudio 0.3.0

## Added
- **Clinical & survival step** (step 18), backed by a new `R/fct_survival.R`:
  Kaplan-Meier curves drawn in ggplot2, the log-rank test, median survival per
  arm, and univariable Cox regression. Stratify by a clinical column or by
  **per-sample cell-type composition** (median / tertile / log-rank-optimal
  cutpoint). The maths uses `survival`, which ships with R, so nothing needs to
  be compiled. The cohort is stored on the shared hub (`rv$clinical`).
- **scop as the engine.** Plotting (`CellDimPlot`, `GroupHeatmap`,
  `DynamicHeatmap`, `VolcanoPlot`, ...) and compute (`Standard_SCP`,
  `Integration_SCP`, `RunDEtest`, `RunEnrichment`, `RunGSEA`, `RunSlingshot`,
  Monocle2/3, PAGA, Palantir, WOT, scVelo, dynamic features) behind wrappers that
  degrade to a readable message when scop is absent.
- Eight further steps: enrichment/GSEA, trajectory, RNA velocity, dynamic
  features, cell cycle & signatures, cell communication (LIANA/CellChat),
  malignancy/CNV (CopyKAT), and the narrated report. Twenty-one steps across
  seven phases.
- Mascarade cluster outlining on the embedding and visualization steps.
- **Import overview.** The Import step previews the raw data before any
  filtering: per-cell QC violins, top-20 expressed genes, a counts-vs-genes
  scatter, the cell metadata table and a slice of the counts matrix. Species is
  detected from gene-name casing so mouse data no longer reads 0% mitochondrial.
- **Shared export menu** in the top bar, available on every step: the object
  (.rds), cell metadata (.csv), counts matrix (.rds) and embeddings (.csv).
- **Bundled pbmc3k.** The real 2,700-cell 10x PBMC dataset ships inside the
  package, so the recommended demo loads instantly and offline — no download and
  no SeuratData.

## Fixed
- `DESCRIPTION` declared `bslib (>= 0.5.0)` while the UI uses
  `input_dark_mode()` (bslib 0.6.0) and `input_task_button()` (bslib 0.7.0);
  installing against an older bslib produced an app that crashed on startup.
  Now `bslib (>= 0.7.0)`.
- Interactive plots asked `ggplotly()` for a `text` tooltip that no layer
  provides, so hovering showed an empty box; and a multi-gene FeaturePlot (a
  patchwork) was handed to `ggplotly()`, which cannot render it. Interactivity is
  now used only where it works, and shows the real values.
- Static assets 404'd without `addResourcePath()`; plot failures surfaced as an
  empty "Plot error:" or an opaque `[object Object]`.
- Visualize / Report / Export never marked themselves complete, so the step
  navigator could not reach the end of the pipeline.
- The Export step could error on its own filename when the format input had not
  been initialised.
- Missing styles for the sidebar's empty-state and dataset readout; the Google
  font is now optional rather than a hard requirement for an offline app.

## Testing
- New guards: every step's UI is built, every module's outputs are forced to
  evaluate (with and without a working object), and every function is checked
  for calls to names that do not resolve. The survival layer has its own unit
  tests. This is the gap that let a stale helper name reach main in the
  [OmicOne](https://github.com/Tianqi-Ma/OmicOne) sibling repo.

## Known limitations
- **Not validated end-to-end** against a live Seurat/scop install with real
  data; some scop argument names still need checking.
- The "optimal" survival cutpoint is exploratory: chosen to maximise separation,
  so its p-value is optimistic.
- `.h5ad` export is best-effort (SeuratDisk); `.rds` is primary.

---

# scStudio 0.2.0

- Dark theme by default with a light toggle; grouped left step navigator with
  per-step status; `input_task_button()` so a running step is visibly running;
  plot-first layout (narrow control rail, large canvas, collapsed explainer).
- Single-language interface (English or 中文, never both) via `data-en`/`data-zh`
  attributes swapped client-side — no server round-trip.
- Startup animation: scattered cells coalescing into UMAP-like clusters.

# scStudio 0.1.0

- First scaffold: localhost-first Shiny app launched with `run_app()`, twelve
  analysis modules following a uniform pattern, MAD-based QC, scDblFinder,
  Harmony, Leiden, UMAP/t-SNE, SingleR/Azimuth, a reproducibility log exportable
  as an R script, a Docker image, and a bilingual README.
