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
  [OMICstudio](https://github.com/Tianqi-Ma/OMICstudio) sibling repo.

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
