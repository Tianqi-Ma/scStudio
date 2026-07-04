# scStudio 0.1.0

First scaffold release.

## Added
- Localhost-first Shiny app launched with a single command, `run_app()`, that
  opens in the browser and computes on the user's own machine.
- Twelve analysis modules following a uniform pattern (beginner explainer,
  method choice, adjustable thresholds, run, result summary, interactive preview
  with hover): Import, QC, Doublet removal, Normalize, Features/PCA, Integrate,
  Cluster, Embed, Markers, Annotate, Visualize, Export.
- Modern-practice methods: MAD-based adaptive QC, scDblFinder, LogNormalize /
  SCTransform, Harmony / CCA / RPCA integration, Leiden clustering, UMAP / t-SNE
  embedding, SingleR / Azimuth annotation.
- Cell-type hover explanations via a built-in plain-language dictionary.
- **Demo data**: a bundled tiny example (instant, offline) plus curated public
  datasets (10x PBMC 1k/5k) and a "from URL" loader.
- Three distribution tiers: `remotes::install_github` + `run_app()`; a Docker
  image (all dependencies baked in); and a planned desktop installer.
- Reproducibility log exportable as an R script.
- Bilingual README (English + 简体中文).

## Known limitations
- Not yet validated end-to-end against a live Seurat/Bioconductor install with
  real data; treat as a working scaffold.
- Heavy analysis packages are in `Suggests` and checked at runtime.
- `.h5ad` export is best-effort (SeuratDisk); `.rds` is primary.
- Pseudobulk differential expression, trajectory, and GRN modules are planned
  for a later release.
