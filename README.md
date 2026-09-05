# scStudio

> 🌐 **Languages:** **English** · [中文版 README](README.zh-CN.md)

**Interactive single-cell RNA-seq analysis on your own machine.** Launch a browser
interface with one command, upload a count matrix, and interactively run a modern
scRNA-seq pipeline — QC, doublet removal, normalization, dimensionality reduction,
integration, clustering, embedding, markers, annotation, enrichment, trajectory,
RNA velocity, cell-cell communication, CNV inference and survival analysis —
using your own computer's CPU and RAM. No cloud server required.

Every step is designed for beginners *and* experts:

- 💡 a plain-language **"What is this step?"** explainer (with a worked example)
- 🔧 a **method choice** (each step offers alternatives — e.g. UMAP *or* t-SNE)
- 🎚️ **adjustable thresholds** with sensible defaults
- 📊 a **result summary** and a large **preview plot**
- 🧾 a **reproducibility log** you can export as an R script or a narrated report
- 🌗 dark / light themes, English / 中文 — switchable from the top bar

> **Status (v0.3.0):** 21 steps across 7 phases, all wired and statically
> validated (parse, UI build, output evaluation, tests, install). It has **not
> yet been run end-to-end against a live Seurat/scop install** — see
> [Caveats](#caveats).
>
> scStudio is also the single-cell half of
> [OMICstudio](https://github.com/Tianqi-Ma/OMICstudio), a multi-omics suite
> built from the same code; shared changes are mirrored between the two.

---

## Why "localhost-first"?

Real scRNA-seq analysis (Seurat/Bioconductor) needs native compute and real RAM.
Browser-only (WASM) apps can't run it. So scStudio runs a **local server + browser
UI**: the interface is a web page, but all computation happens on *your* machine.
This is the same model as `cellxgene launch`.

---

## Three ways to run it

Pick the tier that matches how much you want to install.

### A. You already have R (lightest)
```r
# install.packages("remotes")
remotes::install_github("Tianqi-Ma/scStudio")
scStudio::run_app()   # opens your browser automatically
```
Needs R ≥ 4.1 with **shiny ≥ 1.7.4** and **bslib ≥ 0.7.0** (shiny ≥ 1.8.1 is
recommended), plus the heavy analysis packages for the steps you actually run.
Smallest download, one command.

### B. No R, no dependencies → Docker (recommended for most users)
Everything (R + Seurat + Bioconductor + the app) is baked into one image. You only
need [Docker](https://www.docker.com/products/docker-desktop/).
```bash
docker run --rm -p 3838:3838 -m 16g tianqima/scstudio
# then open http://localhost:3838 in your browser
```
Upload your data through the browser (no volume mounting needed). Give Docker
enough memory (`-m 16g`, and raise Docker Desktop's memory limit for large data).

> Build the image yourself from this repo: `docker build -t scstudio .`

### C. Non-technical, double-click (planned)
A desktop installer (Tauri/Electron/electricShine) that bundles R and the app so
users just double-click — no terminal, no Docker. Planned for a later release.

---

## Try it instantly (no data needed)

On the **Import** step, keep the source on **Demo data** and click *Load demo
data*. Three demos are offered:

| Demo | What it is | Needs |
|---|---|---|
| **PBMC 3k** (default) | the classic 2,700-cell 10x PBMC dataset, **bundled in the package** | nothing — instant, offline |
| **Pancreas** | mouse pancreas with lineage structure and spliced/unspliced layers, for trajectory / velocity | `scop` installed |
| **Tiny example** | 500 × 300 synthetic matrix for a quick UI check | nothing |

You can also paste a URL to any `.rds` / `.h5` file.

## What you can upload

- **RDS** — a saved `Seurat` or `SingleCellExperiment` object (old Seurat objects
  are updated automatically)
- **10x `.h5`** — Cell Ranger HDF5
- **Counts table** — CSV/TSV, genes in rows, cells in columns
- **Clinical table** — CSV/TSV, one row per patient (id, follow-up time, outcome)

Browser uploads are capped high (5 GB by default; change with
`run_app(max_upload_mb = ...)`).

---

## The analysis steps

Seven phases, twenty-one steps. Bold = the default method.

| Phase | # | Step | Methods offered |
|---|---|------|-----------------|
| **Data & QC** | 1 | Import | RDS / 10x .h5 / table / URL / demo |
| | 2 | Quality control | **MAD (adaptive)** / manual; mito, ribo, hemoglobin, dissociation % |
| | 3 | Doublets | **scDblFinder** / DoubletFinder |
| **Preprocess** | 4 | Normalize | **LogNormalize** / SCTransform |
| | 5 | Features / PCA | HVG **vst** / mvp / dispersion |
| | 6 | Integrate | **none** / Harmony / CCA / RPCA / scVI / scanorama / BBKNN |
| **Structure** | 7 | Cluster | **Leiden** / Louvain |
| | 8 | Embed | **UMAP** / t-SNE / PaCMAP / PHATE |
| **Identity** | 9 | Markers | **wilcox** / roc / MAST |
| | 10 | Annotate | **manual** / SingleR / Azimuth / scop KNN prediction |
| | 11 | Enrichment / GSEA | ORA + GSEA (clusterProfiler via scop) |
| **Trajectory** | 12 | Trajectory | **Slingshot** / Monocle2 / Monocle3 / PAGA / Palantir |
| | 13 | RNA velocity | scVelo (steady-state / stochastic / dynamical) |
| | 14 | Dynamic features | scop dynamic features + heatmap |
| **Advanced** | 15 | Cell cycle & signatures | Seurat cell cycle, UCell / AddModuleScore |
| | 16 | Cell communication | LIANA / CellChat |
| | 17 | Malignant / CNV | CopyKAT, stemness |
| | 18 | **Clinical & survival** | Kaplan-Meier + log-rank + univariable Cox; stratify by a clinical column or by per-sample cell-type composition (median / tertile / optimal cutpoint) |
| **Output** | 19 | Visualize | UMAP / violin / dotplot / feature / heatmap |
| | 20 | Report | narrated HTML report from the reproducibility log |
| | 21 | Export | .rds / .h5ad (best-effort) / figures / R script |

A **⤓ Export** menu in the top bar is available on *every* step (object,
metadata, counts, embeddings) so you never have to reach the last step to get
your data out.

Methods and defaults follow current (2023–2025) best practice: MAD-based QC,
scDblFinder, Leiden, Harmony, scop/SCP plotting throughout.

---

## Caveats

- **Not yet validated end-to-end.** Every file parses, every step's UI builds,
  every module's outputs evaluate, and the package installs — but the pipeline
  has **not been run against a live Seurat/scop install with real data**. A few
  scop argument names still need checking on first real run.
- **`scop` is the plotting and compute engine** and is a GitHub package; if your
  machine blocks source builds, use the Docker image, which bakes in scop and a
  pre-built conda environment for the Python-backed steps (scVelo, PAGA,
  Palantir).
- **Heavy packages are optional at install time** (they live in `Suggests`). The
  UI loads without them; each compute step checks for what it needs and tells you
  what to install if it's missing. The Docker image bakes them all in.
- **Some steps need internet** the first time (SingleR/Azimuth reference download).
- **`.h5ad` export** in pure R relies on SeuratDisk (best-effort); `.rds` is the
  primary export format.
- **The "optimal" survival cutpoint is exploratory.** It is chosen to maximise
  separation, so its p-value is optimistic — report the median split, or validate
  the cutpoint in an independent cohort.
- **Memory** scales with cell count. <100k cells is comfortable on 16–32 GB;
  larger needs more. The app warns you and offers downsampling.

---

## Development

R package layout (golem-style): `R/app_ui.R` (shell + the `app_steps()` registry),
`R/app_server.R`, `R/run_app.R`, one `R/mod_*.R` per step, compute wrappers in
`R/fct_compute.R`, scop plotting in `R/fct_plots.R`, survival maths in
`R/fct_survival.R`, shared helpers in `R/utils_*.R`. Run locally with:
```r
pkgload::load_all(); run_app()
```

Before committing, the checks that must pass:
```r
testthat::test_dir("tests/testthat")   # includes a guard that every module's
                                       # outputs evaluate, and that no function
                                       # calls a name that does not exist
```

## License

MIT.
