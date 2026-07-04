# scStudio

> 🌐 **Languages:** **English** · [中文版 README](README.zh-CN.md)

**Interactive single-cell RNA-seq analysis on your own machine.** Launch a browser
interface with one command, upload a count matrix, and interactively run a modern
scRNA-seq pipeline — QC, doublet removal, normalization, dimensionality reduction,
integration, clustering, embedding, marker detection, and cell-type annotation —
using your own computer's CPU and RAM. No cloud server required.

Every step is designed for beginners *and* experts:

- 💡 a plain-language **"What is this step?"** explainer (with a worked example)
- 🔧 a **method choice** (each step offers alternatives — e.g. UMAP *or* t-SNE)
- 🎚️ **adjustable thresholds** with sensible defaults
- 📊 a **result summary** and an **interactive preview plot** with **hover details**
- 🧾 a **reproducibility log** you can export as an R script

> **Status:** early scaffold (v0.1). The UI, module structure, and analysis
> wrappers are in place. It has **not yet been run end-to-end against a live
> Seurat/Bioconductor install** — see [Caveats](#caveats).

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
Requires R plus the heavy analysis packages (Seurat, Bioconductor). Smallest
download, one command.

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
data*. A tiny bundled example loads immediately (offline). You can also pick a
real public dataset (10x PBMC 1k/5k, downloaded on first use) or paste a URL to
any `.rds` / `.h5` file.

## What you can upload

- **RDS** — a saved `Seurat` or `SingleCellExperiment` object (old Seurat objects
  are updated automatically)
- **10x `.h5`** — Cell Ranger HDF5
- **Counts table** — CSV/TSV, genes in rows, cells in columns

Browser uploads are capped high (5 GB by default; change with
`run_app(max_upload_mb = ...)`).

---

## The analysis steps

| # | Step | Methods offered | You control |
|---|------|-----------------|-------------|
| 1 | Import | RDS / 10x .h5 / table | format |
| 2 | QC | **MAD (adaptive)** / manual | MAD multiplier or fixed cutoffs, species |
| 3 | Doublet removal | **scDblFinder** / DoubletFinder | flag vs remove, score cutoff |
| 4 | Normalization | **LogNormalize** / SCTransform | scale factor |
| 5 | Features + PCA | HVG **vst** / mvp / dispersion | #HVGs, #PCs |
| 6 | Integration *(optional)* | **none** / Harmony / CCA / RPCA | batch column |
| 7 | Clustering | **Leiden** / Louvain | resolution(s), dims |
| 8 | Embedding | **UMAP** / t-SNE / PaCMAP | neighbors, min_dist / perplexity |
| 9 | Markers | **wilcox** / roc / MAST | logFC, min.pct, top-N |
| 10 | Annotation | **manual** / SingleR / Azimuth | reference, per-cluster labels |
| 11 | Visualize | UMAP / violin / dotplot / feature / heatmap | genes, grouping |
| 12 | Export | .rds / .h5ad (best-effort) / figures / R script | — |

Methods and defaults follow current (2023–2025) best practice (MAD-based QC,
scDblFinder, Leiden, Harmony, pseudobulk DE in a later release, etc.).

---

## Caveats

- **Not yet validated end-to-end.** Built and parse-checked, but not run against a
  live Seurat/Bioconductor install with real data. Treat v0.1 as a working
  scaffold; expect to fix rough edges on first real run.
- **Heavy packages are optional at install time** (they live in `Suggests`). The
  UI loads without them; each compute step checks for what it needs and tells you
  what to install if it's missing. The Docker image bakes them all in.
- **Some steps need internet** the first time (SingleR/Azimuth reference download).
- **`.h5ad` export** in pure R relies on SeuratDisk (best-effort); `.rds` is the
  primary export format.
- **Memory** scales with cell count. <100k cells is comfortable on 16–32 GB;
  larger needs more. The app warns you and offers downsampling.

---

## Development

R package layout (golem-style): `R/app_ui.R`, `R/app_server.R`, `R/run_app.R`,
one `R/mod_*.R` per step, compute wrappers in `R/fct_compute.R`, shared helpers in
`R/utils_*.R` and `R/fct_*.R`. Run locally during development with:
```r
pkgload::load_all(); run_app()
```

## License

MIT.
