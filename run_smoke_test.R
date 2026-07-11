#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# scStudio end-to-end smoke test (headless — no browser, no clicking)
#
# Runs the real compute chain on the bundled tiny demo dataset:
#   import -> QC metrics -> MAD filter -> doublets -> normalize ->
#   features/PCA -> cluster (Louvain) -> UMAP -> markers
# and prints [OK]/[FAIL] per step. Use it to verify your machine has the
# dependencies working, without launching the app.
#
# Usage (from the repository root, i.e. the folder containing DESCRIPTION):
#   Rscript run_smoke_test.R
# ---------------------------------------------------------------------------

ok <- 0L; fail <- 0L
step <- function(label, expr) {
  res <- tryCatch({ force(expr); TRUE },
                  error = function(e) { cat(sprintf("  [FAIL] %-26s %s\n", label, conditionMessage(e))); FALSE })
  if (isTRUE(res)) { cat(sprintf("  [OK]   %s\n", label)); ok <<- ok + 1L } else { fail <<- fail + 1L }
  invisible(res)
}

cat("== scStudio smoke test ==\n")

# --- load the package sources without building (works under locked-down policy) ---
if (!requireNamespace("pkgload", quietly = TRUE)) {
  install.packages("pkgload", repos = "https://cloud.r-project.org")
}
suppressMessages(pkgload::load_all(".", quiet = TRUE, helpers = FALSE, attach_testthat = FALSE))
ns <- asNamespace("scStudio")

# --- dependency check ---
need <- c("Seurat", "SeuratObject", "Matrix")
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  cat("Missing required packages:", paste(missing, collapse = ", "),
      "\nInstall them first (see README), then re-run.\n")
  quit(status = 1)
}

# --- the pipeline ---
counts <- readRDS(ns$demo_bundled_path())
obj <- NULL; mk <- NULL

step("import (build object)",      { obj <<- ns$as_seurat(counts) })
step("QC metrics",                 { obj <<- ns$qc_add_metrics(obj, "human") })
step("QC MAD filter",              { keep <- ns$qc_mad_keep(obj, 5, 3); obj <<- obj[, keep] })
step("normalize (LogNormalize)",   { obj <<- ns$normalize_obj(obj, "LogNormalize") })
step("features + PCA",             { obj <<- ns$reduce_obj(obj, n_hvg = 300, npcs = 20) })
step("cluster (Louvain)",          { obj <<- ns$cluster_obj(obj, "pca", dims = 15,
                                                            resolutions = 0.5, algorithm = 1) })
step("embed (UMAP)",               { obj <<- ns$embed_obj(obj, "umap", "pca", dims = 15) })
step("marker genes",               { mk  <<- ns$markers_obj(obj, "wilcox", 0.25, 0.1, TRUE) })

# optional steps that need extra packages / internet — reported, not counted as failures
if (requireNamespace("scDblFinder", quietly = TRUE)) {
  step("doublets (scDblFinder)",   { invisible(ns$run_doublets(obj, "scDblFinder")) })
} else cat("  [skip] doublets — install 'scDblFinder' to test this step\n")

cat(sprintf("\n== Result: %d passed, %d failed ==\n", ok, fail))
if (!is.null(obj)) {
  d <- ns$obj_dims(obj)
  cat(sprintf("Final object: %s cells x %s genes; clusters: %s; markers found: %s\n",
              d$cells, d$genes,
              tryCatch(length(unique(SeuratObject::Idents(obj))), error = function(e) NA),
              if (!is.null(mk)) nrow(mk) else NA))
}
quit(status = if (fail > 0) 1 else 0)
