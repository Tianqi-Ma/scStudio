#' Demo datasets for first-time users
#'
#' Three ways to get started:
#' 1. **pbmc3k** — the classic 2,700-cell 10x PBMC dataset (via SeuratData);
#'    the recommended demo for a full run-through.
#' 2. **pancreas_sub** — scop's bundled mouse pancreas dataset with lineage
#'    structure and spliced/unspliced layers; use it for trajectory / RNA velocity.
#' 3. **bundled** — a tiny synthetic matrix shipped in the package; loads instantly
#'    offline for a quick UI smoke test.
#'
#' @name fct_demo
#' @keywords internal
NULL

#' Path to the bundled tiny demo counts matrix (.rds)
#' @keywords internal
demo_bundled_path <- function() app_sys("extdata", "demo_pbmc_small.rds")

#' Curated demo dataset catalogue
#' @return data.frame: id, name, cells, format, source, description.
#' @keywords internal
demo_datasets <- function() {
  data.frame(
    id = c("pbmc3k", "pancreas_sub", "bundled"),
    name = c("PBMC 3k (10x, recommended)",
             "Pancreas (scop, for trajectory/velocity)",
             "Tiny example (instant, offline)"),
    cells = c("~2,700", "~1,000", "~300"),
    format = c("seuratdata", "scop", "rds"),
    source = c("SeuratData::pbmc3k", "scop::pancreas_sub", "bundled"),
    description = c(
      "Classic human PBMC dataset. Best for a full walk-through; needs the SeuratData package (downloads once).",
      "Mouse pancreas with lineage structure and spliced/unspliced counts. Ideal for trajectory and RNA velocity; ships with scop.",
      "Synthetic 500x300 matrix. Loads instantly with no download for a quick UI test."
    ),
    stringsAsFactors = FALSE
  )
}

#' Load a demo dataset, returning a ready object or a (path, format) to import
#'
#' @param id Dataset id from [demo_datasets()].
#' @return list(obj=, path=, format=). Exactly one of obj/path is non-NULL.
#' @keywords internal
fetch_demo <- function(id) {
  ds <- demo_datasets()
  if (!id %in% ds$id) stop("Unknown demo dataset: ", id)

  if (id == "bundled") {
    p <- demo_bundled_path()
    if (!nzchar(p) || !file.exists(p)) stop("Bundled demo not found in package.")
    return(list(obj = NULL, path = p, format = "rds"))
  }

  if (id == "pbmc3k") {
    if (!requireNamespace("SeuratData", quietly = TRUE)) {
      stop("The pbmc3k demo needs the 'SeuratData' package. Install it with:\n",
           "  remotes::install_github('satijalab/seurat-data')\n",
           "then: SeuratData::InstallData('pbmc3k')")
    }
    obj <- tryCatch(SeuratData::LoadData("pbmc3k"), error = function(e) {
      stop("Could not load pbmc3k. Run SeuratData::InstallData('pbmc3k') first. (", conditionMessage(e), ")")
    })
    return(list(obj = obj, path = NULL, format = "object"))
  }

  if (id == "pancreas_sub") {
    if (!requireNamespace("scop", quietly = TRUE)) {
      stop("The pancreas demo ships with 'scop'. Install scop first.")
    }
    e <- new.env()
    utils::data("pancreas_sub", package = "scop", envir = e)
    obj <- get("pancreas_sub", envir = e)
    return(list(obj = obj, path = NULL, format = "object"))
  }

  stop("Unhandled demo id: ", id)
}
