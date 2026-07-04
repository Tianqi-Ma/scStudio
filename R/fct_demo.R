#' Demo datasets for first-time users
#'
#' Two ways to get started without your own data:
#' 1. A tiny **bundled** example shipped inside the package -- loads instantly,
#'    fully offline, ideal for clicking through the whole pipeline in seconds.
#' 2. A curated list of **public** datasets that can be fetched over the network
#'    at runtime (10x-hosted HDF5 files load directly).
#'
#' @name fct_demo
#' @keywords internal
NULL

#' Path to the bundled tiny demo counts matrix (.rds)
#'
#' A small synthetic 500-gene x 300-cell sparse counts matrix with three
#' structured cell groups and a few `MT-` genes so every step (QC included) has
#' something to show. Not biologically real -- it exists so new users can try the
#' interface immediately, offline.
#'
#' @return File path (character), or "" if not found.
#' @keywords internal
demo_bundled_path <- function() {
  app_sys("extdata", "demo_pbmc_small.rds")
}

#' Curated public demo datasets (fetched at runtime)
#'
#' Direct-download files that the app can load. The 10x-hosted HDF5 files load
#' via the "10x .h5" path. URLs point to well-known public resources; they are
#' external and may change -- verify before relying on them.
#'
#' @return data.frame with columns: id, name, cells, format, url, description.
#' @keywords internal
demo_datasets <- function() {
  data.frame(
    id = c("bundled", "pbmc1k_v3", "pbmc5k_v3"),
    name = c(
      "Bundled tiny example (instant, offline)",
      "10x PBMC 1k (v3)",
      "10x PBMC 5k (v3)"
    ),
    cells = c("~300", "~1,000", "~5,000"),
    format = c("rds", "h5", "h5"),
    url = c(
      "",  # bundled -> uses demo_bundled_path()
      "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5",
      "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_5k_v3/pbmc_5k_v3_filtered_feature_bc_matrix.h5"
    ),
    description = c(
      "Synthetic 500x300 matrix with 3 cell groups. Loads instantly, no download.",
      "Real human PBMCs, ~1k cells. Small download, good first real test.",
      "Real human PBMCs, ~5k cells. Larger, more realistic structure."
    ),
    stringsAsFactors = FALSE
  )
}

#' Fetch a demo dataset to a local temp file and return its path + format
#'
#' @param id Dataset id from [demo_datasets()].
#' @return list(path, format) or stops with an informative error.
#' @keywords internal
fetch_demo <- function(id) {
  ds <- demo_datasets()
  row <- ds[ds$id == id, , drop = FALSE]
  if (nrow(row) == 0) stop("Unknown demo dataset: ", id)

  if (id == "bundled") {
    p <- demo_bundled_path()
    if (!nzchar(p) || !file.exists(p)) stop("Bundled demo not found in package.")
    return(list(path = p, format = "rds"))
  }

  ext <- switch(row$format, h5 = ".h5", rds = ".rds", ".dat")
  dest <- tempfile(fileext = ext)
  ok <- tryCatch(
    utils::download.file(row$url, dest, mode = "wb", quiet = TRUE),
    error = function(e) 1
  )
  if (!identical(ok, 0L) || !file.exists(dest) || file.size(dest) == 0) {
    stop("Could not download demo (need internet). URL: ", row$url)
  }
  list(path = dest, format = row$format)
}
