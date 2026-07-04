#' Object helpers shared across modules
#'
#' Thin, defensive wrappers around Seurat/SeuratObject so modules can query the
#' working object without every module repeating null-checks and class handling.
#'
#' @name fct_object
#' @keywords internal
NULL

#' Dimensions (cells x genes) of the working object
#' @param obj A Seurat object.
#' @return list(cells, genes); NA if unknown.
#' @keywords internal
obj_dims <- function(obj) {
  if (is.null(obj)) return(list(cells = NA_integer_, genes = NA_integer_))
  cells <- tryCatch(ncol(obj), error = function(e) NA_integer_)
  genes <- tryCatch(nrow(obj), error = function(e) NA_integer_)
  list(cells = cells, genes = genes)
}

#' Cell metadata as a data.frame
#' @param obj A Seurat object.
#' @return data.frame (0-col if unavailable).
#' @keywords internal
obj_meta <- function(obj) {
  if (is.null(obj)) return(data.frame())
  md <- tryCatch(obj[[]], error = function(e) NULL)
  if (is.null(md)) md <- tryCatch(obj@meta.data, error = function(e) data.frame())
  if (is.null(md)) md <- data.frame()
  md
}

#' Names of metadata columns
#' @param obj A Seurat object.
#' @keywords internal
obj_meta_cols <- function(obj) {
  colnames(obj_meta(obj))
}

#' Available dimensional reductions (e.g. pca, umap, tsne, harmony)
#' @param obj A Seurat object.
#' @keywords internal
obj_reductions <- function(obj) {
  if (is.null(obj)) return(character(0))
  tryCatch(names(obj@reductions), error = function(e) character(0))
}

#' Does the object have a given reduction?
#' @param obj Seurat object; @param red reduction name.
#' @keywords internal
has_reduction <- function(obj, red) {
  red %in% obj_reductions(obj)
}

#' Extract a 2D embedding + metadata as a tidy data.frame for plotting
#'
#' @param obj A Seurat object.
#' @param reduction Reduction name (e.g. "umap", "tsne", "pca").
#' @param dims Which two dims to take. Default c(1, 2).
#' @param color_by Optional metadata column to attach as `color`.
#' @return data.frame with columns dim1, dim2, cell, and (optionally) color.
#' @keywords internal
embedding_df <- function(obj, reduction, dims = c(1, 2), color_by = NULL) {
  emb <- tryCatch(
    SeuratObject::Embeddings(obj, reduction = reduction),
    error = function(e) NULL
  )
  if (is.null(emb)) stop("Reduction '", reduction, "' not found in object.")
  df <- data.frame(
    dim1 = emb[, dims[1]],
    dim2 = emb[, dims[2]],
    cell = rownames(emb),
    stringsAsFactors = FALSE
  )
  if (!is.null(color_by)) {
    md <- obj_meta(obj)
    if (color_by %in% colnames(md)) df$color <- md[df$cell, color_by]
  }
  df
}

#' Coerce an uploaded object/matrix into a Seurat object
#'
#' Accepts: a Seurat object (updated if from an old version), a
#' SingleCellExperiment (converted), or a matrix / dgCMatrix / data.frame of
#' counts (genes x cells).
#'
#' @param x The loaded R object.
#' @param project Project name for a freshly created object.
#' @return A Seurat object.
#' @keywords internal
as_seurat <- function(x, project = "scStudio") {
  if (!has_pkg("Seurat") || !has_pkg("SeuratObject")) {
    stop("Package 'Seurat' is required to build the working object.")
  }
  cls <- class(x)[1]

  # Already a Seurat object -> update if needed
  if (methods::is(x, "Seurat")) {
    return(tryCatch(SeuratObject::UpdateSeuratObject(x), error = function(e) x))
  }

  # SingleCellExperiment -> Seurat
  if (methods::is(x, "SingleCellExperiment")) {
    if (!has_pkg("SingleCellExperiment")) {
      stop("Package 'SingleCellExperiment' is required to convert this object.")
    }
    return(Seurat::as.Seurat(x, counts = "counts", data = NULL))
  }

  # A plain counts matrix / sparse matrix / data.frame (genes x cells)
  if (methods::is(x, "dgCMatrix") || is.matrix(x) || is.data.frame(x)) {
    m <- x
    if (is.data.frame(m)) m <- as.matrix(m)
    return(Seurat::CreateSeuratObject(counts = m, project = project))
  }

  stop("Unsupported object of class '", cls,
       "'. Expected Seurat, SingleCellExperiment, or a counts matrix/table.")
}

#' Read a table (csv/tsv) of counts into a matrix (genes x cells)
#' @param path File path; @param sep field separator.
#' @keywords internal
read_counts_table <- function(path, sep = "\t") {
  df <- utils::read.delim(path, sep = sep, row.names = 1, check.names = FALSE)
  as.matrix(df)
}
