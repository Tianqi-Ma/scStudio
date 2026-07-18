#' Shared "export current data" download handlers
#'
#' Wires the four top-bar download buttons (available on every step) to the
#' current working object in `rv$obj`: the Seurat object, the cell metadata, the
#' expression matrix, and the dimensional-reduction embeddings.
#'
#' @param input,output,session Shiny server context.
#' @param rv Shared reactive hub (uses `rv$obj`).
#' @keywords internal
register_exports <- function(input, output, session, rv) {

  stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
  need_obj <- function() {
    if (is.null(rv$obj)) {
      shiny::showNotification("Load data first, then export.", type = "warning")
      return(FALSE)
    }
    TRUE
  }

  # 1. Full Seurat object (.rds)
  output$dl_rds <- shiny::downloadHandler(
    filename = function() paste0("scstudio_object_", stamp(), ".rds"),
    content = function(file) {
      if (!need_obj()) { saveRDS(NULL, file); return() }
      saveRDS(rv$obj, file)
    }
  )

  # 2. Cell metadata (.csv)
  output$dl_meta <- shiny::downloadHandler(
    filename = function() paste0("scstudio_metadata_", stamp(), ".csv"),
    content = function(file) {
      if (!need_obj()) { utils::write.csv(data.frame(), file); return() }
      utils::write.csv(obj_meta(rv$obj), file, row.names = TRUE)
    }
  )

  # 3. Expression matrix — sparse counts as .rds (keeps gene/cell names, compact)
  output$dl_matrix <- shiny::downloadHandler(
    filename = function() paste0("scstudio_counts_", stamp(), ".rds"),
    content = function(file) {
      if (!need_obj()) { saveRDS(NULL, file); return() }
      m <- tryCatch(SeuratObject::LayerData(rv$obj, layer = "counts"),
                    error = function(e) tryCatch(SeuratObject::GetAssayData(rv$obj, slot = "counts"),
                                                 error = function(e2) NULL))
      if (is.null(m)) { shiny::showNotification("No counts matrix found.", type = "error"); saveRDS(NULL, file); return() }
      saveRDS(m, file)
    }
  )

  # 4. Dimensional-reduction embeddings (.csv), all reductions side by side
  output$dl_embed <- shiny::downloadHandler(
    filename = function() paste0("scstudio_embeddings_", stamp(), ".csv"),
    content = function(file) {
      if (!need_obj()) { utils::write.csv(data.frame(), file); return() }
      reds <- obj_reductions(rv$obj)
      if (!length(reds)) {
        shiny::showNotification("No reductions yet (run PCA/UMAP first).", type = "warning")
        utils::write.csv(data.frame(), file); return()
      }
      parts <- lapply(reds, function(r) {
        e <- tryCatch(SeuratObject::Embeddings(rv$obj, reduction = r), error = function(x) NULL)
        if (is.null(e)) return(NULL)
        colnames(e) <- paste0(r, "_", seq_len(ncol(e)))
        as.data.frame(e)
      })
      parts <- Filter(Negate(is.null), parts)
      out <- do.call(cbind, parts)
      utils::write.csv(out, file, row.names = TRUE)
    }
  )
}
