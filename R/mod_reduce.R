#' Module 5: Feature selection & PCA
#'
#' Pick highly variable genes (HVGs), scale them, and run PCA to compress the data
#' into a handful of informative components used by every later step.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_reduce
NULL

#' @rdname mod_reduce
#' @keywords internal
mod_reduce_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = "Feature selection & PCA",
    what = "Select the most informative genes (highly variable features), then run
            <b>PCA</b> to summarise them as a few principal components.",
    why  = "Most genes vary little between cells and just add noise. Focusing on
            highly variable genes and compressing them with PCA makes clustering
            and embedding faster and cleaner.",
    how  = "Choose how variable genes are ranked, how many to keep, and how many
            principal components to compute. The elbow plot on the right shows how
            much variation each component captures &mdash; keep components before the
            curve flattens.",
    example = "From 20,000 genes you keep ~2,000 variable ones, then summarise them
               as 50 PCs; the first ~20 usually carry the real structure."
  )
  controls <- shiny::tagList(
    label_with_help("HVG method",
                    "How variable genes are ranked. vst = variance-stabilizing (recommended); mvp/dispersion = mean-variance/dispersion based."),
    shiny::selectInput(ns("hvg_method"), NULL,
                       choices = c("vst" = "vst", "mvp" = "mvp",
                                   "dispersion" = "dispersion"),
                       selected = "vst"),
    label_with_help("Number of variable genes",
                    "How many highly variable genes to keep. 2,000 is a common default."),
    shiny::sliderInput(ns("n_hvg"), NULL, min = 500, max = 5000, value = 2000,
                       step = 100),
    label_with_help("Number of principal components",
                    "How many PCs to compute. 50 is a common default; you rarely use them all downstream."),
    shiny::numericInput(ns("npcs"), NULL, value = 50, min = 2, max = 200, step = 1),
    run_button(ns("run"), "Select features & run PCA")
  )
  step_container(explainer, controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview")))
}

#' @rdname mod_reduce
#' @keywords internal
mod_reduce_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(done = FALSE, n_hvg = NA, npcs = NA, stdev = NULL)

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs("Seurat", "Feature selection & PCA")) return(NULL)
      n_hvg <- input$n_hvg
      npcs <- input$npcs
      hvg_method <- input$hvg_method
      obj <- with_progress_notify({
        reduce_obj(rv$obj, n_hvg = n_hvg, npcs = npcs, hvg_method = hvg_method)
      }, message = "Selecting features and running PCA...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$done <- TRUE
      res$n_hvg <- length(tryCatch(Seurat::VariableFeatures(obj),
                                   error = function(e) character(0)))
      res$stdev <- tryCatch(obj@reductions$pca@stdev, error = function(e) NULL)
      res$npcs <- length(res$stdev)
      log_step(log_rv, "Feature selection & PCA",
               params = list(hvg_method = hvg_method, n_hvg = n_hvg, npcs = npcs),
               code = sprintf(paste0(
                 "obj <- Seurat::FindVariableFeatures(obj, selection.method = '%s', nfeatures = %d)\n",
                 "obj <- Seurat::ScaleData(obj)\n",
                 "obj <- Seurat::RunPCA(obj, npcs = %d)"),
                 hvg_method, n_hvg, npcs))
      shiny::showNotification(sprintf("PCA done: %d HVGs, %d PCs.",
                                      res$n_hvg, res$npcs), type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) return(shiny::div(class = "scstudio-placeholder",
                                               "Set parameters and click Select features & run PCA."))
      bslib::layout_columns(
        col_widths = c(6, 6),
        stat_tile("Variable genes", format(res$n_hvg, big.mark = ",")),
        stat_tile("Principal components", format(res$npcs, big.mark = ","))
      )
    })

    output$preview <- render_preview_plot(function() {
      stdev <- res$stdev
      shiny::req(stdev)
      df <- data.frame(PC = seq_along(stdev), stdev = stdev)
      df$text <- sprintf("PC%d\nstdev=%.3f", df$PC, df$stdev)
      ggplot2::ggplot(df, ggplot2::aes(x = PC, y = stdev, text = text)) +
        ggplot2::geom_line(colour = "#7d8b8f", linewidth = 0.4) +
        ggplot2::geom_point(colour = scstudio_palette(1), size = 1.4) +
        ggplot2::labs(x = "Principal component", y = "Standard deviation",
                      title = "PCA elbow plot") +
        scstudio_theme()
    })
  })
}
