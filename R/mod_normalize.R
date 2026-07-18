#' Module 4: Normalization
#'
#' Put cells on a comparable scale so that differences reflect biology rather than
#' sequencing depth. Supports classic LogNormalize (default) or SCTransform.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_normalize
NULL

#' @rdname mod_normalize
#' @keywords internal
mod_normalize_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = "Normalization",
    what = "Adjust each cell's counts so that cells sequenced to different depths
            become comparable.",
    why  = "Raw counts depend on how deeply each cell was sequenced. Without
            normalization, deeper cells look artificially 'more expressing'; the
            comparisons downstream would reflect depth, not biology.",
    how  = "<b>LogNormalize</b> scales each cell to a common total, then log
            transforms (robust default). <b>SCT</b> models counts with a
            regularized negative binomial and often needs no extra scaling.",
    example = "A cell with 20,000 UMIs and one with 5,000 UMIs are put on the same
               scale so a shared marker reads similarly in both."
  )
  controls <- shiny::tagList(
    label_with_help("Method",
                    "LogNormalize = classic log-scaled counts (recommended default). SCT = variance-stabilizing transform (SCTransform)."),
    shiny::selectInput(ns("method"), NULL,
                       choices = c("LogNormalize" = "LogNormalize",
                                   "SCT" = "SCT"),
                       selected = "LogNormalize"),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'LogNormalize'", ns("method")),
      label_with_help("Scale factor",
                      "Common total each cell is scaled to before log. 10,000 is the standard default."),
      shiny::numericInput(ns("scale_factor"), NULL, value = 1e4, min = 1, step = 1e3)
    ),
    run_button(ns("run"), "Normalize")
  )
  step_container(title = list(en = "Normalization", zh = "归一化"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview")))
}

#' @rdname mod_normalize
#' @keywords internal
mod_normalize_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(done = FALSE, method = NULL, scale_factor = NULL,
                                 md = NULL)

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs("Seurat", "Normalization")) return(NULL)
      method <- input$method
      sf <- input$scale_factor
      obj <- with_progress_notify({
        o <- normalize_obj(rv$obj, method = method, scale_factor = sf)
        res$md <- obj_meta(o)
        o
      }, message = "Normalizing counts...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      res$done <- TRUE
      res$method <- method
      res$scale_factor <- sf
      mark_done(rv, "normalize")
      log_step(log_rv, "Normalization",
               params = list(method = method,
                             scale.factor = if (method == "LogNormalize") sf else NA),
               code = if (method == "LogNormalize") {
                 sprintf("obj <- Seurat::NormalizeData(obj, normalization.method = 'LogNormalize', scale.factor = %g)", sf)
               } else {
                 "obj <- Seurat::SCTransform(obj)"
               })
      shiny::showNotification(sprintf("Normalization done (%s).", method),
                              type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) return(shiny::div(class = "scstudio-placeholder",
                                               "Pick a method and click Normalize."))
      bslib::layout_columns(
        col_widths = c(6, 6),
        stat_tile("Method", res$method),
        stat_tile("Scale factor",
                  if (identical(res$method, "LogNormalize"))
                    format(res$scale_factor, big.mark = ",") else "n/a")
      )
    })

    output$preview <- render_preview_plot(function() {
      md <- res$md
      shiny::req(md)
      shiny::req(!is.null(md$nCount_RNA))
      # Robust view of sequencing depth: distribution of per-cell library sizes
      df <- data.frame(cell = rownames(md),
                       nCount = md$nCount_RNA,
                       stringsAsFactors = FALSE)
      df$log_count <- log10(df$nCount + 1)
      df$text <- sprintf("cell=%s\nUMIs=%s", df$cell,
                         format(df$nCount, big.mark = ","))
      ggplot2::ggplot(df, ggplot2::aes(x = log_count, text = text)) +
        ggplot2::geom_histogram(bins = 50, fill = scstudio_palette(1), alpha = 0.85) +
        ggplot2::labs(x = "Library size (log10 UMIs)", y = "Cells",
                      title = "Per-cell sequencing depth") +
        scstudio_theme()
    })
  })
}
