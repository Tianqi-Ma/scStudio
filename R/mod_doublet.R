#' Module 3: Doublet removal (去双细胞)
#'
#' Detect droplets that captured two cells (doublets) with scDblFinder (default)
#' or DoubletFinder, then either flag them or drop them from the object.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_doublet
NULL

#' @rdname mod_doublet
#' @keywords internal
mod_doublet_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = "Doublet removal (去双细胞)",
    what = "Detect and remove <b>doublets</b>: droplets that accidentally captured
            two cells instead of one.",
    why  = "A doublet = one droplet captured two cells, masquerading as a fake
            'intermediate' cell type. Left in, they create spurious clusters and
            confuse downstream annotation; remove them.",
    how  = "The recommended <b>scDblFinder</b> method scores every cell for how
            doublet-like it is. Choose whether to just flag doublets (keep all
            cells but label them) or remove them (drop the flagged cells).",
    example = "Two cells of different types share a droplet and look like a novel
               'hybrid' population &mdash; scDblFinder flags them so you can drop them."
  )
  controls <- shiny::tagList(
    label_with_help("Detection method",
                    "scDblFinder = fast one-click default (recommended). DoubletFinder = classic method needing a per-dataset pK sweep."),
    shiny::selectInput(ns("method"), NULL,
                       choices = c("scDblFinder" = "scDblFinder",
                                   "DoubletFinder" = "DoubletFinder"),
                       selected = "scDblFinder"),
    label_with_help("Action",
                    "Flag = keep every cell but label it doublet/singlet. Remove = drop cells classed as doublets."),
    shiny::radioButtons(ns("action"), NULL,
                        c("Flag only" = "flag", "Remove doublets" = "remove"),
                        selected = "remove", inline = TRUE),
    label_with_help("Custom score threshold (optional)",
                    "Leave blank to use the method's own call. Set a value to class cells with doublet_score above it as doublets."),
    shiny::numericInput(ns("threshold"), NULL, value = NA, min = 0, max = 1, step = 0.05),
    run_button(ns("run"), "Detect doublets")
  )
  step_container(title = list(en = "Doublet removal", zh = "去除双细胞"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview")))
}

#' @rdname mod_doublet
#' @keywords internal
mod_doublet_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(before = NA, after = NA, n_doublet = NA,
                                 md = NULL, action = NULL)

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs(c("Seurat", "scDblFinder"), "Doublet removal")) return(NULL)
      thr <- if (is.null(input$threshold) || is.na(input$threshold)) NULL else input$threshold
      action <- input$action
      obj <- with_progress_notify({
        o <- run_doublets(rv$obj, method = input$method)
        md <- obj_meta(o)
        # Optional custom score threshold overrides the method's own call
        if (!is.null(thr) && !is.null(md$doublet_score)) {
          cls <- ifelse(md$doublet_score > thr, "doublet", "singlet")
          o[["doublet_class"]] <- cls
        }
        res$before <- ncol(o)
        res$md <- obj_meta(o)
        is_doub <- res$md$doublet_class == "doublet"
        res$n_doublet <- sum(is_doub, na.rm = TRUE)
        if (action == "remove") o[, !is_doub] else o
      }, message = "Detecting doublets...")
      if (is.null(obj)) return(NULL)
      res$after <- ncol(obj)
      res$action <- action
      rv$obj <- obj
      mark_done(rv, "doublet")
      log_step(log_rv, "Doublet removal",
               params = list(method = input$method, action = action,
                             threshold = thr),
               code = if (action == "remove") {
                 "obj <- run_doublets(obj); obj <- obj[, obj$doublet_class == 'singlet']"
               } else {
                 "obj <- run_doublets(obj)  # flag only"
               })
      shiny::showNotification(
        sprintf("Doublet detection done: %d doublets (%s).",
                res$n_doublet, if (action == "remove") "removed" else "flagged"),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (is.na(res$before)) return(shiny::div(class = "scstudio-placeholder",
                                               "Pick a method and click Detect doublets."))
      pct <- if (res$before > 0) 100 * res$n_doublet / res$before else 0
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile("Cells in", format(res$before, big.mark = ",")),
        stat_tile("Doublets", sprintf("%s (%.1f%%)",
                                      format(res$n_doublet, big.mark = ","), pct)),
        stat_tile(if (identical(res$action, "remove")) "Kept" else "Cells out",
                  format(res$after, big.mark = ","))
      )
    })

    output$preview <- render_preview_plot(function() {
      md <- res$md
      shiny::req(md)
      shiny::req(!is.null(md$doublet_score))
      md$cellid <- rownames(md)
      md$doublet_class <- ifelse(is.na(md$doublet_class), "singlet", md$doublet_class)
      md$text <- sprintf("score=%.3f\nclass=%s", md$doublet_score, md$doublet_class)
      ggplot2::ggplot(md, ggplot2::aes(x = doublet_score, fill = doublet_class,
                                       text = text)) +
        ggplot2::geom_histogram(bins = 50, alpha = 0.75, position = "identity") +
        ggplot2::scale_fill_manual(values = c(singlet = "#3b6ea5", doublet = "#c1476b"),
                                   name = NULL) +
        ggplot2::labs(x = "Doublet score", y = "Cells",
                      title = "Doublet score distribution") +
        scstudio_theme()
    })
  })
}
