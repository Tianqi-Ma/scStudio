#' Module 2: Quality control
#'
#' Compute per-cell QC metrics and filter low-quality cells using either
#' adaptive MAD-based thresholds (recommended) or manual cutoffs.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_qc
NULL

#' @rdname mod_qc
#' @keywords internal
mod_qc_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = "Quality control",
    what = "Remove low-quality cells: empty droplets, dying cells, and debris.",
    why  = "Dying cells leak cytoplasmic RNA and show high mitochondrial content;
            empty droplets have very few genes. Keeping them adds noise.",
    how  = "The recommended <b>MAD</b> method flags cells that are statistical
            outliers for their own dataset (no guessing fixed numbers). Increase
            the MAD multiplier to keep more cells; decrease to be stricter.",
    example = "A cell with 40% mitochondrial reads is likely dying and gets flagged;
               a healthy cell (~5%) is kept."
  )
  controls <- shiny::tagList(
    label_with_help("Species", "Sets gene-name patterns for mitochondrial/ribosomal/hemoglobin genes."),
    shiny::selectInput(ns("species"), NULL, c("Human" = "human", "Mouse" = "mouse")),
    label_with_help("Threshold method",
                    "MAD = adaptive, data-driven (recommended). Manual = you set fixed cutoffs."),
    shiny::radioButtons(ns("method"), NULL,
                        c("MAD (adaptive)" = "mad", "Manual" = "manual"), inline = TRUE),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'mad'", ns("method")),
      label_with_help("MAD multiplier (library size / genes)",
                      "Higher = more permissive. 5 is a common default."),
      shiny::sliderInput(ns("nmad_lib"), NULL, min = 2, max = 8, value = 5, step = 0.5),
      label_with_help("MAD multiplier (mito %)", "Upper-tail only. 3 is common."),
      shiny::sliderInput(ns("nmad_mt"), NULL, min = 2, max = 8, value = 3, step = 0.5)
    ),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'manual'", ns("method")),
      shiny::numericInput(ns("min_genes"), "Min genes/cell", 200, min = 0),
      shiny::numericInput(ns("max_genes"), "Max genes/cell", 6000, min = 0),
      shiny::numericInput(ns("max_mt"), "Max mito %", 15, min = 0, max = 100)
    ),
    run_button(ns("run"), "Compute & filter")
  )
  step_container(title = list(en = "Quality control", zh = "质量控制"),
                 explainer = explainer, controls = controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview")))
}

#' @rdname mod_qc
#' @keywords internal
mod_qc_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(before = NA, after = NA, keep = NULL, md = NULL)

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs("Seurat", "QC")) return(NULL)
      obj <- with_progress_notify({
        o <- qc_add_metrics(rv$obj, species = input$species)
        keep <- if (input$method == "mad") {
          qc_mad_keep(o, input$nmad_lib, input$nmad_mt)
        } else {
          qc_manual_keep(o, input$min_genes, input$max_genes, input$max_mt)
        }
        res$before <- ncol(o)
        res$md <- obj_meta(o)
        res$md$keep <- keep
        res$keep <- keep
        o[, keep]
      }, message = "Computing QC and filtering...")
      if (is.null(obj)) return(NULL)
      res$after <- ncol(obj)
      rv$obj <- obj
      mark_done(rv, "qc")
      log_step(log_rv, "QC",
               params = list(method = input$method, species = input$species,
                             nmad_lib = input$nmad_lib, nmad_mt = input$nmad_mt),
               code = "obj <- subset(obj, cells = keep_cells)  # MAD/manual QC")
      shiny::showNotification(sprintf("QC done: kept %d of %d cells.",
                                      res$after, res$before), type = "message")
    })

    output$summary <- shiny::renderUI({
      if (is.na(res$before)) return(shiny::div(class = "scstudio-placeholder",
                                               "Set thresholds and click Compute & filter."))
      removed <- res$before - res$after
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile("Before", format(res$before, big.mark = ",")),
        stat_tile("Kept", format(res$after, big.mark = ",")),
        stat_tile("Removed", format(removed, big.mark = ","))
      )
    })

    output$preview <- render_preview_plot(function() {
      md <- res$md
      shiny::req(md)
      md$cellid <- rownames(md)
      md$text <- sprintf("nCount=%s\nnGenes=%s\nmito=%.1f%%\n%s",
                         format(md$nCount_RNA, big.mark = ","),
                         format(md$nFeature_RNA, big.mark = ","),
                         ifelse(is.null(md$percent.mt), NA, md$percent.mt),
                         ifelse(md$keep, "kept", "flagged"))
      ggplot2::ggplot(md, ggplot2::aes(x = nCount_RNA, y = percent.mt,
                                       colour = keep, text = text)) +
        ggplot2::geom_point(size = 0.5, alpha = 0.6) +
        ggplot2::scale_x_log10() +
        ggplot2::scale_colour_manual(values = c(`TRUE` = "#3b6ea5", `FALSE` = "#c1476b"),
                                     labels = c(`TRUE` = "kept", `FALSE` = "flagged"),
                                     name = NULL) +
        ggplot2::labs(x = "UMI count (log10)", y = "Mitochondrial %",
                      title = "QC: cells kept vs flagged") +
        scstudio_theme()
    })
  })
}
