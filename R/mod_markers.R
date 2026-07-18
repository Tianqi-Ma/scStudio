#' Module: Marker genes
#'
#' Find genes that are differentially expressed in each cluster relative to the
#' rest of the cells. These marker genes are what you use to give a cluster a
#' biological identity in the next (annotation) step.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_markers
NULL

#' @rdname mod_markers
#' @keywords internal
mod_markers_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = "Marker genes",
    what = "Find the genes that are specifically up- (or down-) regulated in each
            cluster compared with all other cells.",
    why  = "Clusters are just groups of similar cells until you know what makes
            them different. Marker genes are the evidence you use to call a
            cluster a cell type (e.g. <code>CD3D</code> for T cells).",
    how  = "<b>Wilcoxon</b> is the fast, robust default. Raise the log fold-change
            or min.pct to keep only stronger, more specific markers. Keep
            <b>only positive</b> markers if you only care about what a cluster
            expresses <i>more</i> than others.",
    example = "For a T-cell cluster you would expect markers like <code>CD3D</code>,
               <code>CD3E</code> and <code>TRAC</code> at the top of the list."
  )
  controls <- shiny::tagList(
    label_with_help("Statistical test",
                    "Wilcoxon = fast rank test (default). ROC = ranks genes by classification power. MAST = models dropout (needs the MAST package)."),
    shiny::selectInput(ns("test"), NULL,
                       choices = c("Wilcoxon" = "wilcox", "ROC" = "roc", "MAST" = "MAST"),
                       selected = "wilcox"),
    label_with_help("Log fold-change threshold",
                    "Minimum log2 fold-change to test a gene. Higher = fewer, stronger markers."),
    shiny::numericInput(ns("logfc"), NULL, value = 0.25, min = 0, step = 0.05),
    label_with_help("Min fraction expressing (min.pct)",
                    "A gene must be detected in at least this fraction of cells in one of the two groups."),
    shiny::numericInput(ns("min_pct"), NULL, value = 0.1, min = 0, max = 1, step = 0.05),
    shiny::checkboxInput(ns("only_pos"), "Only positive markers", value = TRUE),
    label_with_help("Top N per cluster",
                    "How many top markers per cluster to show in the preview chart."),
    shiny::numericInput(ns("top_n"), NULL, value = 10, min = 1, max = 50, step = 1),
    run_button(ns("run"), "Find markers")
  )
  step_container(
    title     = list(en = "Marker genes", zh = "标志基因"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::tagList(
      preview_plot_ui(ns("preview")),
      shiny::div(class = "scstudio-table", shiny::uiOutput(ns("table")))
    )
  )
}

#' @rdname mod_markers
#' @keywords internal
mod_markers_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    markers <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs("Seurat", "Marker genes")) return(NULL)
      df <- with_progress_notify({
        markers_obj(rv$obj, test = input$test, logfc = input$logfc,
                    min_pct = input$min_pct, only_pos = input$only_pos)
      }, message = "Finding marker genes...")
      if (is.null(df)) return(NULL)
      markers(df)
      rv$markers <- df
      mark_done(rv, "markers")
      log_step(log_rv, "Markers",
               params = list(test = input$test, logfc = input$logfc,
                             min_pct = input$min_pct, only_pos = input$only_pos),
               code = sprintf(
                 'markers <- Seurat::FindAllMarkers(obj, test.use="%s", logfc.threshold=%s, min.pct=%s, only.pos=%s)',
                 input$test, input$logfc, input$min_pct, input$only_pos))
      shiny::showNotification(sprintf("Found %d markers across %d clusters.",
                                      nrow(df), length(unique(df$cluster))),
                              type = "message")
    })

    output$summary <- shiny::renderUI({
      df <- markers()
      if (is.null(df)) {
        return(shiny::div(class = "scstudio-placeholder",
                          "Set options and click ", shiny::tags$b("Find markers"), "."))
      }
      n_clusters <- length(unique(df$cluster))
      genes <- if ("gene" %in% colnames(df)) length(unique(df$gene)) else NA_integer_
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile("Markers", format(nrow(df), big.mark = ",")),
        stat_tile("Clusters", format(n_clusters, big.mark = ",")),
        stat_tile("Unique genes", format(genes, big.mark = ","))
      )
    })

    # Top-N markers per cluster as a tidy data.frame (shared by chart + table)
    top_markers <- shiny::reactive({
      df <- markers()
      shiny::req(df)
      n <- max(1, as.integer(input$top_n))
      parts <- split(df, df$cluster)
      picked <- lapply(parts, function(d) {
        ord <- if ("avg_log2FC" %in% colnames(d)) order(-d$avg_log2FC) else seq_len(nrow(d))
        utils::head(d[ord, , drop = FALSE], n)
      })
      do.call(rbind, picked)
    })

    output$preview <- render_preview_plot(function() {
      d <- top_markers()
      shiny::req(d)
      d$cluster <- factor(d$cluster)
      lfc <- if ("avg_log2FC" %in% colnames(d)) d$avg_log2FC else rep(NA_real_, nrow(d))
      padj <- if ("p_val_adj" %in% colnames(d)) d$p_val_adj else rep(NA_real_, nrow(d))
      gene <- if ("gene" %in% colnames(d)) d$gene else rownames(d)
      d$gene_lab <- gene
      d$avg_log2FC <- lfc
      d$p_val_adj  <- padj
      d$text <- sprintf("gene: %s\ncluster: %s\navg_log2FC: %.2f\np_val_adj: %.2g",
                        gene, as.character(d$cluster), lfc, padj)
      ggplot2::ggplot(d, ggplot2::aes(x = stats::reorder(gene_lab, avg_log2FC),
                                      y = avg_log2FC, fill = cluster, text = text)) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::facet_wrap(~cluster, scales = "free_y") +
        ggplot2::scale_fill_manual(values = scstudio_palette(nlevels(d$cluster)),
                                   guide = "none") +
        ggplot2::labs(x = NULL, y = "avg log2 fold-change",
                      title = "Top marker genes per cluster") +
        scstudio_theme()
    })

    output$table <- shiny::renderUI({
      shiny::req(markers())
      ns <- session$ns
      if (has_pkg("DT")) {
        DT::dataTableOutput(ns("dt"))
      } else {
        shiny::verbatimTextOutput(ns("txt"))
      }
    })

    if (has_pkg("DT")) {
      output$dt <- DT::renderDataTable({
        df <- markers()
        shiny::req(df)
        DT::datatable(df, filter = "top", rownames = FALSE,
                      options = list(pageLength = 10, scrollX = TRUE))
      })
    } else {
      output$txt <- shiny::renderPrint({
        df <- markers()
        shiny::req(df)
        utils::head(df, 20)
      })
    }
  })
}
