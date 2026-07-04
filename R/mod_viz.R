#' Module: Visualize
#'
#' A free exploration panel. Pick a plot type, a metadata column to group by, and
#' (for expression plots) some genes, then inspect the result and download it.
#' This module never modifies the working object -- it is read-only.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_viz
NULL

#' @rdname mod_viz
#' @keywords internal
mod_viz_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = "Visualize",
    what = "Explore your data freely: colour the embedding by any metadata column,
            or show the expression of specific genes.",
    why  = "A picture is the fastest way to sanity-check clustering, annotation and
            marker genes -- and to build the figures for your report.",
    how  = "Pick a plot type. <b>UMAP by metadata</b> colours cells by a column.
            <b>Violin / Dot / Feature / Heatmap</b> show expression of the genes
            you type (comma separated). Nothing here changes your object.",
    example = "Type <code>CD3D, MS4A1, LYZ</code> and choose 'Feature plot' to see
               where T cells, B cells and monocytes sit on the UMAP."
  )
  controls <- shiny::tagList(
    label_with_help("Plot type",
                    "UMAP colours cells by metadata; the others show gene expression."),
    shiny::selectInput(ns("ptype"), NULL,
                       choices = c("UMAP by metadata" = "umap",
                                   "Violin plot"      = "violin",
                                   "Dot plot"         = "dotplot",
                                   "Feature plot"     = "feature",
                                   "Heatmap"          = "heatmap"),
                       selected = "umap"),
    label_with_help("Group by (metadata column)",
                    "Which metadata column to colour or split cells by (e.g. seurat_clusters, celltype)."),
    shiny::selectInput(ns("meta_col"), NULL, choices = NULL),
    shiny::conditionalPanel(
      sprintf("input['%s'] != 'umap'", ns("ptype")),
      label_with_help("Genes",
                      "Comma-separated gene names for expression plots (violin/dot/feature/heatmap)."),
      shiny::textInput(ns("genes"), NULL, placeholder = "e.g. CD3D, MS4A1, LYZ")
    ),
    label_with_help("Download format", "File type for the downloaded figure."),
    shiny::radioButtons(ns("fmt"), NULL, c("PNG" = "png", "PDF" = "pdf"), inline = TRUE),
    shiny::downloadButton(ns("download"), "Download plot", class = "w-100")
  )
  step_container(
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::uiOutput(ns("plot_slot"))
  )
}

#' @rdname mod_viz
#' @keywords internal
mod_viz_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Keep the metadata-column selector in sync with the current object.
    shiny::observe({
      obj <- rv$obj
      cols <- if (is.null(obj)) character(0) else obj_meta_cols(obj)
      sel <- if ("celltype" %in% cols) "celltype"
             else if ("seurat_clusters" %in% cols) "seurat_clusters"
             else if (length(cols)) cols[1] else NULL
      shiny::updateSelectInput(session, "meta_col", choices = cols, selected = sel)
    })

    parse_genes <- function(txt) {
      if (is.null(txt) || !nzchar(trimws(txt))) return(character(0))
      g <- trimws(strsplit(txt, ",", fixed = TRUE)[[1]])
      g[nzchar(g)]
    }

    # Is this plot type worth making interactive (ggplotly)?
    is_interactive <- shiny::reactive({
      has_pkg("plotly") && input$ptype %in% c("umap", "feature")
    })

    # Build the current plot as a ggplot object (or NULL on failure).
    current_plot <- shiny::reactive({
      obj <- rv$obj
      shiny::req(obj)
      if (!require_pkgs("Seurat", "Visualization")) return(NULL)
      genes <- parse_genes(input$genes)
      red <- if (has_reduction(obj, "umap")) "umap" else obj_reductions(obj)[1]
      tryCatch({
        switch(input$ptype,
          umap = {
            shiny::validate(shiny::need(length(red) && !is.na(red),
                                        "No UMAP/embedding found. Run an embedding first."))
            Seurat::DimPlot(obj, reduction = red, group.by = input$meta_col) +
              scstudio_theme()
          },
          feature = {
            shiny::validate(shiny::need(length(genes) > 0, "Enter at least one gene."))
            shiny::validate(shiny::need(length(red) && !is.na(red),
                                        "No UMAP/embedding found. Run an embedding first."))
            Seurat::FeaturePlot(obj, features = genes, reduction = red) &
              scstudio_theme()
          },
          violin = {
            shiny::validate(shiny::need(length(genes) > 0, "Enter at least one gene."))
            Seurat::VlnPlot(obj, features = genes, group.by = input$meta_col) &
              scstudio_theme()
          },
          dotplot = {
            shiny::validate(shiny::need(length(genes) > 0, "Enter at least one gene."))
            Seurat::DotPlot(obj, features = genes, group.by = input$meta_col) +
              scstudio_theme()
          },
          heatmap = {
            shiny::validate(shiny::need(length(genes) > 0, "Enter at least one gene."))
            Seurat::DoHeatmap(obj, features = genes, group.by = input$meta_col)
          })
      }, error = function(e) {
        shiny::showNotification(paste("Plot error:", conditionMessage(e)),
                                type = "error", duration = 10)
        NULL
      })
    })

    # Choose the correct output widget for the current plot type.
    output$plot_slot <- shiny::renderUI({
      if (is_interactive()) {
        plotly::plotlyOutput(ns("iplot"), height = "480px")
      } else {
        shiny::plotOutput(ns("splot"), height = "480px")
      }
    })

    output$iplot <- plotly::renderPlotly({
      gg <- current_plot()
      shiny::req(gg)
      plotly::ggplotly(gg, tooltip = "text") |>
        plotly::config(displayModeBar = FALSE)
    })

    output$splot <- shiny::renderPlot({
      gg <- current_plot()
      shiny::req(gg)
      gg
    })

    output$summary <- shiny::renderUI({
      obj <- rv$obj
      if (is.null(obj)) {
        return(shiny::div(class = "scstudio-placeholder",
                          "Load and process data first, then explore it here."))
      }
      dims <- obj_dims(obj)
      genes <- parse_genes(input$genes)
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile("Cells", format(dims$cells, big.mark = ",")),
        stat_tile("Plot", input$ptype),
        stat_tile("Genes requested", length(genes))
      )
    })

    output$download <- shiny::downloadHandler(
      filename = function() paste0("scstudio_", input$ptype, "_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".", input$fmt),
      content = function(file) {
        gg <- current_plot()
        shiny::req(gg)
        ggplot2::ggsave(file, plot = gg, width = 8, height = 6, dpi = 300,
                        device = input$fmt)
      }
    )
  })
}
