#' Module 1: Import & inspect
#'
#' Upload a count matrix and turn it into the working object. Accepts a Seurat or
#' SingleCellExperiment `.rds`, a 10x directory/`.h5`, or a plain counts table
#' (csv/tsv, genes x cells). Old Seurat objects are updated automatically.
#'
#' @param id Module id.
#' @name mod_import
NULL

#' @rdname mod_import
#' @keywords internal
mod_import_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = "Import your data",
    what = "Load your single-cell count data into the app as the working object.",
    why  = "Every later step operates on this object. Counts are the raw number of
            transcripts detected per gene per cell.",
    how  = "<b>Just exploring?</b> Choose <b>Demo data</b> and click load to try the
            whole pipeline in seconds. To use your own data, pick <b>Upload file</b>
            and the format that matches it (RDS if it's a saved Seurat object).",
    example = "The bundled demo loads instantly with no download. Or upload
               <code>pbmc.rds</code> (a Seurat object), a 10x <code>.h5</code>, or a
               counts table (genes in rows, cells in columns)."
  )
  demos <- demo_datasets()
  demo_choices <- stats::setNames(demos$id, paste0(demos$name, "  (", demos$cells, " cells)"))

  controls <- shiny::tagList(
    label_with_help("Data source",
                    "New here? Pick 'Demo data' to try the app instantly. Otherwise upload your own file, or fetch one from a URL."),
    shiny::radioButtons(ns("source"), NULL,
                        c("Demo data" = "demo",
                          "Upload file" = "upload",
                          "From URL" = "url"),
                        selected = "demo"),

    # --- Demo data ---
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'demo'", ns("source")),
      label_with_help("Demo dataset",
                      "The bundled example loads instantly offline. The 10x PBMC sets are real data and download the first time (needs internet)."),
      shiny::selectInput(ns("demo_id"), NULL, choices = demo_choices, selected = "bundled"),
      shiny::helpText(shiny::textOutput(ns("demo_desc"), inline = TRUE)),
      run_button(ns("load_demo"), "Load demo data")
    ),

    # --- Upload ---
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'upload'", ns("source")),
      label_with_help("Input format",
                      "RDS = a saved Seurat/SingleCellExperiment object. 10x = Cell Ranger output. Table = a CSV/TSV of counts (genes in rows, cells in columns)."),
      shiny::selectInput(ns("fmt"), NULL,
                         choices = c("RDS (Seurat/SCE)" = "rds",
                                     "10x HDF5 (.h5)"    = "h5",
                                     "Counts table (csv/tsv)" = "table"),
                         selected = "rds"),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'table'", ns("fmt")),
        label_with_help("Separator", "How columns are separated in your table."),
        shiny::selectInput(ns("sep"), NULL,
                           choices = c("Tab" = "\t", "Comma" = ","), selected = "\t")
      ),
      shiny::fileInput(ns("file"), "Choose file",
                       accept = c(".rds", ".h5", ".csv", ".tsv", ".txt", ".gz")),
      run_button(ns("load"), "Load data")
    ),

    # --- From URL ---
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'url'", ns("source")),
      label_with_help("File URL",
                      "Direct link to a .rds or 10x .h5 file. It is downloaded to a temporary file on your machine."),
      shiny::textInput(ns("url"), NULL, placeholder = "https://.../data.h5"),
      shiny::selectInput(ns("url_fmt"), "Format",
                         choices = c("RDS (Seurat/SCE)" = "rds", "10x HDF5 (.h5)" = "h5"),
                         selected = "h5"),
      run_button(ns("load_url"), "Fetch & load")
    )
  )
  step_container(
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = preview_plot_ui(ns("preview"))
  )
}

#' @rdname mod_import
#' @param rv Shared reactiveValues hub (with `$obj`, `$source`).
#' @param log_rv reactiveVal reproducibility log.
#' @param parent Parent session (unused; reserved for nav control).
#' @keywords internal
mod_import_server <- function(id, rv, log_rv, parent = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    # Show the selected demo's description under the dropdown.
    output$demo_desc <- shiny::renderText({
      ds <- demo_datasets()
      row <- ds[ds$id == input$demo_id, , drop = FALSE]
      if (nrow(row)) row$description else ""
    })

    # Shared loader: read a file of a given format -> Seurat object -> hub.
    load_into_hub <- function(path, fmt, sep = "\t", source_label = fmt,
                              log_params = list(), log_file = basename(path)) {
      if (!require_pkgs(c("Seurat", "SeuratObject"), "Import")) return(invisible(NULL))
      obj <- with_progress_notify({
        loaded <- switch(
          fmt,
          rds   = readRDS(path),
          h5    = Seurat::Read10X_h5(path),
          table = read_counts_table(path, sep = sep)
        )
        as_seurat(loaded)
      }, message = "Loading and building object...")
      if (is.null(obj)) return(invisible(NULL))
      rv$obj    <- obj
      rv$source <- source_label
      log_step(log_rv, "Import",
               params = log_params,
               code = sprintf('obj <- %s',
                              switch(fmt,
                                     rds   = sprintf('readRDS("%s")', log_file),
                                     h5    = sprintf('Seurat::Read10X_h5("%s")', log_file),
                                     table = sprintf('read.delim("%s", row.names=1)', log_file))))
      shiny::showNotification("Data loaded.", type = "message")
    }

    # (a) Upload
    shiny::observeEvent(input$load, {
      shiny::req(input$file)
      load_into_hub(input$file$datapath, input$fmt, sep = input$sep,
                    source_label = switch(input$fmt, rds = "RDS", h5 = "10x .h5", table = "Table"),
                    log_params = list(source = "upload", format = input$fmt, file = input$file$name),
                    log_file = input$file$name)
    })

    # (b) Demo data
    shiny::observeEvent(input$load_demo, {
      got <- tryCatch(fetch_demo(input$demo_id),
                      error = function(e) { shiny::showNotification(conditionMessage(e), type = "error", duration = 12); NULL })
      shiny::req(got)
      nm <- demo_datasets()
      label <- nm$name[nm$id == input$demo_id]
      load_into_hub(got$path, got$format,
                    source_label = paste0("Demo: ", label),
                    log_params = list(source = "demo", demo = input$demo_id),
                    log_file = paste0("demo_", input$demo_id))
    })

    # (c) From URL
    shiny::observeEvent(input$load_url, {
      shiny::req(nzchar(input$url))
      ext <- if (input$url_fmt == "h5") ".h5" else ".rds"
      dest <- tempfile(fileext = ext)
      ok <- with_progress_notify(
        tryCatch(utils::download.file(input$url, dest, mode = "wb", quiet = TRUE) == 0,
                 error = function(e) FALSE),
        message = "Downloading...")
      if (!isTRUE(ok) || !file.exists(dest) || file.size(dest) == 0) {
        shiny::showNotification("Download failed (check the URL and your connection).",
                                type = "error", duration = 12)
        return(NULL)
      }
      load_into_hub(dest, input$url_fmt,
                    source_label = "URL",
                    log_params = list(source = "url", url = input$url, format = input$url_fmt),
                    log_file = basename(input$url))
    })

    output$summary <- shiny::renderUI({
      obj <- rv$obj
      if (is.null(obj)) {
        return(shiny::div(class = "scstudio-placeholder",
                          "No data yet. Tip: pick ", shiny::tags$b("Demo data"),
                          " and click ", shiny::tags$b("Load demo data"), " to try it instantly."))
      }
      dims <- obj_dims(obj)
      meta <- obj_meta(obj)
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        stat_tile("Cells", format(dims$cells, big.mark = ",")),
        stat_tile("Genes", format(dims$genes, big.mark = ",")),
        stat_tile("Metadata columns", ncol(meta))
      )
    })

    output$preview <- render_preview_plot(function() {
      obj <- rv$obj
      shiny::req(obj)
      meta <- obj_meta(obj)
      batch <- guess_batch_col(meta)
      if (is.null(batch)) {
        df <- data.frame(sample = "all cells", n = obj_dims(obj)$cells)
      } else {
        tab <- as.data.frame(table(meta[[batch]]), stringsAsFactors = FALSE)
        names(tab) <- c("sample", "n")
        df <- tab
      }
      df$text <- paste0(df$sample, ": ", format(df$n, big.mark = ","), " cells")
      ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(sample, -n), y = n, text = text)) +
        ggplot2::geom_col(fill = scstudio_palette(1)) +
        ggplot2::labs(x = NULL, y = "Cells", title = "Cells per sample") +
        scstudio_theme() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
    })
  })
}
