#' Module: Export & reproducibility
#'
#' Save the processed object for later use and download a runnable R script that
#' reproduces every step you performed, built from the reproducibility log.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_export
NULL

#' @rdname mod_export
#' @keywords internal
mod_export_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = "Export & reproducibility",
    what = "Download your processed object and a script that reproduces the whole
            analysis.",
    why  = "Reproducibility is the point: anyone (including future you) should be
            able to regenerate these results from the raw data and the script.",
    how  = "Choose <b>RDS</b> to reload the object in R/Seurat, or <b>.h5ad</b> for
            Python/Scanpy (best-effort, needs SeuratDisk). The R script lists every
            step in order -- you only need to set the input path where it says so.",
    example = "Re-run with <code>source(\"scstudio_analysis.R\")</code> after
               editing the <code>input_path</code> line at the top."
  )
  controls <- shiny::tagList(
    label_with_help("Object format",
                    "RDS = native R/Seurat. .h5ad = AnnData for Python/Scanpy (needs SeuratDisk). Figures = guidance only."),
    shiny::selectInput(ns("fmt"), NULL,
                       choices = c("RDS (.rds)"   = "rds",
                                   "AnnData (.h5ad)" = "h5ad",
                                   "Figures (note)"  = "figures"),
                       selected = "rds"),
    shiny::downloadButton(ns("download_obj"), "Download object", class = "w-100"),
    shiny::tags$hr(),
    label_with_help("Reproducibility script",
                    "A commented R script rebuilding every step you ran, in order."),
    shiny::downloadButton(ns("download_script"), "Download R script", class = "w-100")
  )
  step_container(
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::uiOutput(ns("preview"))
  )
}

#' @rdname mod_export
#' @keywords internal
mod_export_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {

    # Build the reproducibility R script text from the log entries.
    build_script <- function(entries) {
      header <- c(
        "# scStudio reproducibility script",
        paste0("# Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        "#",
        "# IMPORTANT: set the input path below to your own data before running.",
        'input_path <- "PATH/TO/YOUR/DATA"  # <-- edit this',
        "library(Seurat)",
        ""
      )
      if (is.null(entries) || !length(entries)) {
        return(paste(c(header, "# No steps were recorded."), collapse = "\n"))
      }
      body <- unlist(lapply(seq_along(entries), function(i) {
        e <- entries[[i]]
        params <- if (length(e$params)) {
          paste(vapply(names(e$params), function(k)
            sprintf("%s=%s", k, paste(deparse(e$params[[k]]), collapse = "")),
            character(1)), collapse = ", ")
        } else ""
        c(sprintf("# Step %d: %s", i, e$step),
          if (nzchar(params)) sprintf("#   params: %s", params) else NULL,
          if (length(e$code)) e$code else "# (no code recorded)",
          "")
      }))
      paste(c(header, body), collapse = "\n")
    }

    output$summary <- shiny::renderUI({
      entries <- log_rv()
      if (is.null(entries) || !length(entries)) {
        return(shiny::div(class = "scstudio-placeholder",
                          "No steps recorded yet. Run some analysis steps first."))
      }
      shiny::tagList(
        stat_tile("Steps performed", length(entries)),
        shiny::tags$ol(class = "scstudio-steps",
          lapply(entries, function(e) {
            shiny::tags$li(shiny::tags$b(e$step),
                           shiny::tags$span(class = "scstudio-muted",
                                            paste0("  (", e$time, ")")))
          }))
      )
    })

    output$preview <- shiny::renderUI({
      shiny::div(class = "scstudio-note",
                 "No plot for this step. Use the buttons on the left to download ",
                 "your processed object and the reproducibility script. ",
                 "Remember to set input/output paths when you re-run the script.")
    })

    output$download_obj <- shiny::downloadHandler(
      filename = function() {
        ext <- switch(input$fmt, rds = "rds", h5ad = "h5ad", figures = "txt")
        paste0("scstudio_object_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
      },
      content = function(file) {
        obj <- rv$obj
        if (is.null(obj)) {
          writeLines("No object available to export.", file)
          return(invisible(NULL))
        }
        if (input$fmt == "rds") {
          saveRDS(obj, file)
        } else if (input$fmt == "h5ad") {
          if (!require_pkgs(c("SeuratDisk", "Seurat"), "AnnData (.h5ad) export")) {
            writeLines("SeuratDisk not installed; could not write .h5ad.", file)
            return(invisible(NULL))
          }
          tmp <- tempfile(fileext = ".h5Seurat")
          SeuratDisk::SaveH5Seurat(obj, filename = tmp, overwrite = TRUE)
          SeuratDisk::Convert(tmp, dest = "h5ad", overwrite = TRUE)
          file.copy(sub("\\.h5Seurat$", ".h5ad", tmp), file, overwrite = TRUE)
        } else {
          writeLines(c(
            "Figures are downloaded individually from the Visualize step",
            "using its 'Download plot' button (PNG/PDF)."), file)
        }
      }
    )

    output$download_script <- shiny::downloadHandler(
      filename = function()
        paste0("scstudio_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".R"),
      content = function(file) {
        writeLines(build_script(log_rv()), file)
      }
    )
  })
}
