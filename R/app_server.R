#' Top-level server: shared state + module wiring
#'
#' A single reactive "hub" object (`rv`) holds the working single-cell object and
#' flows through every module. Each analysis module reads `rv$obj`, does its work,
#' and writes the updated object back -- mirroring the sequential pipeline while
#' still allowing users to revisit earlier steps.
#'
#' @param input,output,session Standard Shiny server arguments.
#' @keywords internal
app_server <- function(input, output, session) {

  # --- Shared reactive state -------------------------------------------------
  rv <- shiny::reactiveValues(
    obj    = NULL,   # the working Seurat object
    source = NULL,   # description of where it came from
    log    = NULL    # reproducibility log (set below as reactiveVal)
  )
  log_rv <- shiny::reactiveVal(list())

  # --- Global dataset status card (in the left sidebar) ----------------------
  output$global_status <- shiny::renderUI({
    obj <- rv$obj
    if (is.null(obj)) {
      return(shiny::div(class = "scstudio-status-empty",
                        shiny::icon("circle-dot"), " No data loaded.",
                        shiny::tags$p(class = "text-muted small",
                                      "Start at step 1 (Import).")))
    }
    dims <- obj_dims(obj)
    advice <- memory_advice(dims$cells)
    shiny::tagList(
      stat_line("Cells", format(dims$cells, big.mark = ",")),
      stat_line("Genes", format(dims$genes, big.mark = ",")),
      if (!is.null(rv$source)) stat_line("Source", rv$source),
      if (nzchar(advice))
        shiny::div(class = "scstudio-warn", shiny::icon("triangle-exclamation"), " ", advice)
    )
  })

  # --- Wire up every module (all receive and can update the hub) -------------
  mod_import_server("import", rv, log_rv, parent = session)
  mod_qc_server("qc", rv, log_rv)
  mod_doublet_server("doublet", rv, log_rv)
  mod_normalize_server("normalize", rv, log_rv)
  mod_reduce_server("reduce", rv, log_rv)
  mod_integrate_server("integrate", rv, log_rv)
  mod_cluster_server("cluster", rv, log_rv)
  mod_embed_server("embed", rv, log_rv)
  mod_markers_server("markers", rv, log_rv)
  mod_annotate_server("annotate", rv, log_rv)
  mod_viz_server("viz", rv, log_rv)
  mod_export_server("export", rv, log_rv)
}

#' Small labelled status line for the sidebar
#' @param label,value Character.
#' @keywords internal
stat_line <- function(label, value) {
  shiny::div(class = "scstudio-statline",
             shiny::tags$span(class = "scstudio-statlabel", label),
             shiny::tags$span(class = "scstudio-statvalue", value))
}
