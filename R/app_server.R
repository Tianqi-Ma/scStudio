#' Top-level server: shared state, grouped stepper, navigation
#'
#' A single reactive hub (`rv`) holds the working object, a per-step status map
#' (colours the left navigator), and the reproducibility log. Each module reads
#' `rv$obj`, does its work, writes it back, and marks its step done. Language is
#' handled entirely client-side (app.js); theme via bslib input_dark_mode.
#'
#' @param input,output,session Standard Shiny server arguments.
#' @keywords internal
app_server <- function(input, output, session) {

  rv <- shiny::reactiveValues(obj = NULL, source = NULL, status = list())
  log_rv <- shiny::reactiveVal(list())

  # --- Shared "export current data" (available on every step) ----------------
  register_exports(input, output, session, rv)

  # --- Grouped, status-coloured left navigator -------------------------------
  output$step_nav <- shiny::renderUI({
    current <- input$steps %||% "import"
    status  <- rv$status
    phases  <- app_phases()
    steps   <- app_steps()
    # order phases by first appearance in steps
    phase_order <- unique(vapply(steps, function(s) s$phase, character(1)))
    children <- list()
    for (ph in phase_order) {
      lab <- phases[[ph]]
      children[[length(children) + 1]] <-
        shiny::div(class = "scstudio-phase", i18n(lab$en, lab$zh))
      for (s in Filter(function(x) identical(x$phase, ph), steps)) {
        state <- if (identical(s$v, current)) "current"
                 else if (isTRUE(status[[s$v]])) "done" else "todo"
        children[[length(children) + 1]] <- shiny::tags$a(
          class = paste("scstudio-navitem", state),
          onclick = sprintf("Shiny.setInputValue('goto','%s',{priority:'event'})", s$v),
          shiny::span(class = "scstudio-navdot"),
          shiny::span(class = "scstudio-navnum", s$n),
          shiny::span(class = "scstudio-navlabel", i18n(s$en, s$zh))
        )
      }
    }
    shiny::div(class = "scstudio-nav", children)
  })

  shiny::observeEvent(input$goto, { bslib::nav_select("steps", input$goto) })

  # --- Global dataset status (bottom of sidebar) -----------------------------
  output$global_status <- shiny::renderUI({
    obj <- rv$obj
    if (is.null(obj)) {
      return(shiny::div(class = "scstudio-status-empty",
                        i18n("No data loaded.", "尚未加载数据")))
    }
    dims <- obj_dims(obj)
    advice <- memory_advice(dims$cells)
    shiny::tagList(
      stat_line(i18n("Cells", "细胞"), format(dims$cells, big.mark = ",")),
      stat_line(i18n("Genes", "基因"), format(dims$genes, big.mark = ",")),
      if (nzchar(advice))
        shiny::div(class = "scstudio-warn", shiny::icon("triangle-exclamation"), " ", advice)
    )
  })

  # --- Wire modules ----------------------------------------------------------
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
  mod_enrichment_server("enrichment", rv, log_rv)
  mod_trajectory_server("trajectory", rv, log_rv)
  mod_velocity_server("velocity", rv, log_rv)
  mod_dynamic_server("dynamic", rv, log_rv)
  mod_cellcycle_signatures_server("cellcycle", rv, log_rv)
  mod_cellcomm_server("cellcomm", rv, log_rv)
  mod_malignancy_server("malignancy", rv, log_rv)
  mod_viz_server("viz", rv, log_rv)
  mod_report_server("report", rv, log_rv)
  mod_export_server("export", rv, log_rv)
}

#' NULL-coalescing helper
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Mark a pipeline step completed (colours the navigator)
#' @param rv Shared hub. @param step Step key.
#' @keywords internal
mark_done <- function(rv, step) {
  st <- rv$status; st[[step]] <- TRUE; rv$status <- st; invisible(TRUE)
}

#' Small labelled status line for the sidebar
#' @keywords internal
stat_line <- function(label, value) {
  shiny::div(class = "scstudio-statline",
             shiny::span(class = "scstudio-statlabel", label),
             shiny::span(class = "scstudio-statvalue", value))
}
