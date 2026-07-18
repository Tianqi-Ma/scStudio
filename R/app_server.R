#' Top-level server: shared state, step-status tracking, navigation
#'
#' A single reactive hub (`rv`) holds the working object, a per-step status map
#' (used to colour the left navigator), and the reproducibility log. Each module
#' reads `rv$obj`, does its work, writes the object back, and marks its step done.
#'
#' @param input,output,session Standard Shiny server arguments.
#' @keywords internal
app_server <- function(input, output, session) {

  rv <- shiny::reactiveValues(
    obj    = NULL,
    source = NULL,
    status = list()   # step key -> TRUE once it has run successfully
  )
  log_rv <- shiny::reactiveVal(list())

  # Ordered steps (must match app_ui) for rendering the navigator.
  steps <- list(
    list(v = "import",    n = 1,  en = "Import",       zh = "导入"),
    list(v = "qc",        n = 2,  en = "QC",           zh = "质控"),
    list(v = "doublet",   n = 3,  en = "Doublets",     zh = "去双细胞"),
    list(v = "normalize", n = 4,  en = "Normalize",    zh = "归一化"),
    list(v = "reduce",    n = 5,  en = "Features/PCA", zh = "特征/PCA"),
    list(v = "integrate", n = 6,  en = "Integrate",    zh = "整合"),
    list(v = "cluster",   n = 7,  en = "Cluster",      zh = "聚类"),
    list(v = "embed",     n = 8,  en = "Embed",        zh = "降维图"),
    list(v = "markers",   n = 9,  en = "Markers",      zh = "标志基因"),
    list(v = "annotate",  n = 10, en = "Annotate",     zh = "注释"),
    list(v = "viz",       n = 11, en = "Visualize",    zh = "可视化"),
    list(v = "export",    n = 12, en = "Export",       zh = "导出")
  )

  # --- Left step navigator (status-coloured, clickable) ----------------------
  output$step_nav <- shiny::renderUI({
    current <- input$steps %||% "import"
    status  <- rv$status
    items <- lapply(steps, function(s) {
      state <- if (identical(s$v, current)) "current"
               else if (isTRUE(status[[s$v]])) "done" else "todo"
      shiny::tags$a(
        class = paste("scstudio-navitem", state),
        onclick = sprintf("Shiny.setInputValue('goto','%s',{priority:'event'})", s$v),
        shiny::span(class = "scstudio-navdot"),
        shiny::span(class = "scstudio-navnum", s$n),
        i18n(s$en, s$zh)
      )
    })
    shiny::div(class = "scstudio-nav", items)
  })

  shiny::observeEvent(input$goto, {
    bslib::nav_select("steps", input$goto)
  })

  # --- Global dataset status (bottom of the sidebar) -------------------------
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

  # --- Wire modules; each returns nothing but updates rv (incl. rv$status) ----
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

#' NULL-coalescing helper
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Mark a pipeline step as completed (colours the navigator)
#' @param rv Shared hub. @param step Step key.
#' @keywords internal
mark_done <- function(rv, step) {
  st <- rv$status
  st[[step]] <- TRUE
  rv$status <- st
  invisible(TRUE)
}

#' Small labelled status line for the sidebar
#' @param label,value Character/UI.
#' @keywords internal
stat_line <- function(label, value) {
  shiny::div(class = "scstudio-statline",
             shiny::span(class = "scstudio-statlabel", label),
             shiny::span(class = "scstudio-statvalue", value))
}
