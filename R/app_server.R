#' Top-level server: shared state, grouped stepper, navigation
#'
#' A single reactive hub (`rv`) holds the working object, a per-step status map
#' (colours the left navigator), the clinical/survival cohort, and the
#' reproducibility log. Each module reads `rv$obj`, does its work, writes it
#' back, and marks its step done; the staleness / epoch / checkpoint rules are
#' in fct_state.R. Language is handled client-side (app.js); theme via bslib
#' input_dark_mode.
#'
#' @param input,output,session Standard Shiny server arguments.
#' @keywords internal
app_server <- function(input, output, session) {

  rv <- shiny::reactiveValues(
    obj = NULL, source = NULL, status = list(),
    clinical = NULL,  # cohort loaded by the Clinical & survival step
    epoch_sc = 0L,    # pipeline-state bookkeeping, see fct_state.R
    ckpt = new.env(parent = emptyenv())
  )
  log_rv <- shiny::reactiveVal(list())

  # --- Shared "export current data" (available on every step) ----------------
  register_exports(input, output, session, rv)

  # --- Grouped, status-coloured left navigator (collapsible phase groups) ----
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
      ph_steps <- Filter(function(x) identical(x$phase, ph), steps)
      items <- lapply(ph_steps, function(s) {
        st <- step_state(rv, s$v)
        state <- if (identical(s$v, current)) "current" else st
        shiny::tags$a(
          class = paste("scstudio-navitem", state, if (st == "stale") "is-stale"),
          title = if (st == "stale") "Needs a re-run: an earlier step changed",
          onclick = sprintf("Shiny.setInputValue('goto','%s',{priority:'event'})", s$v),
          shiny::span(class = "scstudio-navdot"),
          shiny::span(class = "scstudio-navnum", s$n),
          shiny::span(class = "scstudio-navlabel", i18n(s$en, s$zh))
        )
      })
      n_done <- sum(vapply(ph_steps, function(s) isTRUE(status[[s$v]]), logical(1)))
      has_current <- any(vapply(ph_steps, function(s) identical(s$v, current),
                                logical(1)))
      children[[length(children) + 1]] <- shiny::tags$details(
        class = "scstudio-phasegroup",
        open = if (has_current || !length(children)) NA else NULL,
        shiny::tags$summary(
          class = "scstudio-phase",
          i18n(lab$en, lab$zh),
          shiny::span(class = "scstudio-phase-count",
                      sprintf("%d/%d", n_done, length(ph_steps)))
        ),
        items
      )
    }
    shiny::div(class = "scstudio-nav", children)
  })

  shiny::observeEvent(input$goto, { bslib::nav_select("steps", input$goto) })

  # --- Progress chip (topbar): done/total + the next unfinished step ---------
  # The "next step" link lives here, in the always-visible top bar, rather than
  # floating over the workspace where it covered the figure-export row.
  output$progress_chip <- shiny::renderUI({
    steps <- app_steps()
    n <- length(steps)
    keys <- vapply(steps, function(s) s$v, character(1))
    is_done <- vapply(keys, function(k) step_state(rv, k) == "done", logical(1))
    done <- sum(is_done)
    nxt <- NULL
    cur <- input$steps
    if (!is.null(cur) && isTRUE(is_done[cur]) && !all(is_done)) {
      i <- match(cur, keys)
      undone <- which(!is_done)
      j <- undone[undone > i][1]
      if (is.na(j)) j <- undone[1]
      s <- steps[[j]]
      nxt <- shiny::tags$a(
        class = "scstudio-nextlink",
        onclick = sprintf("Shiny.setInputValue('goto','%s',{priority:'event'})", s$v),
        shiny::span(class = "scstudio-next-label", i18n("Next:", "下一步：")),
        shiny::strong(i18n(s$en, s$zh)), " \u2192")
    }
    shiny::div(class = "scstudio-progress",
               title = i18n("Steps completed", "已完成步骤"),
               shiny::span(class = "scstudio-progress-num",
                           sprintf("%d/%d", done, n)),
               shiny::div(class = "scstudio-progress-bar",
                          shiny::div(class = "scstudio-progress-fill",
                                     style = sprintf("width:%.0f%%", 100 * done / n))),
               nxt)
  })

  # --- Banner above the workspace when the current step is stale ------------
  output$stale_banner <- shiny::renderUI(stale_banner(rv, input$steps))

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
  mod_clinical_server("clinical", rv, log_rv)
  mod_viz_server("viz", rv, log_rv)
  mod_report_server("report", rv, log_rv)
  mod_export_server("export", rv, log_rv)
}
