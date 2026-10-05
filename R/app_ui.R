#' Top-level UI: plot-first, grouped left stepper, splash animation
#'
#' Layout: a top bar (brand, progress + next step, export menu, EN/中 switch,
#' light/dark toggle), a collapsible left navigator grouped by analysis phase
#' with per-step status colouring, and a main area showing one step at a time.
#' A startup splash animation (cells coalescing into a UMAP) plays on load.
#'
#' @return A [bslib::page_sidebar()] UI wrapped with the splash overlay.
#' @keywords internal
app_ui <- function() {
  # Inter is linked from Google Fonts, but the app must look right offline and
  # where Google is unreachable: the browser then falls through the system UI
  # stack (CJK-capable on every OS) instead of a default serif face.
  ui_font <- app_font_stack()

  theme <- bslib::bs_theme(
    version = 5, preset = "shiny",
    primary = "#2f81c7",
    base_font = ui_font,
    heading_font = ui_font
  )

  panels <- lapply(app_steps(), function(s) {
    bslib::nav_panel(title = s$en, value = s$v, do.call(s$ui, list(s$v)))
  })
  main <- do.call(bslib::navset_hidden, c(list(id = "steps"), panels))

  page <- bslib::page_sidebar(
    theme = theme,
    title = shiny::div(
      class = "scstudio-topbar",
      shiny::span(class = "scstudio-brand",
                  shiny::strong("scStudio"),
                  shiny::span(class = "scstudio-sub", "single-cell, locally")),
      shiny::uiOutput("progress_chip", inline = TRUE),
      shiny::div(
        class = "scstudio-topright",
        bslib::popover(
          shiny::actionLink("export_menu", i18n("⤓ Export", "⤓ 导出"),
                            class = "scstudio-switch"),
          title = i18n("Export current data", "导出当前数据"),
          shiny::uiOutput("export_items")
        ),
        shiny::tags$div(
          class = "scstudio-lang",
          shiny::tags$button(class = "scstudio-lang-btn active", `data-lang` = "en",
                             onclick = "scStudioSetLang('en')", "EN"),
          shiny::tags$button(class = "scstudio-lang-btn", `data-lang` = "zh",
                             onclick = "scStudioSetLang('zh')", "中")
        ),
        # follow the operating system's light/dark preference; the toggle
        # still switches it for the session
        bslib::input_dark_mode(id = "dark")
      )
    ),
    sidebar = bslib::sidebar(
      title = i18n("Workflow", "分析流程"),
      width = 232, open = "open", id = "stepbar",
      shiny::uiOutput("step_nav"),
      shiny::hr(),
      shiny::div(class = "scstudio-mini", shiny::uiOutput("global_status"))
    ),
    shiny::tags$head(
      shiny::tags$link(rel = "stylesheet", type = "text/css", href = "scstudio/custom.css"),
      shiny::tags$script(src = "scstudio/app.js"),
      shiny::tags$script(src = "scstudio/explain.js"),
      shiny::tags$script(src = "scstudio/explain-sc.js")
    ),
    shiny::uiOutput("stale_banner"),
    main
  )

  # Splash overlay (position:fixed) + the app page.
  shiny::tagList(
    shiny::div(
      id = "scstudio-splash",
      shiny::tags$canvas(id = "scstudio-splash-canvas"),
      shiny::div(id = "scstudio-splash-logo",
                 shiny::div(class = "t", "scStudio"),
                 shiny::div(class = "s", "single-cell analysis, on your own machine"))
    ),
    page
  )
}

#' Ordered pipeline steps (the shared single-cell registry, see steps_sc.R)
#' @keywords internal
app_steps <- function() steps_sc()

#' Phase labels (en/zh) for the grouped navigator
#' @keywords internal
app_phases <- function() phases_sc()

#' Every registry, keyed by omics (used by the pipeline-state helpers)
#' @keywords internal
all_step_registries <- function() list(sc = steps_sc())
