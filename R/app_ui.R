#' Top-level UI: dark, plot-first, grouped left stepper, splash animation
#'
#' Layout: a top bar (brand, EN/中 segmented switch, dark/light toggle, Python-env
#' status), a collapsible left navigator grouped by analysis phase with per-step
#' status colouring, and a main area showing one step at a time. A startup splash
#' animation (cells coalescing into a UMAP) plays on load.
#'
#' @return A [bslib::page_sidebar()] UI wrapped with the splash overlay.
#' @keywords internal
app_ui <- function() {
  theme <- bslib::bs_theme(
    version = 5, preset = "shiny",
    primary = "#2f81c7",
    base_font = bslib::font_google("Inter", local = FALSE),
    heading_font = bslib::font_google("Inter", local = FALSE)
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
      shiny::div(
        class = "scstudio-topright",
        shiny::div(class = "scstudio-pydot", id = "py-status",
                   i18n("Python: not set up", "Python：未配置")),
        shiny::tags$div(
          class = "scstudio-lang",
          shiny::tags$button(class = "scstudio-lang-btn active", `data-lang` = "en",
                             onclick = "scStudioSetLang('en')", "EN"),
          shiny::tags$button(class = "scstudio-lang-btn", `data-lang` = "zh",
                             onclick = "scStudioSetLang('zh')", "中")
        ),
        bslib::input_dark_mode(id = "dark", mode = "dark")
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
      shiny::tags$script(src = "scstudio/app.js")
    ),
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

#' Ordered pipeline steps with phase grouping and module UI functions
#'
#' Single source of truth for the navigator and the hidden tabset. New modules
#' are added here as they are implemented.
#' @keywords internal
app_steps <- function() {
  list(
    list(v = "import",    n = 1,  phase = "data",     en = "Import",        zh = "导入",       ui = mod_import_ui),
    list(v = "qc",        n = 2,  phase = "data",     en = "Quality control", zh = "质控",     ui = mod_qc_ui),
    list(v = "doublet",   n = 3,  phase = "data",     en = "Doublets",      zh = "去双细胞", ui = mod_doublet_ui),
    list(v = "normalize", n = 4,  phase = "preproc",  en = "Normalize",     zh = "归一化", ui = mod_normalize_ui),
    list(v = "reduce",    n = 5,  phase = "preproc",  en = "Features / PCA",zh = "特征/PCA",   ui = mod_reduce_ui),
    list(v = "integrate", n = 6,  phase = "preproc",  en = "Integrate",     zh = "整合",       ui = mod_integrate_ui),
    list(v = "cluster",   n = 7,  phase = "structure",en = "Cluster",       zh = "聚类",       ui = mod_cluster_ui),
    list(v = "embed",     n = 8,  phase = "structure",en = "Embed",         zh = "降维图", ui = mod_embed_ui),
    list(v = "markers",   n = 9,  phase = "identity", en = "Markers",       zh = "标志基因", ui = mod_markers_ui),
    list(v = "annotate",  n = 10, phase = "identity", en = "Annotate",      zh = "注释",       ui = mod_annotate_ui),
    list(v = "viz",       n = 11, phase = "output",   en = "Visualize",     zh = "可视化", ui = mod_viz_ui),
    list(v = "export",    n = 12, phase = "output",   en = "Export",        zh = "导出",       ui = mod_export_ui)
  )
}

#' Phase labels (en/zh) for the grouped navigator
#' @keywords internal
app_phases <- function() {
  list(
    data      = list(en = "Data & QC",    zh = "数据与质控"),
    preproc   = list(en = "Preprocess",   zh = "预处理"),
    structure = list(en = "Structure",    zh = "结构"),
    identity  = list(en = "Identity",     zh = "身份"),
    output    = list(en = "Output",       zh = "产出")
  )
}
