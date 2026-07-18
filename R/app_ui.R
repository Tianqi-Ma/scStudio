#' Top-level UI: dark, plot-first, left step navigator
#'
#' Layout: a collapsible left sidebar lists the 12 steps, each coloured by status
#' (blue = done, highlighted = current, grey = not yet run). The main area shows
#' one step at a time via a hidden tabset. Dark theme by default with a light
#' toggle, and an English/中文 language switch (client-side, no round-trip).
#'
#' @return A [bslib::page_sidebar()] UI definition.
#' @keywords internal
app_ui <- function() {
  theme <- bslib::bs_theme(
    version = 5,
    preset = "shiny",
    primary = "#4f8fd0",
    "border-radius" = "0.6rem",
    base_font = bslib::font_google("Inter", local = FALSE),
    heading_font = bslib::font_google("Inter", local = FALSE)
  )

  # The ordered steps: value key, icon, and bilingual label.
  steps <- list(
    list(v = "import",    i = "upload",       en = "Import",       zh = "导入"),
    list(v = "qc",        i = "filter",       en = "QC",           zh = "质控"),
    list(v = "doublet",   i = "clone",        en = "Doublets",     zh = "去双细胞"),
    list(v = "normalize", i = "wave-square",  en = "Normalize",    zh = "归一化"),
    list(v = "reduce",    i = "chart-line",   en = "Features/PCA", zh = "特征/PCA"),
    list(v = "integrate", i = "layer-group",  en = "Integrate",    zh = "整合"),
    list(v = "cluster",   i = "circle-nodes", en = "Cluster",      zh = "聚类"),
    list(v = "embed",     i = "braille",      en = "Embed",        zh = "降维图"),
    list(v = "markers",   i = "dna",          en = "Markers",      zh = "标志基因"),
    list(v = "annotate",  i = "tags",         en = "Annotate",     zh = "注释"),
    list(v = "viz",       i = "chart-simple", en = "Visualize",    zh = "可视化"),
    list(v = "export",    i = "download",     en = "Export",       zh = "导出")
  )
  ui_of <- list(
    import = mod_import_ui("import"), qc = mod_qc_ui("qc"),
    doublet = mod_doublet_ui("doublet"), normalize = mod_normalize_ui("normalize"),
    reduce = mod_reduce_ui("reduce"), integrate = mod_integrate_ui("integrate"),
    cluster = mod_cluster_ui("cluster"), embed = mod_embed_ui("embed"),
    markers = mod_markers_ui("markers"), annotate = mod_annotate_ui("annotate"),
    viz = mod_viz_ui("viz"), export = mod_export_ui("export")
  )
  panels <- lapply(steps, function(s) {
    bslib::nav_panel(title = s$en, value = s$v, ui_of[[s$v]])
  })
  main <- do.call(bslib::navset_hidden, c(list(id = "steps"), panels))

  bslib::page_sidebar(
    theme = theme,
    title = shiny::div(
      class = "scstudio-topbar",
      shiny::span(class = "scstudio-brand",
                  shiny::strong("scStudio"),
                  shiny::span(class = "scstudio-sub", "single-cell, locally")),
      shiny::div(
        class = "scstudio-topright",
        bslib::input_switch("lang_zh", "中文", value = FALSE),
        bslib::input_dark_mode(id = "dark", mode = "dark")
      )
    ),
    sidebar = bslib::sidebar(
      title = i18n("Steps", "分析步骤"),
      width = 220, open = "open", id = "stepbar",
      shiny::uiOutput("step_nav"),
      shiny::hr(),
      shiny::div(class = "scstudio-mini", shiny::uiOutput("global_status"))
    ),

    # head: stylesheet + client-side language toggle (no server round-trip)
    shiny::tags$head(
      shiny::tags$link(rel = "stylesheet", type = "text/css", href = "www/custom.css"),
      shiny::tags$script(shiny::HTML(
        "$(document).on('shiny:inputchanged', function(e){ if(e.name==='lang_zh'){ document.body.classList.toggle('lang-zh', !!e.value); } });"
      ))
    ),
    main
  )
}
