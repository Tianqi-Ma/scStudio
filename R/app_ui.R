#' Top-level UI: a clean, modern, wizard-style navbar
#'
#' The app is organised as an ordered set of steps (Import -> QC -> Doublet ->
#' Normalize -> Features/PCA -> Integrate -> Cluster -> Embed -> Markers ->
#' Annotate -> Visualize -> Export). Each step is a self-contained module that
#' follows the same pattern (explainer, method choice, adjustable thresholds,
#' run, summary, interactive preview with hover).
#'
#' @return A [bslib::page_navbar()] UI definition.
#' @keywords internal
app_ui <- function() {
  theme <- bslib::bs_theme(
    version = 5,
    preset = "shiny",
    primary = "#3b6ea5",
    "border-radius" = "0.6rem",
    base_font = bslib::font_google("Inter", local = FALSE),
    heading_font = bslib::font_google("Inter", local = FALSE)
  )

  bslib::page_navbar(
    title = shiny::tags$span(class = "scstudio-brand",
                             shiny::tags$strong("scStudio"),
                             shiny::tags$span(class = "scstudio-sub", "single-cell, locally")),
    id = "main_nav",
    theme = theme,
    fillable = TRUE,
    header = shiny::tagList(
      shiny::tags$head(
        shiny::tags$link(rel = "stylesheet", type = "text/css", href = "www/custom.css")
      )
    ),
    sidebar = bslib::sidebar(
      title = "Dataset",
      width = 280,
      open = "desktop",
      shiny::uiOutput("global_status"),
      shiny::hr(),
      shiny::div(class = "scstudio-hint",
                 shiny::tags$p("Work through the steps left to right."),
                 shiny::tags$p("Each step explains itself, lets you pick a method and thresholds, then shows a preview."))
    ),

    bslib::nav_panel("1. Import",    shiny::icon("upload"),      mod_import_ui("import")),
    bslib::nav_panel("2. QC",        shiny::icon("filter"),      mod_qc_ui("qc")),
    bslib::nav_panel("3. Doublets",  shiny::icon("clone"),       mod_doublet_ui("doublet")),
    bslib::nav_panel("4. Normalize", shiny::icon("wave-square"), mod_normalize_ui("normalize")),
    bslib::nav_panel("5. Features/PCA", shiny::icon("chart-line"), mod_reduce_ui("reduce")),
    bslib::nav_panel("6. Integrate", shiny::icon("layer-group"), mod_integrate_ui("integrate")),
    bslib::nav_panel("7. Cluster",   shiny::icon("circle-nodes"),mod_cluster_ui("cluster")),
    bslib::nav_panel("8. Embed",     shiny::icon("braille"),     mod_embed_ui("embed")),
    bslib::nav_panel("9. Markers",   shiny::icon("dna"),         mod_markers_ui("markers")),
    bslib::nav_panel("10. Annotate", shiny::icon("tags"),        mod_annotate_ui("annotate")),
    bslib::nav_panel("11. Visualize",shiny::icon("chart-simple"),mod_viz_ui("viz")),
    bslib::nav_panel("12. Export",   shiny::icon("download"),    mod_export_ui("export")),

    bslib::nav_spacer(),
    bslib::nav_item(
      shiny::tags$a(shiny::icon("circle-info"), " Help",
                    href = "https://github.com/Tianqi-Ma/scStudio",
                    target = "_blank")
    )
  )
}
