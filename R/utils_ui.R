#' Reusable UI building blocks for scStudio
#'
#' These helpers implement the uniform per-step interaction pattern used by every
#' analysis module: a beginner-friendly explainer card, a `?` help tooltip, and a
#' two-column "controls / preview" step container. Keeping them in one place keeps
#' every module visually consistent (clean, modern, card-based).
#'
#' @name utils_ui
#' @keywords internal
NULL

#' Collapsible beginner explainer card ("What is this step?")
#'
#' Renders a collapsible card that explains, in plain language, what a step does,
#' why it matters, how to read the result, and how to choose parameters -- with a
#' concrete example. Beginners expand it; experts collapse it.
#'
#' @param title Character. Short step title (e.g. "Doublet removal").
#' @param what,why,how,example Character (may contain inline HTML). The four
#'   explanation facets. `example` is highlighted as a worked example.
#' @param open Logical. Whether the card starts expanded. Default `TRUE`.
#' @return A [bslib::accordion()] element.
#' @keywords internal
explainer_card <- function(title, what, why, how = NULL, example = NULL,
                           open = TRUE) {
  body <- list(
    shiny::tags$p(shiny::tags$strong("What: "), shiny::HTML(what)),
    shiny::tags$p(shiny::tags$strong("Why: "), shiny::HTML(why))
  )
  if (!is.null(how)) {
    body <- c(body, list(shiny::tags$p(shiny::tags$strong("How to read / choose: "),
                                       shiny::HTML(how))))
  }
  if (!is.null(example)) {
    body <- c(body, list(
      shiny::div(class = "scstudio-example",
                 shiny::tags$span(class = "scstudio-example-tag", "Example"),
                 shiny::HTML(example))
    ))
  }
  bslib::accordion(
    class = "scstudio-explainer",
    open = open,
    bslib::accordion_panel(
      title = shiny::tags$span(
        shiny::tags$span("\U0001F4A1 ", `aria-hidden` = "true"),
        paste0("What is this step? — ", title)
      ),
      value = paste0("explainer-", gsub("[^a-z0-9]+", "-", tolower(title))),
      icon = NULL,
      body
    )
  )
}

#' Inline help tooltip on a `?` icon
#'
#' @param text Tooltip text (plain language definition of a term).
#' @return A [bslib::tooltip()]-wrapped `?` badge.
#' @keywords internal
help_tip <- function(text) {
  bslib::tooltip(
    shiny::tags$span(class = "scstudio-help", "?"),
    text,
    placement = "right"
  )
}

#' A labelled control with an attached help tooltip
#'
#' @param label Character label shown before the `?`.
#' @param tip Character tooltip text.
#' @return A shiny tag with label + help icon (place above an input).
#' @keywords internal
label_with_help <- function(label, tip) {
  shiny::tags$label(class = "scstudio-label", label, help_tip(tip))
}

#' Standard two-column step container (controls left, results right)
#'
#' Implements the uniform module layout: explainer on top, then a sidebar of
#' method/threshold controls + a "Run" button, and a main area holding the result
#' summary and the interactive preview plot(s).
#'
#' @param explainer An [explainer_card()] (or NULL).
#' @param controls Tag list of inputs (method selectors, sliders, the Run button).
#' @param summary A UI output slot for the textual/numeric result summary.
#' @param preview A UI output slot for the preview plot(s).
#' @param width_sidebar Sidebar width (CSS unit). Default "340px".
#' @return A [bslib::card()] with a [bslib::layout_sidebar()] inside.
#' @keywords internal
step_container <- function(explainer, controls, summary, preview,
                           width_sidebar = "340px") {
  bslib::card(
    class = "scstudio-step",
    full_screen = TRUE,
    if (!is.null(explainer)) bslib::card_body(explainer, class = "scstudio-explainer-wrap"),
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width = width_sidebar,
        title = "Settings",
        controls
      ),
      bslib::card_body(
        class = "scstudio-results",
        shiny::div(class = "scstudio-summary", summary),
        shiny::div(class = "scstudio-preview", preview)
      )
    )
  )
}

#' A prominent primary "Run" action button used by every module
#'
#' @param id Input id (already namespaced by the caller's `ns()`).
#' @param label Button label. Default "Run this step".
#' @keywords internal
run_button <- function(id, label = "Run this step") {
  shiny::actionButton(
    id, label,
    icon = shiny::icon("play"),
    class = "btn-primary scstudio-run w-100"
  )
}

#' Small helper to render a value_box-style summary tile
#'
#' @param title,value Character/numeric shown in the tile.
#' @param showcase Optional icon.
#' @keywords internal
stat_tile <- function(title, value, showcase = NULL) {
  bslib::value_box(
    title = title,
    value = value,
    showcase = showcase,
    class = "scstudio-tile"
  )
}
