#' Plotting helpers: uniform theme + interactive hover
#'
#' Every preview plot uses the same clean theme and, when {plotly} is available,
#' becomes interactive with hover tooltips (the app's hover-preview feature). If
#' plotly is not installed the app falls back to a static ggplot so nothing breaks.
#'
#' @name fct_plots
#' @keywords internal
NULL

#' The shared ggplot theme (clean, minimal, modern)
#' @keywords internal
scstudio_theme <- function() {
  ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "#e6e6e6"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      strip.text = ggplot2::element_text(face = "bold")
    )
}

#' A brand-neutral categorical palette
#' @param n Number of colours needed.
#' @keywords internal
scstudio_palette <- function(n) {
  base <- c("#3b6ea5", "#e07b39", "#4c9f70", "#c1476b", "#8a6fb0",
            "#d9a441", "#5aa9c9", "#a15a4e", "#7d8b8f", "#b05fa0")
  if (n <= length(base)) return(base[seq_len(n)])
  grDevices::colorRampPalette(base)(n)
}

#' UI output slot for a preview plot (interactive if plotly present)
#' @param id Namespaced output id.
#' @param height CSS height.
#' @keywords internal
preview_plot_ui <- function(id, height = "420px") {
  if (has_pkg("plotly")) {
    plotly::plotlyOutput(id, height = height)
  } else {
    shiny::plotOutput(id, height = height)
  }
}

#' Render a ggplot as an interactive hover plot (or static fallback)
#'
#' @param gg A ggplot object. Map hover text via `text = ` inside `aes()` and it
#'   will be shown on hover when plotly is available.
#' @param tooltip Character vector of aesthetics to show on hover
#'   (default "text").
#' @return A render function suitable for assigning to `output$id`.
#' @keywords internal
render_preview_plot <- function(gg_expr, tooltip = "text") {
  if (has_pkg("plotly")) {
    plotly::renderPlotly({
      gg <- gg_expr()
      shiny::req(gg)
      plotly::ggplotly(gg, tooltip = tooltip) |>
        plotly::config(displayModeBar = FALSE)
    })
  } else {
    shiny::renderPlot({
      gg <- gg_expr()
      shiny::req(gg)
      gg
    })
  }
}
