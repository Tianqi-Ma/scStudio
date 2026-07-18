#' Plotting layer — wraps scop's plotting functions for a unified look
#'
#' All previews use scop's plotting functions (CellDimPlot, FeatureDimPlot,
#' GroupHeatmap, DynamicHeatmap, VolcanoPlot, EnrichmentPlot, ...) with
#' `scop::palette_scp()` colours, so every figure matches the scop/SCP aesthetic.
#' Each wrapper is gated by `require_pkgs("scop")` and wrapped in tryCatch so a
#' missing package or a signature mismatch surfaces as a friendly message rather
#' than crashing the app.
#'
#' NOTE: scop signatures are verified at runtime on the user's machine; a few
#' argument names may need adjustment against the installed scop version.
#'
#' @name fct_plots
#' @keywords internal
NULL

#' Shared minimal ggplot theme (fallback when not using a scop plot)
#' @keywords internal
scstudio_theme <- function() {
  ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "#8b98a533"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}

#' Categorical palette — scop's palette_scp when available, else a curated set
#' @param n Number of colours.
#' @param type "discrete", "continuous", or "diverging".
#' @keywords internal
sc_palette <- function(n = 8, type = "discrete") {
  if (has_pkg("scop")) {
    out <- tryCatch(scop::palette_scp(seq_len(n), n = n, type = type),
                    error = function(e) NULL)
    if (!is.null(out)) return(unname(out))
  }
  base <- c("#2f81c7", "#e4572e", "#3fb37f", "#b5179e", "#f4a261", "#4361ee",
            "#e63946", "#2a9d8f", "#9c6ade", "#ffca3a", "#577590", "#d68fb0",
            "#43aa8b", "#f9844a", "#277da1", "#f94144", "#90be6d", "#845ec2",
            "#ff9f1c", "#4d908e", "#c9184a", "#00b4d8", "#bc6c25", "#606c38")
  if (n <= length(base)) return(base[seq_len(n)])
  grDevices::colorRampPalette(base)(n)
}

#' UI slot for a (large) preview plot. scop plots are static high-res images.
#' @param id Namespaced output id. @param height CSS height.
#' @keywords internal
preview_plot_ui <- function(id, height = "100%") {
  shiny::plotOutput(id, height = height)
}

#' Render a scop/ggplot/ComplexHeatmap object to a Shiny plot output
#'
#' Accepts whatever a scop plotting function returns: a ggplot/patchwork object
#' (printed) or a ComplexHeatmap (drawn). `plot_expr` is a function returning the
#' plot object; it runs inside tryCatch so failures show as a message.
#' @keywords internal
render_scop_plot <- function(plot_expr) {
  shiny::renderPlot({
    # Build the plot object. req()/validate() (no data yet) stay silent.
    p <- tryCatch(
      plot_expr(),
      shiny.silent.error = function(e) NULL,
      error = function(e) structure(list(msg = conditionMessage(e)), class = "scstudio_plot_error"))
    shiny::req(!is.null(p))

    # Draw it. Any drawing error is turned into a readable message ON the canvas
    # (and a toast) so the user never sees an opaque "[object Object]".
    show_err <- function(msg) {
      shiny::showNotification(paste("Plot error:", msg), type = "error", duration = 12)
      op <- graphics::par(mar = c(0, 0, 0, 0)); on.exit(graphics::par(op), add = TRUE)
      graphics::plot.new()
      graphics::text(0.5, 0.5, paste0("Plot error:\n", msg), col = "#c1476b", cex = 1.1)
    }
    if (inherits(p, "scstudio_plot_error")) { show_err(p$msg); return(invisible()) }
    tryCatch({
      if (methods::is(p, "Heatmap") || methods::is(p, "HeatmapList")) {
        if (has_pkg("ComplexHeatmap")) ComplexHeatmap::draw(p) else print(p)
      } else {
        print(p)
      }
    }, error = function(e) show_err(conditionMessage(e)))
  })
}

# Backward-compatible alias used by older modules (renders a ggplot expr).
#' @keywords internal
render_preview_plot <- function(gg_expr, tooltip = "text") {
  render_scop_plot(gg_expr)
}

# ---- scop plotting wrappers -------------------------------------------------

#' Dimensional-reduction scatter (clusters / metadata), scop::CellDimPlot,
#' with optional mascarade cell-type outlines.
#' @keywords internal
sc_dimplot <- function(srt, group_by, reduction = NULL, mask = FALSE,
                       palette = "Paired", label = TRUE, ...) {
  if (!require_pkgs("scop", "Dimension plot")) return(NULL)
  p <- scop::CellDimPlot(srt, group.by = group_by, reduction = reduction,
                         palette = palette, label = label, ...)
  if (isTRUE(mask) && has_pkg("mascarade")) {
    p <- tryCatch(add_mascarade(p, srt, group_by, reduction),
                  error = function(e) p)
  }
  p
}

#' Feature (gene / score) on a reduction, scop::FeatureDimPlot
#' @keywords internal
sc_featureplot <- function(srt, features, reduction = NULL, ...) {
  if (!require_pkgs("scop", "Feature plot")) return(NULL)
  scop::FeatureDimPlot(srt, features = features, reduction = reduction, ...)
}

#' Grouped mean-expression heatmap, scop::GroupHeatmap (signature figure)
#' @keywords internal
sc_groupheatmap <- function(srt, features, group_by, ...) {
  if (!require_pkgs("scop", "GroupHeatmap")) return(NULL)
  scop::GroupHeatmap(srt, features = features, group.by = group_by, ...)
}

#' Dynamic (pseudotime) heatmap, scop::DynamicHeatmap
#' @keywords internal
sc_dynamicheatmap <- function(srt, lineages, ...) {
  if (!require_pkgs("scop", "DynamicHeatmap")) return(NULL)
  scop::DynamicHeatmap(srt, lineages = lineages, ...)
}

#' Composition / statistics plot, scop::CellStatPlot
#' @keywords internal
sc_cellstat <- function(srt, stat_by, group_by = NULL, plot_type = "bar", ...) {
  if (!require_pkgs("scop", "Cell statistics")) return(NULL)
  scop::CellStatPlot(srt, stat.by = stat_by, group.by = group_by,
                     plot_type = plot_type, ...)
}

#' Per-group feature distribution, scop::FeatureStatPlot
#' @keywords internal
sc_featurestat <- function(srt, stat_by, group_by, plot_type = "violin", ...) {
  if (!require_pkgs("scop", "Feature statistics")) return(NULL)
  scop::FeatureStatPlot(srt, stat.by = stat_by, group.by = group_by,
                        plot_type = plot_type, ...)
}

#' Volcano plot of DE results, scop::VolcanoPlot
#' @keywords internal
sc_volcano <- function(srt, group_by, ...) {
  if (!require_pkgs("scop", "Volcano plot")) return(NULL)
  scop::VolcanoPlot(srt, group_by = group_by, ...)
}

#' Enrichment plot (GO/KEGG/...), scop::EnrichmentPlot
#' @keywords internal
sc_enrichplot <- function(srt, group_by, plot_type = "bar", ...) {
  if (!require_pkgs("scop", "Enrichment plot")) return(NULL)
  scop::EnrichmentPlot(srt, group_by = group_by, plot_type = plot_type, ...)
}

#' GSEA running-score plot, scop::GSEAPlot
#' @keywords internal
sc_gseaplot <- function(srt, ...) {
  if (!require_pkgs("scop", "GSEA plot")) return(NULL)
  scop::GSEAPlot(srt, ...)
}

#' RNA-velocity stream/grid, scop::VelocityPlot
#' @keywords internal
sc_velocityplot <- function(srt, reduction = NULL, ...) {
  if (!require_pkgs("scop", "Velocity plot")) return(NULL)
  scop::VelocityPlot(srt, reduction = reduction, ...)
}

#' PAGA graph on an embedding, scop::PAGAPlot
#' @keywords internal
sc_pagaplot <- function(srt, ...) {
  if (!require_pkgs("scop", "PAGA plot")) return(NULL)
  scop::PAGAPlot(srt, ...)
}

#' Mascarade cell-type outlines overlaid on a scop dim plot
#'
#' Uses mascarade::generateMask() to compute polygon outlines around each group
#' on the 2D embedding and overlays them on an existing ggplot dim plot.
#' @keywords internal
add_mascarade <- function(p, srt, group_by, reduction = NULL) {
  if (!require_pkgs("mascarade", "Cell-type outlines")) return(p)
  emb <- SeuratObject::Embeddings(srt, reduction = reduction %||% SeuratObject::DefaultDimReduc(srt))
  labels <- obj_meta(srt)[[group_by]]
  mask <- mascarade::generateMask(dims = emb[, 1:2], cluster = labels)
  p + ggplot2::geom_path(
    data = mask,
    ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
    colour = "grey20", linewidth = 0.4, inherit.aes = FALSE
  )
}
