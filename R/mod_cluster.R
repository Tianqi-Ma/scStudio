#' Module: Clustering
#'
#' Group cells into clusters (candidate cell populations) with a
#' neighbor-graph + community-detection approach, at one or more resolutions.
#' Leiden (algorithm 4) is the modern default; Louvain (algorithm 1) is classic.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_cluster
NULL

#' @rdname mod_cluster
#' @keywords internal
mod_cluster_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = "Clustering",
    what = "Partition cells into clusters that likely correspond to distinct
            cell types or states.",
    why  = "Clusters are the units you annotate and compare. Good clustering
            separates real populations without over-splitting noise.",
    how  = "Higher <b>resolution</b> = more, smaller clusters. Try several
            resolutions and compare. Use the reduction you want to cluster on
            (PCA, or an integrated embedding like Harmony).",
    example = "At resolution 0.2 you may get 6 broad clusters; at 1.0 they split
               into finer subtypes.<br><b>Note:</b> Leiden (algorithm 4) needs the
               Python <code>leidenalg</code> (or <code>leidenbase</code>) backend;
               if unavailable, switch to Louvain."
  )
  controls <- shiny::tagList(
    label_with_help("Method",
                    "Leiden (algorithm 4) is the modern default; Louvain (algorithm 1) is classic."),
    shiny::selectInput(ns("method"), NULL,
                       c("Leiden (algorithm 4)" = "leiden",
                         "Louvain (algorithm 1)" = "louvain")),
    label_with_help("Reduction",
                    "Which dimensional reduction to build the neighbor graph on."),
    shiny::uiOutput(ns("reduction_ui")),
    label_with_help("Dimensions", "Number of leading dimensions to use (e.g. PCs)."),
    shiny::numericInput(ns("dims"), NULL, value = 30, min = 2, max = 100),
    label_with_help("Neighbors (k)", "Neighbors used to build the graph."),
    shiny::numericInput(ns("neighbors"), NULL, value = 20, min = 2, max = 100),
    label_with_help("Resolutions",
                    "Comma-separated list; each is clustered. The last one is used for the summary/preview."),
    shiny::textInput(ns("resolutions"), NULL, value = "0.2,0.5,0.8,1.0"),
    run_button(ns("run"), "Run clustering")
  )
  step_container(explainer, controls,
                 summary = shiny::uiOutput(ns("summary")),
                 preview = preview_plot_ui(ns("preview")))
}

#' @rdname mod_cluster
#' @keywords internal
mod_cluster_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- shiny::reactiveValues(done = FALSE, col = NULL, res_used = NULL,
                                 n_clusters = NA_integer_)

    # Offer reductions suitable for clustering (pca / harmony if present).
    output$reduction_ui <- shiny::renderUI({
      reds <- obj_reductions(rv$obj)
      choices <- intersect(c("pca", "harmony"), reds)
      if (length(choices) == 0) choices <- reds
      if (length(choices) == 0) {
        return(shiny::div(class = "scstudio-placeholder",
                          "No reductions yet — run PCA first."))
      }
      selected <- if ("harmony" %in% choices) "harmony" else choices[1]
      shiny::selectInput(session$ns("reduction"), NULL,
                         choices = choices, selected = selected)
    })

    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      shiny::req(input$reduction)
      if (!require_pkgs("Seurat", "Clustering")) return(NULL)
      # Parse the comma-separated resolutions into a numeric vector.
      resolutions <- suppressWarnings(as.numeric(
        trimws(strsplit(input$resolutions, ",", fixed = TRUE)[[1]])))
      resolutions <- resolutions[!is.na(resolutions)]
      if (length(resolutions) == 0) {
        shiny::showNotification("Enter at least one valid resolution.",
                                type = "error")
        return(NULL)
      }
      algorithm <- if (input$method == "leiden") 4 else 1
      reduction <- input$reduction
      dims <- input$dims
      obj <- with_progress_notify({
        cluster_obj(rv$obj, reduction = reduction, dims = dims,
                    resolutions = resolutions, algorithm = algorithm)
      }, message = "Building graph and clustering...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      # Seurat stores the last FindClusters result in `seurat_clusters`.
      last_res <- resolutions[length(resolutions)]
      clusters <- obj_meta(obj)$seurat_clusters
      res$done       <- TRUE
      res$col        <- "seurat_clusters"
      res$res_used   <- last_res
      res$n_clusters <- if (is.null(clusters)) NA_integer_ else nlevels(factor(clusters))
      log_step(log_rv, "Clustering",
               params = list(method = input$method, algorithm = algorithm,
                             reduction = reduction, dims = dims,
                             resolutions = resolutions),
               code = sprintf(
                 "obj <- cluster_obj(obj, reduction = '%s', dims = %d, resolutions = c(%s), algorithm = %d)",
                 reduction, dims, paste(resolutions, collapse = ", "), algorithm))
      shiny::showNotification(
        sprintf("Clustering done: %d clusters at resolution %s.",
                res$n_clusters, format(last_res)),
        type = "message")
    })

    output$summary <- shiny::renderUI({
      if (!isTRUE(res$done)) {
        return(shiny::div(class = "scstudio-placeholder",
                          "Set parameters and click Run clustering."))
      }
      bslib::layout_columns(
        col_widths = c(6, 6),
        stat_tile("Clusters", format(res$n_clusters)),
        stat_tile("Resolution", format(res$res_used))
      )
    })

    output$preview <- render_preview_plot(function() {
      shiny::req(res$done)
      md <- obj_meta(rv$obj)
      shiny::req(res$col %in% colnames(md))
      cl <- factor(md[[res$col]])
      counts <- as.data.frame(table(cluster = cl), stringsAsFactors = FALSE)
      names(counts) <- c("cluster", "n")
      counts$cluster <- factor(counts$cluster, levels = levels(cl))
      counts$text <- sprintf("Cluster %s\n%s cells",
                             counts$cluster, format(counts$n, big.mark = ","))
      ggplot2::ggplot(counts, ggplot2::aes(x = cluster, y = n,
                                           fill = cluster, text = text)) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_manual(values = scstudio_palette(nlevels(counts$cluster)),
                                   guide = "none") +
        ggplot2::labs(x = "Cluster", y = "Cells",
                      title = sprintf("Cluster sizes (resolution %s)",
                                      format(res$res_used))) +
        scstudio_theme()
    })
  })
}
