#' Launch the scStudio interactive analysis app
#'
#' Starts a local Shiny server and (by default) opens it in your browser. All
#' computation runs on your own machine; nothing is uploaded to any server.
#'
#' Typical use:
#' \preformatted{
#'   scStudio::run_app()
#' }
#'
#' @param launch.browser Logical; open the app in a browser automatically.
#'   Default `TRUE`.
#' @param port Integer port, or NULL to pick a free one.
#' @param host Host to bind. Default "127.0.0.1" (localhost only). Use
#'   "0.0.0.0" when running inside Docker so the browser on the host can reach it.
#' @param max_upload_mb Maximum upload size in MB for the file inputs. Default
#'   5120 (5 GB) so users can upload real count matrices via the browser.
#' @return Called for its side effect (runs the app). Does not return.
#' @export
run_app <- function(launch.browser = TRUE,
                    port = getOption("shiny.port"),
                    host = "127.0.0.1",
                    max_upload_mb = 5120) {
  # Allow large browser uploads (default shiny cap is only 5 MB).
  old <- options(shiny.maxRequestSize = max_upload_mb * 1024^2)
  on.exit(options(old), add = TRUE)

  app <- shiny::shinyApp(ui = app_ui, server = app_server)
  shiny::runApp(app, launch.browser = launch.browser, port = port, host = host)
}

#' Internal: path to an app resource under inst/
#' @param ... Path components under inst/.
#' @keywords internal
app_sys <- function(...) {
  system.file(..., package = "scStudio")
}
