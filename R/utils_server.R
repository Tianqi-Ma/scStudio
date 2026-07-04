#' Reusable server-side helpers for scStudio
#'
#' @name utils_server
#' @keywords internal
NULL

#' Check that suggested packages are installed, notify if not
#'
#' Heavy compute packages (Seurat, scDblFinder, harmony, SingleR, ...) live in
#' Suggests so the app installs and the UI loads even without them. Each compute
#' step calls this first and refuses to run (with a helpful install hint) if a
#' required package is missing.
#'
#' @param pkgs Character vector of package names required by the step.
#' @param what Character. Short description of the step, used in the message.
#' @return `TRUE` if all present; otherwise `FALSE` (and shows a notification).
#' @keywords internal
require_pkgs <- function(pkgs, what = "this step") {
  missing <- pkgs[!vapply(pkgs, function(p) {
    requireNamespace(p, quietly = TRUE)
  }, logical(1))]
  if (length(missing) == 0) return(TRUE)
  msg <- sprintf(
    "%s needs package(s) not installed: %s. Install them, then retry.",
    what, paste(missing, collapse = ", ")
  )
  if (shiny::isRunning()) {
    shiny::showNotification(msg, type = "error", duration = 10)
  } else {
    warning(msg, call. = FALSE)
  }
  FALSE
}

#' Is a package available?
#' @param pkg Package name.
#' @keywords internal
has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

#' Run an expression with a Shiny progress bar and graceful error notify
#'
#' Wraps a (possibly slow) compute call so the UI shows progress and any error
#' surfaces as a notification instead of crashing the session.
#'
#' @param expr Expression to evaluate.
#' @param message Progress message shown to the user.
#' @param session Shiny session (for progress).
#' @return The value of `expr`, or `NULL` on error.
#' @keywords internal
with_progress_notify <- function(expr, message = "Working...", session = shiny::getDefaultReactiveDomain()) {
  prog <- shiny::Progress$new(session)
  prog$set(message = message, value = 0.1)
  on.exit(prog$close(), add = TRUE)
  out <- tryCatch({
    prog$set(value = 0.5)
    force(expr)
  }, error = function(e) {
    shiny::showNotification(paste("Error:", conditionMessage(e)),
                            type = "error", duration = 12)
    NULL
  })
  prog$set(value = 1)
  out
}

#' Append a reproducibility log entry (step + parameters + equivalent R code)
#'
#' @param log_rv A `reactiveVal` holding a list of log entries.
#' @param step Character step name.
#' @param params Named list of chosen parameters.
#' @param code Character vector of equivalent R code lines.
#' @keywords internal
log_step <- function(log_rv, step, params = list(), code = character(0)) {
  entry <- list(
    step = step,
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    params = params,
    code = code
  )
  cur <- log_rv()
  cur[[length(cur) + 1]] <- entry
  log_rv(cur)
  invisible(entry)
}

#' Guess the sample/batch column in an object's metadata
#'
#' @param meta A data.frame of cell metadata.
#' @return A best-guess column name, or NULL.
#' @keywords internal
guess_batch_col <- function(meta) {
  if (is.null(meta) || ncol(meta) == 0) return(NULL)
  candidates <- c("sample", "orig.ident", "batch", "donor", "patient",
                  "Sample", "Batch", "sampleID", "sample_id")
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit)) return(hit[1])
  NULL
}

#' Human-readable size of the current dataset for memory guardrails
#'
#' @param n_cells Integer number of cells.
#' @return A short advisory string, or "" if within comfortable range.
#' @keywords internal
memory_advice <- function(n_cells) {
  if (is.null(n_cells) || is.na(n_cells)) return("")
  if (n_cells > 5e5) {
    return(paste0("This dataset has ", format(n_cells, big.mark = ","),
                  " cells. Expect high memory use (>32 GB) and slow steps; ",
                  "consider downsampling for exploration."))
  }
  if (n_cells > 1e5) {
    return(paste0(format(n_cells, big.mark = ","),
                  " cells: plan for 16-32 GB RAM; heavy steps may take minutes."))
  }
  ""
}
