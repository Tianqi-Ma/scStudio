#' Module: Export report
#'
#' Generate a narrated, scCancer-style HTML report of the analysis: the methods
#' used, their parameters and equivalent code, assembled from the reproducibility
#' log. No heavy compute happens here -- it only reads what previous steps logged.
#'
#' @param id Module id. @param rv shared hub. @param log_rv repro log.
#' @name mod_report
NULL

# Report sections offered to the user. Each maps to keyword(s) matched (case-
# insensitively) against a logged step's name, so the report only narrates the
# sections the user asked for.
.report_sections <- list(
  qc            = list(en = "Quality control",   zh = "质量控制",   kw = c("qc", "quality")),
  doublet       = list(en = "Doublets",          zh = "双细胞",     kw = c("doublet")),
  normalize     = list(en = "Normalization",     zh = "标准化",     kw = c("normal")),
  reduce        = list(en = "Features & PCA",    zh = "特征与 PCA", kw = c("pca", "feature", "reduc")),
  integrate     = list(en = "Integration",       zh = "整合",       kw = c("integr")),
  cluster       = list(en = "Clustering",        zh = "聚类",       kw = c("cluster")),
  embed         = list(en = "Embedding",         zh = "降维可视化", kw = c("embed", "umap", "tsne")),
  markers       = list(en = "Marker genes",      zh = "标志基因",   kw = c("marker")),
  annotation    = list(en = "Annotation",        zh = "注释",       kw = c("annot")),
  trajectory    = list(en = "Trajectory",        zh = "拟时序",     kw = c("trajector", "pseudotime", "dynamic")),
  enrichment    = list(en = "Enrichment / GSEA", zh = "富集 / GSEA", kw = c("enrich", "gsea")),
  malignancy    = list(en = "Malignant / CNV",   zh = "恶性 / CNV", kw = c("cnv", "malignant", "stemness"))
)

#' @rdname mod_report
#' @keywords internal
mod_report_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Export report", zh = "导出报告"),
    what = list(
      en = "Generate a narrated HTML report of your analysis: the methods,
            figures, and parameters for each step you ran.",
      zh = "生成一份叙述式 HTML 分析报告：包含你运行的每一步的方法、图表和参数。"),
    why  = list(
      en = "A shareable, human-readable summary (scCancer-style) documents what
            was done and how -- for collaborators, supervisors, or a methods
            section.",
      zh = "一份可分享、易读的摘要（scCancer 风格）记录了做了什么以及如何做 —— 供合作者、导师或方法学部分使用。"),
    how  = list(
      en = "Tick the sections to include and give the report a title, then click
            <b>Download report</b>. The report is built from the reproducibility
            log, so it always matches what you actually ran.",
      zh = "勾选要包含的章节并为报告命名，然后点击<b>下载报告</b>。报告依据可复现日志生成，因此始终与你实际运行的步骤一致。"),
    example = list(
      en = "A report titled <i>PBMC 3k analysis</i> with QC, Clustering and
               Annotation sections, ready to attach to an email.",
      zh = "一份题为 <i>PBMC 3k analysis</i> 的报告，包含质量控制、聚类和注释章节，可直接作为邮件附件。")
  )
  controls <- shiny::tagList(
    label_with_help("Report title", "Shown as the report heading.",
                    label_zh = "报告标题", tip_zh = "作为报告的标题显示。"),
    shiny::textInput(ns("title"), NULL, value = "scStudio analysis report"),
    label_with_help("Sections to include",
                    "Only the ticked sections are narrated (matched against the steps you ran).",
                    label_zh = "包含的章节",
                    tip_zh = "只有勾选的章节会被写入报告（与你运行过的步骤匹配）。"),
    shiny::checkboxGroupInput(
      ns("sections"), NULL,
      choices  = stats::setNames(names(.report_sections),
                                 vapply(.report_sections, `[[`, character(1), "en")),
      selected = names(.report_sections)),
    shiny::downloadButton(ns("download_report"),
                          i18n("Download report", "下载报告"), class = "w-100")
  )
  step_container(
    title     = list(en = "Export report", zh = "导出报告"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = shiny::uiOutput(ns("preview"))
  )
}

#' @rdname mod_report
#' @keywords internal
mod_report_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {

    # Keep only the logged steps that match the ticked sections. A step whose
    # name matches no known section keyword is always kept (never silently lost).
    filter_entries <- function(entries, sections) {
      if (is.null(entries) || !length(entries)) return(entries)
      all_kw  <- unlist(lapply(.report_sections, `[[`, "kw"))
      sel_kw  <- unlist(lapply(.report_sections[sections], `[[`, "kw"))
      Filter(function(e) {
        s <- tolower(e$step)
        known <- any(vapply(all_kw, function(k) grepl(k, s, fixed = TRUE), logical(1)))
        selected <- any(vapply(sel_kw, function(k) grepl(k, s, fixed = TRUE), logical(1)))
        selected || !known
      }, entries)
    }

    # Format one step's parameters as "k=v, k=v".
    params_text <- function(params) {
      if (!length(params)) return("")
      paste(vapply(names(params), function(k)
        sprintf("%s = %s", k, paste(deparse(params[[k]]), collapse = "")),
        character(1)), collapse = ", ")
    }

    # Assemble the report body as an R Markdown string (used by rmarkdown, or
    # rendered to a self-contained HTML fallback).
    build_rmd <- function(title, entries) {
      lines <- c(
        "---",
        sprintf('title: "%s"', gsub('"', "'", title)),
        sprintf('date: "%s"', format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        "output:",
        "  html_document:",
        "    self_contained: true",
        "    toc: true",
        "---",
        "",
        "> Generated by scStudio from the reproducibility log.",
        "")
      if (is.null(entries) || !length(entries)) {
        return(paste(c(lines, "No analysis steps were recorded yet."), collapse = "\n"))
      }
      body <- unlist(lapply(seq_along(entries), function(i) {
        e <- entries[[i]]
        p <- params_text(e$params)
        c(sprintf("## %d. %s", i, e$step),
          sprintf("*Run at %s.*", e$time %||% "unknown time"),
          "",
          if (nzchar(p)) c("**Parameters:**", "", sprintf("`%s`", p), "") else NULL,
          if (length(e$code)) c("**Code:**", "", "```r", e$code, "```", "") else NULL)
      }))
      paste(c(lines, body), collapse = "\n")
    }

    # Self-contained HTML fallback when rmarkdown/pandoc is unavailable.
    build_html <- function(title, entries) {
      esc <- function(x) {
        x <- gsub("&", "&amp;", x, fixed = TRUE)
        x <- gsub("<", "&lt;",  x, fixed = TRUE)
        gsub(">", "&gt;", x, fixed = TRUE)
      }
      css <- paste(
        "body{font-family:system-ui,Arial,sans-serif;max-width:860px;margin:2rem auto;",
        "padding:0 1rem;color:#1f2933;line-height:1.6}",
        "h1{border-bottom:2px solid #3b6ea5;padding-bottom:.3rem}",
        "h2{color:#3b6ea5;margin-top:2rem}",
        ".muted{color:#66788a;font-size:.9em}",
        "pre{background:#f4f6f8;padding:.8rem;border-radius:6px;overflow:auto}",
        "code{background:#f4f6f8;padding:.1rem .3rem;border-radius:4px}",
        collapse = "")
      head <- c(
        "<!DOCTYPE html><html><head><meta charset='utf-8'>",
        sprintf("<title>%s</title>", esc(title)),
        sprintf("<style>%s</style></head><body>", css),
        sprintf("<h1>%s</h1>", esc(title)),
        sprintf("<p class='muted'>Generated by scStudio &middot; %s</p>",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      if (is.null(entries) || !length(entries)) {
        body <- "<p>No analysis steps were recorded yet.</p>"
      } else {
        body <- unlist(lapply(seq_along(entries), function(i) {
          e <- entries[[i]]
          p <- params_text(e$params)
          c(sprintf("<h2>%d. %s</h2>", i, esc(e$step)),
            sprintf("<p class='muted'>Run at %s.</p>", esc(e$time %||% "unknown time")),
            if (nzchar(p)) sprintf("<p><b>Parameters:</b> <code>%s</code></p>", esc(p)) else "",
            if (length(e$code))
              sprintf("<p><b>Code:</b></p><pre><code>%s</code></pre>",
                      esc(paste(e$code, collapse = "\n"))) else "")
        }))
      }
      paste(c(head, body, "</body></html>"), collapse = "\n")
    }

    output$summary <- shiny::renderUI({
      entries <- log_rv()
      if (is.null(entries) || !length(entries)) {
        return(shiny::div(class = "scstudio-placeholder",
                          i18n("No steps recorded yet. Run some analysis steps first.",
                               "尚未记录任何步骤。请先运行一些分析步骤。")))
      }
      shiny::tagList(
        stat_tile(i18n("Steps performed", "已执行步骤数"), length(entries)),
        shiny::tags$ol(class = "scstudio-steps",
          lapply(entries, function(e) {
            shiny::tags$li(shiny::tags$b(e$step),
                           shiny::tags$span(class = "scstudio-muted",
                                            paste0("  (", e$time, ")")))
          }))
      )
    })

    output$preview <- shiny::renderUI({
      shiny::div(class = "scstudio-note",
                 i18n(paste0("No plot for this step. Tick the sections to include, ",
                             "set a title, then click Download report to save a ",
                             "narrated HTML summary of your analysis."),
                      paste0("这一步没有图表。请勾选要包含的章节、设置标题，",
                             "然后点击“下载报告”以保存一份叙述式的 HTML 分析摘要。")))
    })

    output$download_report <- shiny::downloadHandler(
      filename = function()
        paste0("scstudio_report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html"),
      content = function(file) {
        title   <- if (nzchar(trimws(input$title %||% ""))) input$title
                   else "scStudio analysis report"
        entries <- filter_entries(log_rv(), input$sections)

        # Prefer a proper rmarkdown render; fall back to a self-contained HTML
        # string if rmarkdown/pandoc is unavailable or rendering fails.
        rendered <- FALSE
        if (has_pkg("rmarkdown") && rmarkdown::pandoc_available()) {
          rendered <- tryCatch({
            rmd <- tempfile(fileext = ".Rmd")
            writeLines(build_rmd(title, entries), rmd)
            out <- rmarkdown::render(rmd, output_format = "html_document",
                                     output_file = basename(tempfile(fileext = ".html")),
                                     output_dir = dirname(rmd), quiet = TRUE)
            file.copy(out, file, overwrite = TRUE)
            TRUE
          }, error = function(e) FALSE)
        }
        if (!rendered) {
          writeLines(build_html(title, entries), file)
        }
      }
    )
  })
}
