# Smoke tests for the app shell: the UI must build and the full server graph
# (all 12 modules) must initialise without data and without the heavy Suggests
# packages installed. These guard against wiring regressions.

test_that("app_ui builds and renders to HTML", {
  ui <- app_ui()
  html <- as.character(ui)
  expect_true(nchar(html) > 1000)
  # the twelve numbered steps are present in the navbar
  for (lbl in c("Import", "QC", "Doublets", "Normalize", "Cluster",
                "Embed", "Markers", "Annotate", "Visualize", "Export")) {
    expect_match(html, lbl, fixed = TRUE)
  }
})

test_that("full server graph initialises with no data loaded", {
  # Should not error even though Seurat/plotly/etc. may be absent, because
  # compute is gated behind Run buttons and outputs handle the NULL object.
  expect_no_error(
    shiny::testServer(app_server, {
      session$flushReact()
    })
  )
})

test_that("import module server initialises", {
  rv <- shiny::reactiveValues(obj = NULL, source = NULL)
  log_rv <- shiny::reactiveVal(list())
  expect_no_error(
    shiny::testServer(mod_import_server, args = list(rv = rv, log_rv = log_rv), {
      session$flushReact()
    })
  )
})
