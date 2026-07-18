# Smoke tests for the app shell: the UI must build and the full server graph
# (all 12 modules) must initialise without data and without the heavy Suggests
# packages installed. These guard against wiring regressions.

test_that("app_ui builds and renders to HTML", {
  ui <- app_ui()
  rt <- htmltools::renderTags(ui)
  html <- paste(as.character(rt$html), as.character(rt$head))
  expect_true(nchar(html) > 2000)
  # chrome: stylesheet + script served from the registered resource path,
  # the splash overlay, the language switch, and the bilingual data attributes
  expect_match(html, "scstudio/custom.css", fixed = TRUE)
  expect_match(html, "scstudio/app.js", fixed = TRUE)
  expect_match(html, "scstudio-splash", fixed = TRUE)
  expect_match(html, "scStudioSetLang", fixed = TRUE)
  expect_match(html, "data-zh=", fixed = TRUE)
})

test_that("app_steps and app_phases are consistent", {
  steps <- app_steps()
  phases <- app_phases()
  expect_true(length(steps) >= 12)
  # every step's phase must have a label
  for (s in steps) {
    expect_true(!is.null(phases[[s$phase]]), info = s$v)
    expect_true(is.function(s$ui))
  }
  # step values are unique
  vals <- vapply(steps, function(s) s$v, character(1))
  expect_equal(anyDuplicated(vals), 0L)
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
