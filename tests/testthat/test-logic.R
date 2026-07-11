# Unit tests for the pure-logic helpers (no Seurat / no network required).
# These run anywhere, including CI without the heavy Suggests packages.

test_that("cell-type dictionary explains common labels (and tolerates variants)", {
  expect_match(explain_celltype("T cell"), "immune", ignore.case = TRUE)
  # tolerant of markers / plurals / case
  expect_match(explain_celltype("CD8+ T cells"), "[Cc]ytotoxic")
  expect_equal(explain_celltype("totally-unknown-type-xyz"), "")
})

test_that("guess_batch_col finds common sample columns", {
  expect_equal(guess_batch_col(data.frame(orig.ident = 1, x = 2)), "orig.ident")
  expect_equal(guess_batch_col(data.frame(sample = 1)), "sample")
  expect_null(guess_batch_col(data.frame(a = 1, b = 2)))
  expect_null(guess_batch_col(data.frame()))
})

test_that("memory_advice warns only for large datasets", {
  expect_equal(memory_advice(1e4), "")
  expect_gt(nchar(memory_advice(2e5)), 0)
  expect_gt(nchar(memory_advice(2e6)), 0)
})

test_that("demo dataset catalogue is well-formed", {
  ds <- demo_datasets()
  expect_s3_class(ds, "data.frame")
  expect_true(all(c("id", "name", "format", "url", "description") %in% names(ds)))
  expect_true("bundled" %in% ds$id)
  # the bundled entry needs no URL; network ones must have one
  expect_equal(ds$url[ds$id == "bundled"], "")
  expect_true(all(nzchar(ds$url[ds$id != "bundled"])))
})

test_that("MAD-based QC keeps normal cells and flags outliers", {
  set.seed(1)
  md <- data.frame(
    nCount_RNA   = c(rpois(60, 1000), 50000, 30),   # last two: huge / tiny
    nFeature_RNA = c(rpois(60, 500),  20,    15),
    percent.mt   = c(runif(60, 1, 8), 60,    55)     # last two: high mito
  )
  keep <- qc_mad_keep_from_meta(md, nmads_lib = 5, nmads_mt = 3)
  expect_length(keep, 62)
  expect_false(keep[61])          # extreme high-count / high-mito cell
  expect_false(keep[62])          # extreme low-count cell
  expect_gt(mean(keep[1:60]), 0.9)  # the bulk of normal cells survive
})

test_that("manual QC thresholds exclude out-of-range cells", {
  md <- data.frame(
    nCount_RNA   = c(1000, 1000, 1000),
    nFeature_RNA = c(500,  50,   500),   # 2nd below min
    percent.mt   = c(5,    5,    40)      # 3rd above max
  )
  keep <- qc_manual_keep_from_meta(md, min_genes = 100, max_genes = 2000, max_mt = 15)
  expect_equal(keep, c(TRUE, FALSE, FALSE))
})

test_that("read_counts_table parses a simple matrix", {
  tf <- tempfile(fileext = ".tsv")
  on.exit(unlink(tf))
  writeLines(c("gene\tc1\tc2", "GENE1\t3\t0", "GENE2\t1\t5"), tf)
  m <- read_counts_table(tf, sep = "\t")
  expect_equal(dim(m), c(2, 2))
  expect_equal(rownames(m), c("GENE1", "GENE2"))
  expect_equal(unname(m["GENE2", "c2"]), 5)
})
