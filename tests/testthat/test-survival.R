# The survival layer is pure data.frame maths, so it can be tested properly
# (unlike the scop/Seurat wrappers, which need packages that are not installed
# here). These tests pin the behaviour the modules rely on.

make_cohort <- function(n = 80, seed = 42) {
  set.seed(seed)
  data.frame(
    sample    = paste0("S", seq_len(n)),
    os_months = round(stats::rexp(n, 1 / 24), 1),
    os_status = sample(c("Dead", "Alive"), n, TRUE, prob = c(0.6, 0.4)),
    stage     = sample(c("I", "II", "III"), n, TRUE),
    age       = round(stats::rnorm(n, 62, 10)),
    stringsAsFactors = FALSE
  )
}

test_that("encode_event handles the codings people actually use", {
  expect_equal(encode_event(c(0, 1, 1, 0)), c(0L, 1L, 1L, 0L))
  expect_equal(encode_event(c(TRUE, FALSE)), c(1L, 0L))
  # survival's other convention: 1 = censored, 2 = event
  expect_equal(encode_event(c(1, 2, 2, 1)), c(0L, 1L, 1L, 0L))
  expect_equal(encode_event(c("Dead", "Alive", "yes", "no")), c(1L, 0L, 1L, 0L))
  expect_equal(encode_event(c("DECEASED", " living ")), c(1L, 0L))
  # unrecognised coding -> NA, so the module can ask which level means "event"
  expect_true(all(is.na(encode_event(c("groupA", "groupB")))))
  # ...and an explicit choice resolves it
  expect_equal(encode_event(c("groupA", "groupB"), positive = "groupA"), c(1L, 0L))
})

test_that("normalise_clinical converts time units and drops unusable rows", {
  cl <- make_cohort(20)
  d <- normalise_clinical(cl, "sample", "os_months", "os_status")
  expect_true(all(c(".id", ".time", ".event") %in% names(d)))
  expect_equal(nrow(d), 20)
  expect_equal(attr(d, "dropped"), 0)
  expect_true(all(d$.event %in% c(0L, 1L)))

  # days -> months
  cl$os_days <- cl$os_months * 30.4375
  dd <- normalise_clinical(cl, "sample", "os_days", "os_status", time_unit = "days")
  expect_equal(dd$.time, d$.time, tolerance = 1e-6)

  # years -> months
  cl$os_years <- cl$os_months / 12
  dy <- normalise_clinical(cl, "sample", "os_years", "os_status", time_unit = "years")
  expect_equal(dy$.time, d$.time, tolerance = 1e-6)

  # bad rows are removed, not passed to survival to fail on
  bad <- cl
  bad$os_months[1:3] <- c(NA, -5, NA)
  bad$os_status[4] <- "unknown-code"
  db <- normalise_clinical(bad, "sample", "os_months", "os_status")
  expect_equal(nrow(db), 16)
  expect_equal(attr(db, "dropped"), 4)

  expect_error(normalise_clinical(cl, "nope", "os_months", "os_status"),
               "Column not found")
})

test_that("split_numeric produces the groups each method promises", {
  x <- 1:100
  m <- split_numeric(x, "median")
  expect_equal(levels(m), c("Low", "High"))
  expect_equal(as.integer(table(m)), c(50L, 50L))

  tt <- split_numeric(x, "tertile")
  expect_true(sum(is.na(tt)) > 0)              # middle third is dropped
  expect_equal(as.integer(table(tt)), c(34L, 34L))

  set.seed(7)
  time  <- c(stats::rexp(50, 1 / 5), stats::rexp(50, 1 / 40))
  event <- rep(1L, 100)
  score <- c(stats::runif(50, 0, 1), stats::runif(50, 1, 2))   # high score = long survival
  op <- split_numeric(score, "optimal", time, event)
  expect_equal(levels(op), c("Low", "High"))
  # the min_frac guard must hold: neither arm may collapse
  expect_true(min(table(op)) >= 20)
  expect_true(!is.null(attr(op, "cutpoint")))

  # too few distinct values -> falls back to the median split rather than erroring
  expect_equal(levels(split_numeric(c(1, 1, 2), "optimal", c(1, 2, 3), c(1, 1, 0))),
               c("Low", "High"))
  expect_error(split_numeric(1:10, "optimal"), "needs the survival time")
})

test_that("composition_by_sample returns fractions that sum to 1 per sample", {
  meta <- data.frame(
    smp = rep(c("A", "B", "C"), times = c(10, 20, 30)),
    ct  = c(rep("T", 5), rep("B", 5),
            rep("T", 10), rep("Mono", 10),
            rep("T", 15), rep("B", 5), rep("Mono", 10)),
    stringsAsFactors = FALSE)
  comp <- composition_by_sample(meta, "smp", "ct")
  expect_equal(comp$.id, c("A", "B", "C"))
  expect_equal(comp$n_obs, c(10L, 20L, 30L))
  fr <- comp[, setdiff(names(comp), c(".id", "n_obs"))]
  expect_equal(unname(rowSums(fr)), rep(1, 3))
  expect_equal(comp$T[1], 0.5)
  expect_error(composition_by_sample(meta, "nope", "ct"), "not found")
})

test_that("km_fit / logrank_test / km_medians agree on a two-arm cohort", {
  set.seed(3)
  n <- 60
  d <- data.frame(
    .id    = paste0("S", seq_len(n)),
    .time  = c(stats::rexp(n / 2, 1 / 6), stats::rexp(n / 2, 1 / 30)),
    .event = rep(1L, n),
    .group = factor(rep(c("Low", "High"), each = n / 2), levels = c("Low", "High")),
    stringsAsFactors = FALSE)

  fit <- km_fit(d)
  expect_s3_class(fit, "survfit")

  lr <- logrank_test(d)
  expect_true(is.list(lr) && all(c("chisq", "df", "p") %in% names(lr)))
  expect_equal(lr$df, 1)
  expect_lt(lr$p, 0.001)                       # arms differ five-fold; must separate

  med <- km_medians(fit)
  expect_equal(nrow(med), 2)
  expect_setequal(med$group, c("Low", "High"))
  expect_equal(sum(med$n), n)
  expect_gt(med$median[med$group == "High"], med$median[med$group == "Low"])

  # one arm only -> nothing to compare, and that must be NULL not an error
  d1 <- d; d1$.group <- factor("only")
  expect_null(logrank_test(d1))
  expect_s3_class(km_fit(d1[, setdiff(names(d1), ".group")]), "survfit")
})

test_that("cox_univariable screens each variable on its own", {
  cl <- make_cohort(90)
  d <- normalise_clinical(cl, "sample", "os_months", "os_status")
  out <- cox_univariable(d, c("stage", "age"))
  expect_true(all(c("variable", "level", "n", "HR", "lower", "upper", "p") %in% names(out)))
  expect_setequal(unique(out$variable), c("stage", "age"))
  expect_equal(nrow(out[out$variable == "stage", ]), 2)   # 3 levels -> 2 contrasts
  expect_true(all(out$HR > 0))
  expect_true(all(out$lower <= out$HR & out$HR <= out$upper))
  expect_false(is.unsorted(out$p))                        # ordered by p

  # a constant variable cannot be modelled: skipped, not fatal
  d$constant <- "same"
  expect_equal(nrow(cox_univariable(d, "constant")), 0)
  expect_s3_class(cox_univariable(d, "constant"), "data.frame")
})

test_that("km_tidy anchors every curve at (0, 1)", {
  set.seed(5)
  d <- data.frame(.time = stats::rexp(40, 1 / 10), .event = rbinom(40, 1, 0.7),
                  .group = factor(rep(c("a", "b"), 20)))
  td <- km_tidy(km_fit(d))
  starts <- td[td$time == 0, ]
  expect_equal(nrow(starts), 2)
  expect_true(all(starts$surv == 1))
  expect_true(all(td$surv >= 0 & td$surv <= 1))
  # monotone non-increasing within each curve
  for (g in unique(td$group)) {
    s <- td$surv[td$group == g]
    expect_false(is.unsorted(rev(s)))
  }
})
