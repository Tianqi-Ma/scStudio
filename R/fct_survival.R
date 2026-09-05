#' Shared survival / prognosis layer
#'
#' Omics-agnostic survival analysis used by every pipeline that has clinical
#' follow-up: single-cell (cell-type composition vs outcome), WES (mutation vs
#' outcome), bulk (molecular subtype vs outcome) and multi-omics integration
#' (integrative cluster vs outcome). The normalised clinical table lives on the
#' shared hub as `rv$clinical`, so it is loaded once and reused everywhere.
#'
#' Everything here is a pure function of data frames -- no Seurat, no Shiny --
#' so it is directly testable, and the maths uses `survival`, which ships with R
#' and therefore installs on any machine (no compiled GitHub packages, no
#' `survminer`; the Kaplan-Meier curve is drawn with ggplot2 here).
#'
#' @name fct_survival
#' @keywords internal
NULL

# ---- clinical table ---------------------------------------------------------

#' Read a clinical table (csv/tsv) with the sample ids as an ordinary column
#' @param path File path. @param sep Field separator.
#' @keywords internal
read_clinical_table <- function(path, sep = ",") {
  utils::read.delim(path, sep = sep, check.names = FALSE,
                    stringsAsFactors = FALSE)
}

#' Normalise a clinical table to the columns the survival code expects
#'
#' Produces `.id` / `.time` / `.event` alongside the untouched original columns,
#' so downstream code never has to care what the user's columns were called.
#' Rows with a missing/negative time or an unusable event value are dropped --
#' `survival` would otherwise fail with a much less obvious message.
#'
#' @param df Raw clinical data.frame.
#' @param id_col,time_col,event_col Column names chosen by the user.
#' @param event_positive Value of `event_col` that means "event occurred"
#'   (death / progression). Numeric 1 and the usual text codings are detected
#'   automatically when this is NULL.
#' @param time_unit One of "days", "months", "years": the unit of `time_col`.
#'   Times are converted to months, which is what the plots label.
#' @return data.frame with `.id`, `.time` (months), `.event` (0/1) plus the
#'   original columns; attribute "dropped" records how many rows were removed.
#' @keywords internal
normalise_clinical <- function(df, id_col, time_col, event_col,
                               event_positive = NULL,
                               time_unit = c("months", "days", "years")) {
  time_unit <- match.arg(time_unit)
  stopifnot(is.data.frame(df))
  for (cl in c(id_col, time_col, event_col)) {
    if (!cl %in% names(df)) stop("Column not found in the clinical table: ", cl)
  }

  time <- suppressWarnings(as.numeric(df[[time_col]]))
  time <- switch(time_unit, months = time, days = time / 30.4375, years = time * 12)

  event <- encode_event(df[[event_col]], event_positive)

  out <- df
  out$.id    <- as.character(df[[id_col]])
  out$.time  <- time
  out$.event <- event

  keep <- !is.na(out$.id) & !is.na(out$.time) & out$.time >= 0 & !is.na(out$.event)
  dropped <- sum(!keep)
  out <- out[keep, , drop = FALSE]
  attr(out, "dropped") <- dropped
  out
}

#' Coerce an event column to 0/1
#'
#' Accepts 0/1, TRUE/FALSE, 1/2 (the `survival` convention), and the common text
#' codings. `positive` forces a specific level to mean "event".
#'
#' @param x The raw event column. @param positive Level meaning "event", or NULL.
#' @return Integer vector of 0/1 (NA where undecidable).
#' @keywords internal
encode_event <- function(x, positive = NULL) {
  if (!is.null(positive) && nzchar(as.character(positive))) {
    return(as.integer(as.character(x) == as.character(positive)))
  }
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) {
    u <- sort(unique(stats::na.omit(x)))
    # survival's other convention is 1 = censored, 2 = event
    if (length(u) == 2 && all(u == c(1, 2))) return(as.integer(x == 2))
    return(as.integer(x > 0))
  }
  s <- tolower(trimws(as.character(x)))
  dead  <- c("1", "true", "yes", "y", "dead", "deceased", "death", "event",
             "progressed", "progression", "recurrence", "relapse")
  alive <- c("0", "false", "no", "n", "alive", "living", "censored", "censor",
             "no event", "disease free", "disease-free")
  out <- rep(NA_integer_, length(s))
  out[s %in% dead]  <- 1L
  out[s %in% alive] <- 0L
  out
}

# ---- grouping ---------------------------------------------------------------

#' Split a numeric variable into survival groups
#'
#' @param x Numeric vector (e.g. a cell-type proportion or an expression score).
#' @param method "median", "tertile" (top vs bottom third, middle dropped), or
#'   "optimal" (cutpoint maximising the log-rank statistic).
#' @param time,event Needed only by "optimal".
#' @param min_frac Smallest share of samples allowed in a group by "optimal".
#' @return A factor the same length as `x`; NA marks samples not assigned.
#' @keywords internal
split_numeric <- function(x, method = c("median", "tertile", "optimal"),
                          time = NULL, event = NULL, min_frac = 0.2) {
  method <- match.arg(method)
  x <- as.numeric(x)

  if (method == "median") {
    cut <- stats::median(x, na.rm = TRUE)
    return(factor(ifelse(x > cut, "High", "Low"), levels = c("Low", "High")))
  }

  if (method == "tertile") {
    q <- stats::quantile(x, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
    g <- rep(NA_character_, length(x))
    g[x <= q[1]] <- "Low"
    g[x >= q[2]] <- "High"
    return(factor(g, levels = c("Low", "High")))
  }

  # "optimal": scan candidate cutpoints, keep the one with the largest log-rank
  # chi-square. Done by hand rather than pulling in maxstat, and deliberately
  # bounded by min_frac so it cannot "find" a split of three samples.
  if (is.null(time) || is.null(event)) {
    stop("The optimal cutpoint needs the survival time and event columns.")
  }
  ok <- !is.na(x) & !is.na(time) & !is.na(event)
  cand <- sort(unique(x[ok]))
  if (length(cand) < 3) return(split_numeric(x, "median"))
  lo <- stats::quantile(x[ok], min_frac, na.rm = TRUE)
  hi <- stats::quantile(x[ok], 1 - min_frac, na.rm = TRUE)
  cand <- cand[cand >= lo & cand < hi]
  if (!length(cand)) return(split_numeric(x, "median"))

  best <- NULL; best_stat <- -Inf
  for (cut in cand) {
    g <- factor(ifelse(x > cut, "High", "Low"), levels = c("Low", "High"))
    st <- tryCatch(
      survival::survdiff(survival::Surv(time[ok], event[ok]) ~ g[ok])$chisq,
      error = function(e) NA_real_)
    if (!is.na(st) && st > best_stat) { best_stat <- st; best <- cut }
  }
  if (is.null(best)) return(split_numeric(x, "median"))
  out <- factor(ifelse(x > best, "High", "Low"), levels = c("Low", "High"))
  attr(out, "cutpoint") <- best
  out
}

#' Per-sample composition (fraction of each group within each sample)
#'
#' The bridge from a per-observation table to the per-sample table survival
#' analysis needs: for single-cell that is cells -> cell-type fractions per
#' sample; the same function serves any per-observation categorical breakdown.
#'
#' @param meta data.frame with one row per observation (e.g. per cell).
#' @param sample_col Column holding the sample/patient id.
#' @param group_col Column holding the category (cluster, cell type, ...).
#' @return data.frame: `.id`, one column per category (fractions summing to 1),
#'   and `n_obs`.
#' @keywords internal
composition_by_sample <- function(meta, sample_col, group_col) {
  if (!all(c(sample_col, group_col) %in% names(meta))) {
    stop("Sample or group column not found in the object metadata.")
  }
  s <- as.character(meta[[sample_col]])
  g <- as.character(meta[[group_col]])
  ok <- !is.na(s) & !is.na(g)
  tab <- table(s[ok], g[ok])
  frac <- as.data.frame.matrix(tab / rowSums(tab))
  frac$.id <- rownames(frac)
  frac$n_obs <- as.integer(rowSums(tab))
  rownames(frac) <- NULL
  frac[, c(".id", setdiff(names(frac), c(".id", "n_obs")), "n_obs"), drop = FALSE]
}

# ---- models -----------------------------------------------------------------

#' Kaplan-Meier fit for a normalised clinical table
#' @param df Output of [normalise_clinical()], optionally with a `.group` column.
#' @keywords internal
km_fit <- function(df) {
  if (".group" %in% names(df) && length(unique(stats::na.omit(df$.group))) > 1) {
    survival::survfit(survival::Surv(.time, .event) ~ .group, data = df)
  } else {
    survival::survfit(survival::Surv(.time, .event) ~ 1, data = df)
  }
}

#' Log-rank test across `.group`
#' @param df Normalised clinical table with a `.group` column.
#' @return list(chisq, df, p) or NULL when there is nothing to compare.
#' @keywords internal
logrank_test <- function(df) {
  if (!".group" %in% names(df)) return(NULL)
  d <- df[!is.na(df$.group), , drop = FALSE]
  if (length(unique(d$.group)) < 2) return(NULL)
  sd <- tryCatch(survival::survdiff(survival::Surv(.time, .event) ~ .group, data = d),
                 error = function(e) NULL)
  if (is.null(sd)) return(NULL)
  k <- length(sd$n) - 1
  list(chisq = unname(sd$chisq), df = k,
       p = stats::pchisq(sd$chisq, df = k, lower.tail = FALSE))
}

#' Median survival per group (months)
#' @param fit A [km_fit()] result.
#' @return data.frame: group, n, events, median, lower, upper.
#' @keywords internal
km_medians <- function(fit) {
  s <- summary(fit)$table
  if (is.null(dim(s))) s <- t(as.matrix(s))
  gp <- rownames(s)
  if (is.null(gp)) gp <- "All"
  pick <- function(nm) if (nm %in% colnames(s)) unname(s[, nm]) else rep(NA_real_, nrow(s))
  data.frame(
    group  = sub("^\\.group=", "", gp),
    n      = pick("records"),
    events = pick("events"),
    median = pick("median"),
    lower  = pick("0.95LCL"),
    upper  = pick("0.95UCL"),
    stringsAsFactors = FALSE
  )
}

#' Univariable Cox regression, one model per variable
#'
#' Univariable (not multivariable) on purpose: it is the screen users want at
#' this point, and it stays interpretable when the cohort is small.
#'
#' @param df Normalised clinical table.
#' @param vars Column names to test.
#' @return data.frame: variable, level, n, HR, lower, upper, p (ordered by p).
#' @keywords internal
cox_univariable <- function(df, vars) {
  rows <- lapply(vars, function(v) {
    d <- df[, c(".time", ".event", v), drop = FALSE]
    names(d)[3] <- "x"
    if (is.character(d$x)) d$x <- factor(d$x)
    d <- d[!is.na(d$x), , drop = FALSE]
    if (nrow(d) < 5) return(NULL)
    if (is.factor(d$x) && nlevels(droplevels(d$x)) < 2) return(NULL)
    if (is.factor(d$x)) d$x <- droplevels(d$x)
    fit <- tryCatch(survival::coxph(survival::Surv(.time, .event) ~ x, data = d),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    sm <- summary(fit)
    cf <- sm$coefficients
    ci <- sm$conf.int
    data.frame(
      variable = v,
      level    = sub("^x", "", rownames(cf)),
      n        = nrow(d),
      HR       = unname(cf[, "exp(coef)"]),
      lower    = unname(ci[, "lower .95"]),
      upper    = unname(ci[, "upper .95"]),
      p        = unname(cf[, ncol(cf)]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(out)) {
    return(data.frame(variable = character(0), level = character(0), n = integer(0),
                      HR = numeric(0), lower = numeric(0), upper = numeric(0),
                      p = numeric(0), stringsAsFactors = FALSE))
  }
  out[order(out$p), , drop = FALSE]
}

# ---- plot -------------------------------------------------------------------

#' Tidy a survfit into a step-plottable data.frame (with the t=0 anchor)
#' @param fit A [km_fit()] result.
#' @keywords internal
km_tidy <- function(fit) {
  strata <- if (is.null(fit$strata)) rep("All", length(fit$time))
            else rep(sub("^\\.group=", "", names(fit$strata)), fit$strata)
  d <- data.frame(time = fit$time, surv = fit$surv, n_censor = fit$n.censor,
                  group = strata, stringsAsFactors = FALSE)
  # every curve starts at (0, 1), otherwise the step plot floats
  anchors <- do.call(rbind, lapply(unique(d$group), function(g)
    data.frame(time = 0, surv = 1, n_censor = 0, group = g,
               stringsAsFactors = FALSE)))
  d <- rbind(anchors, d)
  d[order(d$group, d$time), , drop = FALSE]
}

#' Kaplan-Meier plot (ggplot2, no survminer)
#'
#' @param fit A [km_fit()] result.
#' @param lr Optional [logrank_test()] result, annotated on the panel.
#' @param title Plot title.
#' @param time_label X axis label.
#' @return A ggplot object.
#' @keywords internal
km_plot <- function(fit, lr = NULL, title = "Overall survival",
                    time_label = "Months") {
  d <- km_tidy(fit)
  groups <- unique(d$group)
  cens <- d[d$n_censor > 0, , drop = FALSE]

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$time, y = .data$surv,
                                       colour = .data$group)) +
    ggplot2::geom_step(linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = sc_palette(length(groups)), name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = function(v) paste0(v * 100, "%")) +
    ggplot2::labs(x = time_label, y = "Survival probability", title = title) +
    scstudio_theme()

  if (nrow(cens)) {
    p <- p + ggplot2::geom_point(data = cens, shape = 3, size = 1.8,
                                 show.legend = FALSE)
  }
  if (length(groups) < 2) p <- p + ggplot2::theme(legend.position = "none")

  if (!is.null(lr) && is.finite(lr$p)) {
    lab <- if (lr$p < 0.0001) "log-rank p < 0.0001"
           else sprintf("log-rank p = %.4g", lr$p)
    p <- p + ggplot2::annotate("text", x = 0, y = 0.04, hjust = 0,
                               label = lab, size = 4, colour = "#8b98a5")
  }
  p
}
