# Input processing: response type, standard errors, weights, exclusions.

# Resolves a column given as a bare name or a string.
jp_col <- function(expr, data, arg) {
  if (is.null(expr)) {
    return(NULL)
  }
  if (is.character(expr)) {
    nm <- expr
  } else if (is.name(expr)) {
    nm <- as.character(expr)
  } else {
    stop("`", arg, "` must be a column name.", call. = FALSE)
  }
  if (!all(nm %in% names(data))) {
    stop("Column `", paste(setdiff(nm, names(data)), collapse = "`, `"),
         "` not found in `data`.", call. = FALSE)
  }
  nm
}

# Recommended maximum number of joinpoints for n observations: one more
# joinpoint for every five observations beyond six, at most seven.
jp_recommended_kmax <- function(n) {
  as.integer(pmax(0, pmin(7, floor((n - 2) / 5))))
}

# Builds the per cohort analysis frame. Returns a data.frame with x, the
# original response r, the model response y, model scale variance (NA when
# constant variance), inclusion flag and correction flag.
jp_prepare <- function(x, r, se, count, population, type, variance, log,
                       rate_per, exclude, shift) {
  n <- length(x)
  corrected <- rep(FALSE, n)
  if (anyNA(x) || !all(is.finite(x))) {
    stop("Missing or infinite values in the independent variable.",
         call. = FALSE)
  }
  # Exclusions are known from x alone; the checks of the response and its
  # transformation apply to the included observations only.
  included <- !(x %in% exclude)
  if (sum(included) < 2) {
    stop("At least two observations must remain after exclusions.",
         call. = FALSE)
  }
  if (variance == "poisson" && !is.null(se)) {
    stop("Poisson variance is calculated from the counts; do not supply ",
         "`se` with `variance = \"poisson\"`.", call. = FALSE)
  }
  if (type == "crude_rate" && !is.null(r) && variance == "poisson" &&
      !is.null(count)) {
    # A supplied crude rate with its counts: the Poisson standard error has
    # the relative size 1 / sqrt(count) whatever the rate multiplier.
    if (any(count[included] < 0, na.rm = TRUE)) {
      stop("`count` must be nonnegative.", call. = FALSE)
    }
    se <- r / sqrt(pmax(count, 0.5))
  } else if (type %in% c("crude_rate", "proportion", "percent") &&
             is.null(r)) {
    if (is.null(count) || is.null(population)) {
      stop("`count` and `population` are required to calculate a ", type,
           ".", call. = FALSE)
    }
    ci <- count[included]
    pi <- population[included]
    if (any(pi <= 0, na.rm = TRUE)) {
      stop("`population` must be positive.", call. = FALSE)
    }
    if (any(ci < 0, na.rm = TRUE)) {
      stop("`count` must be nonnegative.", call. = FALSE)
    }
    if (all(ci == 0, na.rm = TRUE)) stop("All counts are zero.", call. = FALSE)
    if (variance == "glm" && (any(ci != round(ci), na.rm = TRUE) ||
        (type != "crude_rate" && any(pi != round(pi), na.rm = TRUE)))) {
      stop("Poisson and binomial fits need whole-number counts and ",
           "binomial totals.", call. = FALSE)
    }
    cnt <- count
    if (type == "crude_rate" &&
        (variance == "poisson" || (variance == "constant" && log))) {
      corrected <- cnt == 0
      cnt[corrected] <- 0.5
    }
    if (type == "crude_rate") {
      r <- rate_per * cnt / population
      if (is.null(se)) se <- rate_per * sqrt(pmax(cnt, 0)) / population
    } else {
      p <- cnt / population
      if (any(p[included] > 1, na.rm = TRUE)) {
        stop("Proportions must lie in [0, 1].", call. = FALSE)
      }
      mult <- if (type == "percent") 100 else 1
      r <- mult * p
      if (is.null(se)) se <- mult * sqrt(pmax(p * (1 - p), 0) / population)
    }
  } else if (type == "count") {
    if (is.null(r)) r <- count
    if (all(r == 0, na.rm = TRUE)) stop("All counts are zero.", call. = FALSE)
    if (variance == "poisson" || (variance == "constant" && log)) {
      corrected <- r == 0
      r[corrected] <- 0.5
    }
    if (variance == "poisson") se <- sqrt(r)
  }
  if (is.null(r)) {
    stop("A response `y` (or `count` and `population`) is required.",
         call. = FALSE)
  }
  if (variance == "poisson") {
    if (!(type == "count" || (type == "crude_rate" && !is.null(count)))) {
      stop("Poisson variance requires a count response or a crude rate ",
           "calculated from `count` and `population`.", call. = FALSE)
    }
  }
  # Excluded rows are only shown: values that cannot be computed (such as a
  # rate with a zero population) are missing.
  r[!included & !is.finite(r)] <- NA
  if (!is.null(se)) se[!included & !is.finite(se)] <- NA
  ri <- r[included]
  if (variance == "se") {
    if (is.null(se)) {
      stop("`variance = \"se\"` requires standard errors (`se`).",
           call. = FALSE)
    }
    if (any(!is.finite(se[included]) | se[included] <= 0)) {
      stop("Standard errors must be positive and finite.", call. = FALSE)
    }
  }
  if (anyNA(ri)) {
    stop("Missing values in the dependent variable.", call. = FALSE)
  }
  if (!all(is.finite(ri))) {
    stop("Infinite values in the dependent variable.", call. = FALSE)
  }
  upper <- switch(type, proportion = 1, percent = 100, Inf)
  if (type != "other" && any(ri < 0 | ri > upper)) {
    stop("Values of a ", type, " must lie in [0, ", upper, "].",
         call. = FALSE)
  }
  if (anyDuplicated(x)) {
    stop("Duplicated values of the independent variable within a cohort.",
         call. = FALSE)
  }
  glm <- variance == "glm"
  if (log && !glm && any(ri <= 0)) {
    stop("The log transformation requires a positive response.",
         call. = FALSE)
  }
  # Excluded observations that cannot be transformed get a missing model
  # scale value; they are only shown.
  y <- if (log && !glm) {
    suppressWarnings(ifelse(!is.na(r) & r > 0, base::log(r), NA_real_))
  } else {
    r
  }
  v <- if (variance %in% c("se", "poisson")) {
    if (log) (se / r)^2 else se^2
  } else {
    rep(NA_real_, n)
  }
  v[!included] <- ifelse(is.finite(v[!included]), v[!included], NA_real_)
  o <- order(x)
  data.frame(
    x = (x + shift)[o], r = r[o], se = if (is.null(se)) NA_real_ else se[o],
    y = y[o], v = v[o], included = included[o], corrected = corrected[o],
    row = o, count = if (is.null(count)) NA_real_ else count[o],
    population = if (is.null(population)) NA_real_ else population[o]
  )
}
