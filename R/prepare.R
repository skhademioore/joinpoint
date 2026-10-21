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
  if (type %in% c("crude_rate", "proportion", "percent") && is.null(r)) {
    if (is.null(count) || is.null(population)) {
      stop("`count` and `population` are required to calculate a ", type,
           ".", call. = FALSE)
    }
    if (any(population <= 0, na.rm = TRUE)) {
      stop("`population` must be positive.", call. = FALSE)
    }
    if (any(count < 0, na.rm = TRUE)) {
      stop("`count` must be nonnegative.", call. = FALSE)
    }
    if (all(count == 0, na.rm = TRUE)) stop("All counts are zero.", call. = FALSE)
    if (variance == "glm" && (any(count != round(count), na.rm = TRUE) ||
        (type != "crude_rate" &&
         any(population != round(population), na.rm = TRUE)))) {
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
      if (is.null(se)) se <- rate_per * sqrt(cnt) / population
    } else {
      p <- cnt / population
      if (any(p > 1, na.rm = TRUE)) {
        stop("Proportions must lie in [0, 1].", call. = FALSE)
      }
      mult <- if (type == "percent") 100 else 1
      r <- mult * p
      if (is.null(se)) se <- mult * sqrt(p * (1 - p) / population)
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
  if (variance == "se") {
    if (is.null(se)) {
      stop("`variance = \"se\"` requires standard errors (`se`).",
           call. = FALSE)
    }
    if (any(!is.finite(se) | se <= 0)) {
      stop("Standard errors must be positive and finite.", call. = FALSE)
    }
  }
  if (anyNA(x) || anyNA(r)) {
    stop("Missing values in the independent or dependent variable.",
         call. = FALSE)
  }
  if (!all(is.finite(x)) || !all(is.finite(r))) {
    stop("Infinite values in the independent or dependent variable.",
         call. = FALSE)
  }
  upper <- switch(type, proportion = 1, percent = 100, Inf)
  if (type != "other" && any(r < 0 | r > upper)) {
    stop("Values of a ", type, " must lie in [0, ", upper, "].",
         call. = FALSE)
  }
  if (anyDuplicated(x)) {
    stop("Duplicated values of the independent variable within a cohort.",
         call. = FALSE)
  }
  glm <- variance == "glm"
  if (log && !glm && any(r <= 0)) {
    stop("The log transformation requires a positive response.",
         call. = FALSE)
  }
  y <- if (log && !glm) base::log(r) else r
  v <- if (variance %in% c("se", "poisson")) {
    if (log) (se / r)^2 else se^2
  } else {
    rep(NA_real_, n)
  }
  included <- !(x %in% exclude)
  if (sum(included) < 2) {
    stop("At least two observations must remain after exclusions.",
         call. = FALSE)
  }
  o <- order(x)
  data.frame(
    x = (x + shift)[o], r = r[o], se = if (is.null(se)) NA_real_ else se[o],
    y = y[o], v = v[o], included = included[o], corrected = corrected[o],
    row = o, count = if (is.null(count)) NA_real_ else count[o],
    population = if (is.null(population)) NA_real_ else population[o]
  )
}
