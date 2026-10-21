#' Directly age adjusted rates
#'
#' Computes directly age standardized rates and their standard errors from
#' age specific counts and populations. With delay factors the result is
#' the delay adjusted age adjusted rate.
#'
#' For age groups `a` with counts `C_a`, populations `P_a`, delay factors
#' `D_a` and standard population weights `w_a = S_a / sum(S)`, the rate is
#' `M * sum(w_a * D_a * C_a / P_a)` with variance
#' `M^2 * sum(w_a^2 * D_a^2 * C_a / P_a^2)`, assuming Poisson counts.
#'
#' @param data A data frame with one row per age group, time point and
#'   cohort.
#' @param x Column of the independent variable (time).
#' @param count,population Columns of counts and populations.
#' @param age Column of age groups. For clustering their order matters: the
#'   levels of a factor, numeric order, or the leading number of each label
#'   (for example \code{"0-4"}, \code{"5-9"}, \code{"85+"}).
#' @param standard Either a column of `data` holding the standard
#'   population of each row, or a numeric vector named by age group (for
#'   example from [std_populations]).
#' @param by Optional character vector of by-variables.
#' @param delay Optional column of delay factors (1 to 20).
#' @param rate_per Rate multiplier.
#' @return A data frame with the by-variables, `x`, `rate`, `se`, `count`
#'   and `population`, ready for [joinpoint()] with `y = rate, se = se`.
#' @export
#' @examples
#' us <- subset(std_populations, standard == "segi1960")
#' std <- stats::setNames(us$population, us$age_group)
#' d <- expand.grid(year = 2000:2004, age = names(std))
#' d$count <- 10 + seq_len(nrow(d)) %% 7
#' d$pop <- 1e5
#' age_adjust(d, year, count, pop, age, std)
age_adjust <- function(data, x, count, population, age, standard, by = NULL,
                       delay = NULL, rate_per = 100000) {
  xn <- jp_col(substitute(x), data, "x")
  cn <- jp_col(substitute(count), data, "count")
  pn <- jp_col(substitute(population), data, "population")
  an <- jp_col(substitute(age), data, "age")
  dn <- jp_col(substitute(delay), data, "delay")
  std <- jp_std_rows(substitute(standard), standard, data, an)
  cnt <- data[[cn]]
  pop <- data[[pn]]
  dl <- if (is.null(dn)) rep(1, nrow(data)) else data[[dn]]
  jp_check_age_data(data, xn, an, by, cnt, pop, std, dl, rate_per)
  if (is.numeric(standard) && !is.null(names(standard))) {
    extra <- setdiff(names(standard), as.character(data[[an]]))
    if (length(extra)) {
      message("Standardized to the ", length(standard) - length(extra),
              " age groups in the data (of ", length(standard),
              " in `standard`).")
    }
  }
  keys <- data[c(by, xn)]
  grp <- interaction(keys, drop = TRUE, lex.order = TRUE)
  out <- lapply(split(seq_len(nrow(data)), grp), function(i) {
    w <- std[i] / sum(std[i])
    data.frame(
      keys[i[1L], , drop = FALSE],
      rate = rate_per * sum(w * dl[i] * cnt[i] / pop[i]),
      se = rate_per * sqrt(sum(w^2 * dl[i]^2 * cnt[i] / pop[i]^2)),
      count = sum(cnt[i]), population = sum(pop[i])
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

# Standard population of each row: `standard` names a column of `data` or
# is a vector named by age group.
jp_std_rows <- function(expr, standard, data, age) {
  if (is.name(expr) && as.character(expr) %in% names(data)) {
    return(data[[as.character(expr)]])
  }
  if (is.character(standard) && length(standard) == 1L &&
      standard %in% names(data)) {
    return(data[[standard]])
  }
  if (is.null(names(standard))) {
    stop("`standard` must be a column of `data` or a named vector.",
         call. = FALSE)
  }
  s <- standard[as.character(data[[age]])]
  if (anyNA(s)) {
    stop("Age groups without a standard population: ",
         paste(unique(data[[age]][is.na(s)]), collapse = ", "), call. = FALSE)
  }
  unname(s)
}

# Age specific data must have complete keys, positive populations and
# standard populations, nonnegative counts, delay factors in [1, 20], one
# row per age group and time, the same age groups at every time and in every
# cohort, and a standard population that does not change over time.
jp_check_age_data <- function(data, xn, an, by, cnt, pop, std, dl, rate_per) {
  if (anyNA(data[c(by, xn, an)])) {
    stop("Missing time, age group or by-variable values.", call. = FALSE)
  }
  if (!(is.numeric(rate_per) && length(rate_per) == 1L && is.finite(rate_per) &&
        rate_per > 0)) {
    stop("`rate_per` must be a positive number.", call. = FALSE)
  }
  vals <- cbind(cnt, pop, std, dl)
  if (anyNA(vals) || any(!is.finite(vals)) || any(pop <= 0) || any(cnt < 0) ||
      any(std <= 0) || any(dl < 1 | dl > 20)) {
    stop("Counts, populations, standard populations or delay factors are ",
         "missing or invalid.", call. = FALSE)
  }
  cohort <- if (length(by)) {
    interaction(data[by], drop = TRUE)
  } else {
    factor(rep(1L, nrow(data)))
  }
  if (anyDuplicated(data.frame(cohort, data[[xn]], data[[an]]))) {
    stop("Duplicated age group and time combinations.", call. = FALSE)
  }
  ages <- sort(unique(as.character(data[[an]])))
  for (g in split(seq_len(nrow(data)), cohort)) {
    tab <- table(data[[xn]][g], factor(data[[an]][g], ages))
    if (any(tab == 0)) {
      stop("Every age group needs a row at every time in every cohort.",
           call. = FALSE)
    }
    if (any(tapply(std[g], as.character(data[[an]][g]),
                   function(v) diff(range(v))) > 0)) {
      stop("The standard population of an age group varies over time.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}
