# Accessors and S3 methods for joinpoint fits.

# Picks the model for a cohort: "final" or a number of joinpoints.
jp_pick <- function(f, model) {
  if (identical(model, "final")) model <- f$k_selected
  m <- f$models[[as.character(model)]]
  if (is.null(m)) {
    stop("No fitted model with ", model, " joinpoints.", call. = FALSE)
  }
  m
}

# Row binds a per cohort table with the cohort keys.
jp_collect <- function(object, fun, model = "final") {
  out <- lapply(seq_along(object$fits), function(g) {
    f <- object$fits[[g]]
    if (is.null(f)) {
      return(NULL)
    }
    tab <- fun(f, if (identical(model, "all")) NULL else jp_pick(f, model))
    if (is.null(tab) || !nrow(tab)) {
      return(NULL)
    }
    key <- object$cohorts[rep(g, nrow(tab)), , drop = FALSE]
    rownames(key) <- NULL
    cbind(key, tab)
  })
  out <- do.call(rbind, out)
  if (!is.null(out)) rownames(out) <- NULL
  out
}

jp_model_cols <- function(f, m) {
  data.frame(
    model = m$k, final = m$k == f$k_selected,
    ci_method = m$ci_method
  )
}

#' Trend estimates from a joinpoint fit
#'
#' `apc()` returns the annual percent change of each segment (the slope of
#' each segment for linear models), `aapc()` the average annual percent
#' change over the requested ranges, and `joinpoints()` the estimated
#' joinpoint locations with confidence intervals.
#'
#' For linear models (`log = FALSE`) `aapc()` returns the average annual
#' change (AAC), the length weighted mean of the segment slopes over the
#' range, with the same interval construction as the AAPC.
#'
#' @param object A `joinpoint` object.
#' @param model `"final"` for the final selected model, a number of
#'   joinpoints, or `"all"` for every fitted model.
#' @param ranges Optional list of `c(start, end)` ranges; recomputes the
#'   AAPC with parametric intervals for these ranges, clipped to the included
#'   data.
#' @return A data frame with one row per segment, range or joinpoint and the
#'   cohort keys.
#' @export
#' @examples
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, year, rate, kmax = 3, ci_method = "parametric")
#' apc(fit)
#' aapc(fit)
#' joinpoints(fit)
apc <- function(object, model = "final") {
  stopifnot(inherits(object, "joinpoint"))
  jp_collect(object, function(f, m) {
    ms <- if (is.null(m)) f$models else list(m)
    do.call(rbind, lapply(ms, function(mm) cbind(jp_model_cols(f, mm), mm$apc)))
  }, model)
}

#' @rdname apc
#' @export
aapc <- function(object, model = "final", ranges = NULL) {
  stopifnot(inherits(object, "joinpoint"))
  ctrl <- object$settings$control
  if (!is.null(ranges)) {
    if (!is.list(ranges)) ranges <- list(ranges)
    ctrl <- do.call(jp_control, utils::modifyList(unclass(ctrl),
                                                  list(aapc_ranges = ranges)))
  }
  jp_collect(object, function(f, m) {
    ms <- if (is.null(m)) f$models else list(m)
    rr <- if (!is.null(ranges)) {
      jp_aapc_ranges(f$data$x, f$data$included, FALSE, NULL, ranges)
    }
    do.call(rbind, lapply(ms, function(mm) {
      if (!is.null(ranges)) mm$ci_method <- "parametric"
      tab <- if (is.null(ranges)) {
        mm$aapc
      } else {
        jp_aapc_table(mm, rr, ctrl$alpha_aapc, "parametric")
      }
      if (is.null(tab)) NULL else cbind(jp_model_cols(f, mm), tab)
    }))
  }, model)
}

#' @rdname apc
#' @export
joinpoints <- function(object, model = "final") {
  stopifnot(inherits(object, "joinpoint"))
  jp_collect(object, function(f, m) {
    ms <- if (is.null(m)) f$models else list(m)
    do.call(rbind, lapply(ms, function(mm) {
      if (is.null(mm$joinpoints)) NULL else {
        cbind(jp_model_cols(f, mm), mm$joinpoints)
      }
    }))
  }, model)
}

#' Model selection results
#'
#' Returns, for each cohort, the selection statistics of all fitted models
#' (SSE, BIC, BIC3 and, depending on the method, the WBIC weight and value
#' or the modified BIC; for Poisson and binomial fits the deviance and
#' dispersion, the criteria using the deviance divided by the dispersion),
#' R squared and the Durbin-Watson statistic (undefined for constant series
#' and exact fits), the permutation or likelihood ratio tests performed, the
#' data dependent selection details and the MADWD steps.
#'
#' @param object A `joinpoint` object.
#' @return A list with data frames `criteria`, `tests`, `dds` and `madwd`
#'   (the last three may be `NULL`), and `selected` with the selected number
#'   of joinpoints per cohort.
#' @examples
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
#' model_selection(fit)$criteria
#' @export
model_selection <- function(object) {
  stopifnot(inherits(object, "joinpoint"))
  crit <- jp_collect(object, function(f, m) {
    tab <- f$criteria
    tab$df <- vapply(f$models, function(mm) {
      mm$n - 2 * (mm$k + 1) - length(mm$extra_names)
    }, numeric(1))
    if ("sse" %in% names(tab)) {
      tab$mse <- tab$sse / tab$df
    } else {
      tab$dispersion <- vapply(f$models, `[[`, numeric(1), "mse")
    }
    # Weighted R^2 and the Durbin-Watson statistic of the standardized
    # residuals (values well below 2 suggest positive autocorrelation); both
    # are undefined for a constant series or an exact fit.
    U <- jp_whiten_mat(f$S)
    tol <- jp_tiny(U, f$y)
    sst <- sum(qr.resid(qr(jp_apply_whiten(U, rep(1, length(f$y)))),
                        drop(jp_apply_whiten(U, f$y)))^2)
    tab$r_squared <- vapply(f$models, function(mm) {
      if (is.finite(mm$df) && sst > tol) 1 - mm$sse / sst else NA_real_
    }, numeric(1))
    tab$durbin_watson <- vapply(f$models, function(mm) {
      e2 <- sum(mm$wresid^2)
      if (e2 > tol) sum(diff(mm$wresid)^2) / e2 else NA_real_
    }, numeric(1))
    tab$selected <- tab$joinpoints == f$k_selected
    tab
  })
  tests <- jp_collect(object, function(f, m) f$permutation$table)
  dds <- jp_collect(object, function(f, m) {
    if (is.null(f$dds)) NULL else as.data.frame(f$dds)
  })
  madwd <- jp_collect(object, function(f, m) f$madwd$table)
  selected <- jp_collect(object, function(f, m) {
    data.frame(joinpoints = f$k_selected, statistical = f$k_statistical,
               method = f$method, kmin = f$kmin, kmax = f$kmax, rho = f$rho)
  })
  list(criteria = crit, tests = tests, dds = dds, madwd = madwd,
       selected = selected)
}

#' Model estimates
#'
#' Regression coefficients of the selected (or a specified) model in the
#' general changepoint parameterization: segment intercepts, segment slopes,
#' slope changes and any jump and covariate coefficients, with standard
#' errors from the unconstrained fit (Wald standard errors for Poisson and
#' binomial fits), test statistics and p values.
#'
#' @inheritParams apc
#' @return A data frame.
#' @examples
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
#' model_estimates(fit)
#' @export
model_estimates <- function(object, model = "final") {
  stopifnot(inherits(object, "joinpoint"))
  jp_collect(object, function(f, m) {
    ms <- if (is.null(m)) f$models else list(m)
    do.call(rbind, lapply(ms, function(mm) {
      k <- mm$k
      V <- mm$vcov_gp
      si <- 2L * seq_len(k + 1L)
      est <- c(mm$intercepts, mm$slopes, if (k) diff(mm$slopes))
      se <- c(mm$se_intercepts, mm$se_slopes, if (k) {
        vapply(seq_len(k), function(j) {
          sqrt(V[si[j], si[j]] + V[si[j + 1L], si[j + 1L]] -
                 2 * V[si[j], si[j + 1L]])
        }, numeric(1))
      })
      lab <- c(paste("Intercept", seq_len(k + 1L)),
               paste("Slope", seq_len(k + 1L)),
               if (k) paste0("Slope ", 2:(k + 1L), " - Slope ", 1:k))
      # Jump and covariate coefficients follow the segment parameters.
      en <- mm$extra_names
      lab <- c(lab, ifelse(en == "jump", "Jump", en))
      est <- c(est, mm$coef[en])
      se <- c(se, sqrt(diag(V))[2L * (k + 1L) + seq_along(en)])
      tt <- est / se
      out <- data.frame(
        parameter = lab, estimate = est, se = se, statistic = tt,
        p_value = if (mm$df > 0) 2 * stats::pt(-abs(tt), mm$df) else NA_real_,
        df = mm$df
      )
      cbind(jp_model_cols(f, mm), out)
    }))
  }, model)
}

#' Observed and modeled values
#'
#' Observed values, standard errors, the supplied counts and populations,
#' modeled values (the fitted trend, with covariates at their average) for
#' the selected model (or all models) and flags for excluded points, zero
#' counts replaced by 0.5 in the observed value and joinpoints.
#'
#' @param x A `joinpoint` object.
#' @param row.names,optional Ignored.
#' @param model `"final"`, a number of joinpoints, or `"all"`.
#' @param ... Ignored.
#' @return A data frame.
#' @examples
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
#' head(as.data.frame(fit))
#' @export
as.data.frame.joinpoint <- function(x, row.names = NULL, optional = FALSE,
                                    model = "final", ...) {
  object <- x
  jp_collect(object, function(f, m) {
    ms <- if (is.null(m)) f$models else list(m)
    do.call(rbind, lapply(ms, function(mm) {
      d <- f$data
      fit <- f$fitted[[as.character(mm$k)]]
      seg <- findInterval(d$x, mm$tau, left.open = TRUE) + 1L
      data.frame(
        model = mm$k, final = mm$k == f$k_selected, x = d$x, observed = d$r,
        se = d$se, count = d$count, population = d$population,
        modeled = fit, included = d$included,
        zero_count = d$corrected, joinpoint = d$x %in% mm$tau, segment = seg,
        apc = mm$apc$estimate[pmin(seg, nrow(mm$apc))]
      )
    }))
  }, model)
}

#' @export
fitted.joinpoint <- function(object, model = "final", ...) {
  stats::predict(object, model = model)$fit
}

#' @export
residuals.joinpoint <- function(object, model = "final",
                                type = c("response", "model", "weighted"),
                                ...) {
  type <- match.arg(type)
  unlist(lapply(object$fits, function(f) {
    if (is.null(f)) {
      return(NULL)
    }
    m <- jp_pick(f, model)
    switch(type,
      response = {
        inc <- f$data$included
        (f$data$r - jp_fitted_at(f, m, f$data$x, f$covariates))[inc]
      },
      model = m$resid,
      weighted = m$wresid
    )
  }), use.names = FALSE)
}

# Applies `fun` to the chosen model of each fitted cohort; one cohort gives
# the value itself, several a list named by cohort.
jp_by_cohort <- function(object, model, fun) {
  ok <- which(!vapply(object$fits, is.null, logical(1)))
  out <- lapply(ok, function(g) fun(jp_pick(object$fits[[g]], model)))
  names(out) <- vapply(ok, function(g) jp_cohort_label(object, g), "")
  if (length(out) == 1L) out[[1L]] else out
}

#' @export
coef.joinpoint <- function(object, model = "final", ...) {
  jp_by_cohort(object, model, function(m) m$coef)
}

#' Covariance of the segment parameters
#'
#' Covariance matrix of the segment intercepts and slopes and of any jump and
#' covariate coefficients, from the unconstrained segment fit used for the
#' standard errors (the Wald covariance for Poisson and binomial fits).
#'
#' @param object A `joinpoint` object.
#' @param model `"final"` or a number of joinpoints.
#' @param ... Ignored.
#' @return A matrix, or a list of matrices named by cohort.
#' @examples
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
#' vcov(fit)
#' @export
vcov.joinpoint <- function(object, model = "final", ...) {
  jp_by_cohort(object, model, jp_vcov)
}

jp_vcov <- function(m) {
  V <- m$vcov_gp
  nseg <- m$k + 1L
  nm <- as.vector(rbind(paste("Intercept", seq_len(nseg)),
                        paste("Slope", seq_len(nseg))))
  dimnames(V) <- rep(list(c(nm, m$extra_names)), 2L)
  V
}

#' @export
nobs.joinpoint <- function(object, ...) {
  sum(vapply(object$fits, function(f) {
    if (is.null(f)) 0L else length(f$x)
  }, integer(1)))
}

#' Predictions from a joinpoint fit
#'
#' Predictions include the effects of the covariates at their values in the
#' data (or in `newdata`). The fitted trend at the average of the
#' covariates, which the plots show, is the `modeled` column of
#' [as.data.frame.joinpoint()].
#'
#' @param object A `joinpoint` object.
#' @param newdata Optional data frame with the independent variable (and
#'   by-variables for multiple cohorts, and covariates of a fit with
#'   covariates). Defaults to the fitted data.
#' @param model `"final"` or a number of joinpoints.
#' @param scale `"response"` (original scale) or `"link"` (model scale).
#' @param ... Ignored.
#' @return A data frame with cohort keys, `x` and `fit`.
#' @examples
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
#' predict(fit, newdata = data.frame(year = c(1950, 1990)))
#' @export
predict.joinpoint <- function(object, newdata = NULL, model = "final",
                              scale = c("response", "link"), ...) {
  scale <- match.arg(scale)
  xn <- object$settings$x
  out <- lapply(seq_along(object$fits), function(g) {
    f <- object$fits[[g]]
    if (is.null(f)) {
      return(NULL)
    }
    m <- jp_pick(f, model)
    cv <- object$settings$covariates
    Z <- f$covariates
    xv <- if (is.null(newdata)) {
      f$data$x
    } else {
      nd <- newdata
      need <- c(xn, object$by, if (!is.null(cv)) all.vars(cv))
      if (!all(need %in% names(nd))) {
        stop("`newdata` needs the columns ", paste(need, collapse = ", "),
             ".", call. = FALSE)
      }
      if (!is.null(object$by)) {
        key <- object$cohorts[g, , drop = FALSE]
        keep <- Reduce(`&`, lapply(names(key), function(v) {
          nd[[v]] == key[[v]]
        }))
        nd <- nd[keep, , drop = FALSE]
      }
      if (!is.null(cv)) {
        ci <- f$covariate_info
        mf <- stats::model.frame(ci$terms, nd, xlev = ci$xlev,
                                 na.action = stats::na.pass)
        Z <- stats::model.matrix(ci$terms, mf, contrasts.arg = ci$contrasts)
        Z <- sweep(Z[, names(ci$center), drop = FALSE], 2L, ci$center)
      }
      nd[[xn]] + object$settings$shift
    }
    if (!length(xv)) {
      return(NULL)
    }
    fit <- if (scale == "response") {
      jp_fitted_at(f, m, xv, Z)
    } else {
      jp_linear(m, xv, Z)
    }
    key <- object$cohorts[rep(g, length(xv)), , drop = FALSE]
    rownames(key) <- NULL
    cbind(key, data.frame(x = xv, fit = fit))
  })
  do.call(rbind, out)
}

jp_cohort_label <- function(object, g) {
  if (is.null(object$by)) {
    return("All")
  }
  paste(vapply(object$cohorts[g, , drop = FALSE], as.character, character(1)),
        collapse = ", ")
}

# Labels of the trend statistics: percent changes for log-linear models,
# absolute changes (slopes) for linear models.
jp_stat_names <- function(interval, log = TRUE, family = "gaussian") {
  ab <- toupper(substr(interval, 1L, 1L))
  if (family == "binomial") {
    list(apc = paste(interval, "Percent Change in the Odds"),
         apc_abbr = paste0(ab, "PC (odds)"),
         aapc = paste("Average", interval, "Percent Change in the Odds"),
         aapc_abbr = paste0("A", ab, "PC (odds)"))
  } else if (log) {
    list(apc = paste(interval, "Percent Change"), apc_abbr = paste0(ab, "PC"),
         aapc = paste("Average", interval, "Percent Change"),
         aapc_abbr = paste0("A", ab, "PC"))
  } else {
    list(apc = paste(interval, "Change"), apc_abbr = "Slope",
         aapc = paste("Average", interval, "Change"),
         aapc_abbr = paste0("A", ab, "C"))
  }
}

#' @export
print.joinpoint <- function(x, digits = 3, ...) {
  s <- x$settings
  nm <- jp_stat_names(s$interval, s$log, s$family)
  cat("Joinpoint regression (", switch(s$family,
    gaussian = paste0(if (s$log) "log-linear" else "linear", " model, ",
                      s$variance, " variance"),
    poisson = "Poisson regression, log link",
    binomial = "binomial regression, logit link"), ")\n", sep = "")
  cat("Model selection: ", toupper(s$method), "; confidence intervals: ",
      s$ci_method, "\n", sep = "")
  for (g in seq_along(x$fits)) {
    f <- x$fits[[g]]
    cat("\nCohort: ", jp_cohort_label(x, g), "\n", sep = "")
    if (is.null(f)) {
      cat("  not fitted\n")
      next
    }
    m <- jp_pick(f, "final")
    cat("  Joinpoints: ", m$k, if (m$k) {
      paste0(" (", paste(format(m$tau), collapse = ", "), ")")
    }, "\n", sep = "")
    a <- m$apc
    lab <- nm$apc_abbr
    for (i in seq_len(nrow(a))) {
      cat(sprintf("  %s-%s  %s = %s%s  [%s, %s]\n", format(a$start[i]),
                  format(a$end[i]), lab,
                  jp_fmt(a$estimate[i], digits),
                  if (a$significant[i]) "*" else " ",
                  jp_fmt(a$lower[i], digits),
                  jp_fmt(a$upper[i], digits)))
    }
    if (!is.null(m$aapc)) {
      for (i in seq_len(nrow(m$aapc))) {
        r <- m$aapc[i, ]
        cat(sprintf("  %s %s = %s%s  [%s, %s]\n", nm$aapc_abbr, r$range,
                    jp_fmt(r$estimate, digits),
                    if (r$significant) "*" else " ",
                    jp_fmt(r$lower, digits),
                    jp_fmt(r$upper, digits)))
      }
    }
    for (nt in f$notes) cat("  Note: ", nt, "\n", sep = "")
  }
  cat("\n* Significantly different from zero at alpha = ",
      s$control$alpha_apc, ".\n", sep = "")
  invisible(x)
}

#' @export
summary.joinpoint <- function(object, ...) {
  structure(list(
    object = object, apc = apc(object), aapc = aapc(object),
    joinpoints = joinpoints(object), selection = model_selection(object),
    estimates = model_estimates(object)
  ), class = "summary.joinpoint")
}

#' @export
print.summary.joinpoint <- function(x, ...) {
  print(x$object)
  cat("\nModel selection:\n")
  print(x$selection$criteria, digits = 4, row.names = FALSE)
  if (!is.null(x$selection$tests)) {
    cat(if (x$object$settings$method == "lrt") "\nMonte Carlo likelihood ratio tests:\n"
        else "\nPermutation tests:\n")
    print(x$selection$tests, digits = 4, row.names = FALSE)
  }
  cat("\nModel estimates:\n")
  print(x$estimates, digits = 4, row.names = FALSE)
  if (!is.null(x$joinpoints)) {
    cat("\nJoinpoints:\n")
    print(x$joinpoints, digits = 6, row.names = FALSE)
  }
  invisible(x)
}

#' Publication table of trends
#'
#' Formats the APCs of the segments (slopes for linear models) and the
#' AAPCs (AACs) of the final, or a given, model as a character table with
#' estimates, confidence intervals and p values, ready for
#' `knitr::kable()` or other table packages.
#'
#' @inheritParams apc
#' @param digits Decimals of estimates and confidence limits.
#' @param aapc Append the AAPC rows.
#' @return A data frame of character columns.
#' @export
#' @examples
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, year, rate, kmax = 3, ci_method = "parametric")
#' jp_table(fit)
jp_table <- function(object, model = "final", digits = 2, aapc = TRUE) {
  stopifnot(inherits(object, "joinpoint"))
  nm <- jp_stat_names(object$settings$interval, object$settings$log,
                      object$settings$family)
  pv <- function(p) ifelse(is.na(p), "", sub("^= ", "", jp_pformat(p)))
  est <- function(d) {
    paste0(jp_fmt(d$estimate, digits), ifelse(d$significant, "*", ""), " (",
           jp_fmt(d$lower, digits), ", ", jp_fmt(d$upper, digits), ")")
  }
  keys <- object$by
  a <- apc(object, model)
  out <- data.frame(
    a[keys], Trend = paste("Segment", a$segment),
    Period = paste0(format(a$start), "-", format(a$end)),
    Estimate = est(a), p = pv(a$p_value), check.names = FALSE
  )
  if (aapc) {
    b <- aapc(object, model)
    if (!is.null(b)) {
      out <- rbind(out, data.frame(
        b[keys], Trend = nm$aapc_abbr,
        Period = paste0(format(b$start), "-", format(b$end)),
        Estimate = est(b), p = pv(b$p_value), check.names = FALSE
      ))
    }
  }
  ctl <- object$settings$control
  ci <- function(a) paste0(" (", 100 * (1 - a), "% CI)")
  first <- if (object$settings$log) nm$apc_abbr else "Slope"
  names(out)[names(out) == "Estimate"] <- if (!aapc) {
    paste0(first, ci(ctl$alpha_apc))
  } else if (ctl$alpha_apc == ctl$alpha_aapc) {
    paste0(first, " / ", nm$aapc_abbr, ci(ctl$alpha_apc))
  } else {
    paste0(first, ci(ctl$alpha_apc), " / ", nm$aapc_abbr, ci(ctl$alpha_aapc))
  }
  if (length(keys)) out <- out[do.call(order, unname(out[keys])), ]
  rownames(out) <- NULL
  out
}

#' Jump and comparability ratio estimates
#'
#' For jump models (see the `jump` argument of [joinpoint()]), the estimated
#' jump on the log scale and the comparability ratio `exp(jump)` with
#' standard errors, confidence intervals and tests (null values 0 and 1).
#'
#' @inheritParams apc
#' @return A data frame, or `NULL` when no jump model was fitted.
#' @examples
#' heart <- subset(us_mortality, cause == "Heart Disease" & year >= 1950)
#' fit <- joinpoint(heart, year, rate, kmax = 2, jump = 1978,
#'                  ci_method = "parametric")
#' jump_estimates(fit)
#' @export
jump_estimates <- function(object, model = "final") {
  stopifnot(inherits(object, "joinpoint"))
  jp_collect(object, function(f, m) {
    ms <- if (is.null(m)) f$models else list(m)
    do.call(rbind, lapply(ms, function(mm) {
      if (is.null(mm$jump_table)) NULL else {
        mm$ci_method <- "parametric"
        cbind(jp_model_cols(f, mm), mm$jump_table)
      }
    }))
  }, model)
}

# p values with three decimals, never shown as zero.
jp_pformat <- function(p) {
  ifelse(p < 0.001, "< 0.001", paste("=", formatC(p, digits = 3, format = "f")))
}
