# APC, AAPC and joinpoint location estimates with confidence intervals.

jp_pct <- function(b, log) if (log) 100 * expm1(b) else b

# Segment overlap weights of the range [a, b] for the segment bounds.
jp_range_weights <- function(bounds, a, b) {
  lo <- utils::head(bounds, -1L)
  hi <- utils::tail(bounds, -1L)
  len <- pmax(0, pmin(b, hi) - pmax(a, lo))
  len / (b - a)
}

# APC (log models) or slope (linear models) table for one model.
jp_apc_table <- function(model, alpha, ci_method, eq = NULL) {
  k <- model$k
  log <- model$log
  b <- model$slopes
  se <- model$se_slopes
  df <- model$df
  out <- data.frame(
    segment = seq_len(k + 1L),
    start = utils::head(model$bounds, -1L),
    end = utils::tail(model$bounds, -1L),
    estimate = jp_pct(b, log)
  )
  if (ci_method == "empirical" && !is.null(eq)) {
    q <- apply(eq$slopes, 2L, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2),
               na.rm = TRUE, names = FALSE)
    out$lower <- jp_pct(q[1L, ], log)
    out$upper <- jp_pct(q[2L, ], log)
    out$statistic <- NA_real_
    out$p_value <- apply(eq$slopes, 2L, jp_eq_pvalue)
  } else {
    tq <- if (!is.na(df) && df > 0) stats::qt(1 - alpha / 2, df) else NA
    out$lower <- jp_pct(b - tq * se, log)
    out$upper <- jp_pct(b + tq * se, log)
    out$statistic <- b / se
    out$p_value <- if (!is.na(df) && df > 0) {
      2 * stats::pt(-abs(b / se), df)
    } else {
      NA_real_
    }
  }
  out$significant <- !is.na(out$p_value) & out$p_value < alpha
  out
}

# AAPC table for the requested ranges (list of c(start, end) with labels).
jp_aapc_table <- function(model, ranges, alpha, ci_method, eq = NULL) {
  if (!length(ranges)) {
    return(NULL)
  }
  log <- model$log
  b <- model$slopes
  V <- model$vcov_slopes
  rows <- lapply(seq_along(ranges), function(i) {
    a <- ranges[[i]][1L]
    e <- ranges[[i]][2L]
    w <- jp_range_weights(model$bounds, a, e)
    mu <- sum(w * b)
    single <- sum(w > 1e-12) == 1L
    row <- data.frame(
      range = names(ranges)[i], start = a, end = e,
      estimate = jp_pct(mu, log), within_segment = single
    )
    if (ci_method == "empirical" && !is.null(eq)) {
      th <- vapply(seq_len(nrow(eq$slopes)), function(r) {
        bd <- c(model$bounds[1L], eq$tau[r, ], utils::tail(model$bounds, 1L))
        sum(jp_range_weights(bd, a, e) * eq$slopes[r, ])
      }, numeric(1))
      q <- stats::quantile(th, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE,
                           names = FALSE)
      row$lower <- jp_pct(q[1L], log)
      row$upper <- jp_pct(q[2L], log)
      row$statistic <- NA_real_
      row$p_value <- jp_eq_pvalue(th)
      row$se_log <- stats::sd(th, na.rm = TRUE)
    } else {
      # Only segments inside the range enter the variance.
      nz <- w > 1e-12
      se <- sqrt(drop(t(w[nz]) %*% V[nz, nz, drop = FALSE] %*% w[nz]))
      if (single) {
        df <- model$df
        crit <- if (df > 0) stats::qt(1 - alpha / 2, df) else NA_real_
        p <- if (df > 0) 2 * stats::pt(-abs(mu / se), df) else NA_real_
      } else {
        crit <- stats::qnorm(1 - alpha / 2)
        p <- 2 * stats::pnorm(-abs(mu / se))
      }
      row$lower <- jp_pct(mu - crit * se, log)
      row$upper <- jp_pct(mu + crit * se, log)
      row$statistic <- mu / se
      row$p_value <- p
      row$se_log <- se
    }
    row$log_estimate <- mu
    row
  })
  out <- do.call(rbind, rows)
  out$significant <- !is.na(out$p_value) & out$p_value < alpha
  out
}

# Joinpoint locations with confidence intervals. The parametric interval
# (Lerman 1980, equation 6) contains the candidates x with
# Q(x, j, k) <= Q_k (1 + k F / (n - p)), F = F(1 - alpha; k, n - p), where
# Q(x, j, k) is the minimum SSE with the j-th joinpoint at x, Q_k its minimum
# and p = 2k + 2 plus any jump and covariate coefficients; it is reported as
# the range of the accepted candidates and the estimate; for continuous
# joinpoints the limits are then located inside the gaps
# (jp_continuous_limits). Fits without a maximum likelihood estimate get
# no interval. The "literal"
# variant Q_k (1 + k / (n - p)) F is far more conservative. For Poisson and
# binomial fits Q is the profile deviance, and candidates within
# phi * chi-square(1 - alpha; k) of its minimum are accepted.
jp_tau_table <- function(model, alpha, ci_method, eq = NULL, cand = NULL,
                         form = "lerman") {
  k <- model$k
  if (k == 0) {
    return(NULL)
  }
  out <- data.frame(joinpoint = seq_len(k), estimate = model$tau,
                    lower = NA_real_, upper = NA_real_)
  if (ci_method == "empirical" && !is.null(eq)) {
    q <- apply(eq$tau, 2L, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2),
               na.rm = TRUE, names = FALSE)
    out$lower <- q[1L, ]
    out$upper <- q[2L, ]
  } else if (!is.null(model$profile) && !isFALSE(model$converged)) {
    crit <- jp_tau_crit(model, alpha, form)
    if (is.finite(crit)) {
      at <- model$profile_at %||% matrix(cand, length(cand), k)
      for (j in seq_len(k)) {
        ok <- at[model$profile[, j] <= crit * (1 + 1e-12), j]
        out$lower[j] <- min(ok, model$tau[j])
        out$upper[j] <- max(ok, model$tau[j])
      }
    }
  }
  out
}

# Cutoff of the profile SSE (deviance) for the parametric joinpoint interval.
jp_tau_crit <- function(model, alpha, form = "lerman") {
  k <- model$k
  df <- model$n - 2 * (k + 1) - length(model$extra_names)
  q <- min(model$profile)
  if (df <= 0 || !is.finite(q)) {
    return(NA_real_)
  }
  fq <- stats::qf(1 - alpha, k, df)
  if (!is.finite(model$df)) {
    q + model$mse * stats::qchisq(1 - alpha, k)
  } else if (form == "literal") {
    q * (1 + k / df) * fq
  } else {
    q * (1 + k * fq / df)
  }
}

# Requested AAPC ranges: full range, last X observations and custom ranges.
jp_aapc_ranges <- function(x_all, included, full = TRUE, last = NULL,
                           ranges = NULL) {
  xs <- sort(x_all)
  inc <- included[order(x_all)]
  xin <- xs[inc]
  out <- list()
  if (isTRUE(full)) {
    out[[sprintf("Full Range (%s-%s)", format(min(xin)), format(max(xin)))]] <-
      c(min(xin), max(xin))
  }
  for (L in last) {
    end_pos <- max(which(inc))
    start_pos <- end_pos - L + 1L
    if (start_pos < 1L || xs[start_pos] < min(xin)) next
    lab <- sprintf("Last %d Obs (%s-%s)", as.integer(L), format(xs[start_pos]),
                   format(xs[end_pos]))
    out[[lab]] <- c(xs[start_pos], xs[end_pos])
  }
  # Custom ranges are clipped to the included data; a range without overlap
  # is dropped.
  for (r in ranges) {
    lo <- max(r[1L], min(xin))
    hi <- min(r[2L], max(xin))
    if (lo >= hi) {
      warning("Ignoring AAPC range outside the data: ",
              paste(r, collapse = "-"), call. = FALSE)
      next
    }
    out[[sprintf("%s-%s", format(lo), format(hi))]] <- c(lo, hi)
  }
  out
}
