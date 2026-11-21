# Permutation tests and empirical quantile resampling.

# Runs the compiled grid search over the columns of Y in chunks.
jp_grid_batch <- function(x, Y, S, cand, cons, kmax, extra = NULL,
                          chunk = 512L) {
  W <- jp_weights(S)
  Y <- as.matrix(Y)
  nd <- ncol(Y)
  sse <- matrix(NA_real_, nd, kmax + 1L)
  idx <- lapply(0:kmax, function(k) matrix(NA_integer_, nd, k))
  for (s in seq(1L, nd, by = chunk)) {
    cols <- s:min(nd, s + chunk - 1L)
    r <- jp_grid_search(x, Y[, cols, drop = FALSE], W, cand, cons, kmax,
                        extra)
    # The compiled search lowers kmax to the number of candidates.
    kk <- ncol(r$sse) - 1L
    sse[cols, seq_len(kk + 1L)] <- r$sse
    for (k in seq_len(kk)) idx[[k + 1L]][cols, ] <- r$idx[[k + 1L]]
  }
  list(sse = sse, idx = idx)
}

# Evaluates `expr` with a fixed random seed, restoring the global RNG state.
jp_with_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(expr)
  }
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  on.exit({
    if (is.null(old)) {
      rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", old, envir = globalenv())
    }
  })
  set.seed(seed)
  expr
}

# Sequential permutation test (Kim et al. 2000). Each test of H0: k = a
# against H1: k = b uses the significance level alpha / (MAX - a), starting
# with a = MIN and b = MAX.
jp_permutation_select <- function(x, S, models, cand, cons, kmin, kmax,
                                  alpha, n_perm, extra = NULL,
                                  n_exceed = NULL) {
  ks <- vapply(models, `[[`, numeric(1), "k")
  U <- jp_whiten_mat(S)
  tol <- jp_tiny(U, models[[1L]]$fitted + models[[1L]]$resid, x, extra)
  a <- kmin
  b <- kmax
  rows <- list()
  while (a < b) {
    m0 <- models[[which(ks == a)]]
    m1 <- models[[which(ks == b)]]
    stat <- jp_ratio(m0$sse, m1$sse, tol)
    ew <- drop(jp_apply_whiten(U, m0$resid))
    mc <- jp_mc_pvalue(stat, function(nb) {
      E <- vapply(seq_len(nb), function(i) ew[sample.int(length(ew))],
                  numeric(length(ew)))
      g <- jp_grid_batch(x, m0$fitted + jp_unwhiten(U, E), S, cand, cons, b,
                         extra)
      jp_ratio(g$sse[, a + 1L], g$sse[, b + 1L], tol)
    }, n_perm, n_exceed, batch = if (is.null(n_exceed)) 512L else 32L)
    lev <- alpha / (kmax - a)
    jp_mc_check(mc, a, b)
    reject <- isTRUE(mc[["p"]] < lev)
    rows[[length(rows) + 1L]] <- data.frame(
      null = a, alternative = b, statistic = stat, p_value = mc[["p"]],
      alpha = lev, reject = reject, simulations = mc[["n"]]
    )
    if (reject) a <- a + 1 else b <- b - 1
  }
  tab <- if (length(rows)) {
    do.call(rbind, rows)
  } else {
    data.frame(null = integer(0), alternative = integer(0),
               statistic = numeric(0), p_value = numeric(0),
               alpha = numeric(0), reject = logical(0),
               simulations = numeric(0))
  }
  list(k = a, table = tab)
}

# Monte Carlo p value from at most B simulated statistics, drawn in batches
# by `draw(nb)` (missing values are failed simulations and are left out):
# (r + 1) / (B + 1), where r statistics reach the observed one, or, with
# `h`, h / L as soon as the h-th does so at simulation L (Besag and
# Clifford 1991). Returns the p value and the number of simulations used.
jp_mc_pvalue <- function(stat, draw, B, h = NULL, batch = 512L) {
  tried <- 0L
  used <- 0L
  r <- 0L
  while (tried < B) {
    nb <- min(batch, B - tried)
    tried <- tried + nb
    t <- draw(nb)
    t <- t[!is.na(t)]
    hit <- jp_reaches(t, stat)
    if (!is.null(h) && r + sum(hit) >= h) {
      L <- used + which(r + cumsum(hit) >= h)[1L]
      return(c(p = h / L, n = L, failed = tried - used - length(t)))
    }
    r <- r + sum(hit)
    used <- used + length(t)
  }
  # Without a usable simulated statistic the test is undecided.
  c(p = if (used > 0L) (r + 1) / (used + 1) else NA_real_, n = used,
    failed = tried - used)
}

# Warns about simulated data sets that could not be fitted (left out of the
# Monte Carlo p value) and about tests left undecided.
jp_mc_check <- function(mc, a, b) {
  if (is.na(mc[["p"]])) {
    warning("No simulated data set of the test of ", a, " against ", b,
            " joinpoints could be fitted; the test is undecided and the ",
            "null model is kept.", call. = FALSE)
  } else if (mc[["failed"]] > 0) {
    warning(mc[["failed"]], " simulated data sets of the test of ", a,
            " against ", b, " joinpoints could not be fitted and were left ",
            "out of its p value.", call. = FALSE)
  }
}

# Simulated statistics that reach the observed one, ties included: the
# tolerance is relative to the observed statistic and never raises the
# threshold, whatever the sign of the statistic.
jp_reaches <- function(t, stat) {
  thr <- if (is.finite(stat)) stat - 1e-10 * abs(stat) else stat
  t >= thr
}

# SSE ratio of a null and an alternative model. SSEs below `tol` are exact
# fits: a perfect alternative fit gives an infinite ratio unless the null
# model is also perfect, which gives no evidence of a difference.
jp_ratio <- function(sse0, sse1, tol = 0) {
  sse0[sse0 < tol] <- 0
  sse1[sse1 < tol] <- 0
  ifelse(sse1 > 0, sse0 / sse1, ifelse(sse0 > 0, Inf, 1))
}

# Sums of squares below this tolerance are numerically zero: a small
# fraction of the whitened sum of squares left by the design common to all
# models, or the rounding error of the response when that is zero. It does
# not change when a component of the common design (a constant, a linear
# trend, a covariate or a jump) is added to the response.
jp_tiny <- function(U, y, x = NULL, extra = NULL, rel = 1e-10) {
  # The numerical zero is the rounding error of a residual vector of the
  # (whitened) responses. For simulated responses, whose sums of squares
  # come from the compiled search on cross products, it also includes the
  # fraction `rel` of the variation left after the base model that every
  # candidate contains (intercept, slope given `x`, and the `extra` columns
  # of covariates and jumps), the precision of that search; the observed
  # responses are searched to full precision and use `rel = 0`.
  yw <- drop(jp_apply_whiten(U, y))
  B <- if (is.null(x) || !isTRUE(stats::sd(x) > 0)) {
    matrix(1, length(y), 1L)
  } else {
    cbind(1, (x - mean(x)) / stats::sd(x))
  }
  if (!is.null(extra)) B <- cbind(B, extra)
  dev <- qr.resid(qr(jp_apply_whiten(U, B)), yw)
  max(rel * sum(dev^2),
      100 * length(y) * (.Machine$double.eps * max(abs(yw)))^2,
      .Machine$double.xmin)
}

# Resampled residuals for the empirical quantile method (Kim et al. 2017,
# method 2): the inverse of the empirical distribution function of the
# residuals, linearly interpolated between the plotting positions i/(n + 1),
# is applied to uniform random numbers. The support is extended by
# z_0 = z_(1) - D and z_(n+1) = z_(n) + D with D = ln(3 + ln n) * IQR(z).
# The draws are centered at the mean of this distribution,
# (sum(z) + (z_(1) + z_(n)) / 2) / (n + 1): with weights or correlated
# errors the whitened residuals need not average zero, and an offset would
# shift the resampled trends once unwhitened.
jp_eq_residuals <- function(z, B) {
  n <- length(z)
  zs <- sort(z)
  D <- log(3 + log(n)) * stats::IQR(z)
  ext <- c(zs[1L] - D, zs, zs[n] + D)
  u <- matrix(stats::runif(n * B), n, B)
  i <- floor(u * (n + 1))
  up <- (u - i / (n + 1)) * (n + 1)
  lo <- ext[i + 1L]
  hi <- ext[i + 2L]
  centre <- (sum(zs) + (zs[1L] + zs[n]) / 2) / (n + 1)
  matrix(lo + (hi - lo) * up, n, B) - centre
}

# Refits the model for each column of Y given joinpoints per column and
# returns slopes and joinpoints for each resample.
jp_refit_many <- function(x, Y, S, taus, k, extra = NULL) {
  B <- ncol(Y)
  U <- jp_whiten_mat(S)
  Yw <- jp_apply_whiten(U, Y)
  ctr <- mean(x)
  slopes <- matrix(NA_real_, B, k + 1L)
  coefs <- matrix(NA_real_, B, 2L + k + if (is.null(extra)) 0L else ncol(extra))
  coefs_c <- coefs
  key <- if (k) apply(taus, 1L, paste, collapse = "|") else rep("", B)
  for (kk in unique(key)) {
    cols <- which(key == kk)
    tau <- if (k) taus[cols[1L], ] else numeric(0)
    if (anyNA(tau)) next
    X <- jp_design(x - ctr, tau - ctr, extra)
    Xw <- jp_apply_whiten(U, X)
    b <- qr.coef(qr(Xw), Yw[, cols, drop = FALSE])
    b[is.na(b)] <- 0
    b <- matrix(b, nrow = ncol(X))
    coefs_c[cols, ] <- t(b)
    b[1L, ] <- b[1L, ] - b[2L, ] * ctr
    coefs[cols, ] <- t(b)
    slopes[cols, ] <- t(apply(b, 2L, jp_slopes, k = k))
  }
  list(slopes = slopes, coef = coefs, coef_c = coefs_c, ctr = ctr,
       tau = taus)
}

# Empirical quantile resampling for a fitted model.
jp_eq_resample <- function(x, S, model, cand, cons, n_res, extra = NULL) {
  U <- jp_whiten_mat(S)
  k <- model$k
  ew <- drop(jp_apply_whiten(U, model$resid))
  E <- jp_eq_residuals(ew, n_res)
  Ystar <- model$fitted + jp_unwhiten(U, E)
  if (k > 0) {
    g <- jp_grid_batch(x, Ystar, S, cand, cons, k, extra)
    taus <- matrix(cand[g$idx[[k + 1L]]], n_res, k)
  } else {
    taus <- matrix(numeric(0), n_res, 0)
  }
  jp_refit_many(x, Ystar, S, taus, k, extra)
}

# Empirical quantile p value for H0: theta = 0: twice the smaller share of
# resampled estimates at or beyond zero on either side, at most 1. Estimates
# within `tol` of zero (roundoff) are ties and count on both sides, so a
# distribution concentrated at zero gives p = 1.
jp_eq_pvalue <- function(theta, tol = 0) {
  theta <- theta[is.finite(theta)]
  B <- length(theta)
  if (!B) {
    return(NA_real_)
  }
  zero <- abs(theta) <= tol
  # Without any resample on the other side of zero the p value is below
  # the resolution 2 / B of the resamples; it is reported as 1 / B rather
  # than 0.
  max(1 / B, min(1, 2 * min(sum(theta < 0 | zero), sum(theta > 0 | zero)) / B))
}

# Rounding error bound of the least squares slopes of the model, one per
# segment (`segment = NULL`: over the whole range, for AAPCs). A slope is a
# linear combination sum(c_i y_i) of the responses of its segment, with
# c_i = (x_i - mean(x)) / sum((x - mean(x))^2), so storing and resampling
# the responses in double precision perturbs it by at most a small multiple
# of eps * sum(|c_i| |y_i|). Resampled estimates within this bound of zero
# are ties at zero (as for an exactly flat series); a slope that is large
# against the bound is resolved whatever the level (intercept) of the
# response.
jp_slope_tol <- function(model, x = NULL, segment = NULL) {
  y <- model$fitted + model$resid
  if (is.null(x)) {
    x <- seq(min(model$bounds), max(model$bounds), length.out = length(y))
  }
  bound <- function(i) {
    if (sum(i) < 2L) {
      return(0)
    }
    xc <- x[i] - mean(x[i])
    sxx <- sum(xc^2)
    if (!is.finite(sxx) || sxx <= 0) {
      return(0)
    }
    64 * .Machine$double.eps * sum(abs(xc) * abs(y[i])) / sxx
  }
  if (is.null(segment)) {
    return(bound(rep(TRUE, length(x))))
  }
  b <- model$bounds
  vapply(segment, function(j) {
    bound(x >= b[j] & x <= b[j + 1L])
  }, numeric(1))
}
