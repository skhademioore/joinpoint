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
    sse[cols, ] <- r$sse[, seq_len(kmax + 1L), drop = FALSE]
    for (k in seq_len(kmax)) idx[[k + 1L]][cols, ] <- r$idx[[k + 1L]]
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
  tol <- jp_tiny(U, models[[1L]]$fitted + models[[1L]]$resid)
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
    }, n_perm, n_exceed)
    lev <- alpha / (kmax - a)
    reject <- mc[["p"]] < lev
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
    hit <- t >= stat * (1 - 1e-10)
    if (!is.null(h) && r + sum(hit) >= h) {
      L <- used + which(r + cumsum(hit) >= h)[1L]
      return(c(p = h / L, n = L))
    }
    r <- r + sum(hit)
    used <- used + length(t)
  }
  c(p = (r + 1) / (used + 1), n = used)
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
# fraction of the whitened sum of squares about the (generalized least
# squares) mean, or the rounding error of the response when that is zero.
# It does not change when a constant is added to the response.
jp_tiny <- function(U, y) {
  yw <- drop(jp_apply_whiten(U, y))
  dev <- qr.resid(qr(jp_apply_whiten(U, rep(1, length(y)))), yw)
  max(1e-10 * sum(dev^2), 100 * length(y) * (.Machine$double.eps * max(abs(yw)))^2,
      .Machine$double.xmin)
}

# Resampled residuals for the empirical quantile method (Kim et al. 2017,
# method 2): the inverse of the empirical distribution function of the
# residuals, linearly interpolated between the plotting positions i/(n + 1),
# is applied to uniform random numbers. The support is extended by
# z_0 = z_(1) - D and z_(n+1) = z_(n) + D with D = ln(3 + ln n) * IQR(z).
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
  matrix(lo + (hi - lo) * up, n, B)
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
    b[1L, ] <- b[1L, ] - b[2L, ] * ctr
    coefs[cols, ] <- t(b)
    slopes[cols, ] <- t(apply(b, 2L, jp_slopes, k = k))
  }
  list(slopes = slopes, coef = coefs, tau = taus)
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
# resampled estimates on either side of zero.
jp_eq_pvalue <- function(theta) {
  theta <- theta[is.finite(theta)]
  B <- length(theta)
  2 * min(sum(theta < 0), sum(theta >= 0)) / B
}
