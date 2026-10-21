# Fitting of joinpoint models with given joinpoints and their parametric
# inference.

# Hinge design: intercept, x, (x - tau_j)_+, plus optional extra columns.
jp_design <- function(x, tau, extra = NULL) {
  X <- cbind(1, x, outer(x, tau, function(a, b) pmax(a - b, 0)))
  colnames(X) <- c("intercept", "slope", if (length(tau)) {
    paste0("delta", seq_along(tau))
  })
  if (!is.null(extra)) X <- cbind(X, extra)
  X
}

# Generalized least squares fit with fixed joinpoints. The x variable is
# centered internally and coefficients are returned on the original x scale.
jp_fit_fixed <- function(x, y, S, tau, extra = NULL, U = jp_whiten_mat(S)) {
  ctr <- mean(x)
  X <- jp_design(x - ctr, tau - ctr, extra)
  Xw <- jp_apply_whiten(U, X)
  yw <- drop(jp_apply_whiten(U, y))
  qx <- qr(Xw)
  b <- qr.coef(qx, yw)
  b[is.na(b)] <- 0
  fitted <- drop(X %*% b)
  rw <- drop(yw - Xw %*% b)
  b[1L] <- b[1L] - b[2L] * ctr
  names(b) <- colnames(X)
  list(
    coef = b, fitted = fitted, resid = y - fitted, wresid = rw,
    sse = sum(rw^2), tau = tau
  )
}

# Slopes of the k + 1 segments from the standard parameterization.
jp_slopes <- function(coef, k) {
  unname(cumsum(c(coef[[2L]], if (k) coef[2L + seq_len(k)])))
}

# Segment intercepts of the general changepoint parameterization.
jp_intercepts <- function(coef, tau) {
  k <- length(tau)
  coef[[1L]] - cumsum(c(0, if (k) coef[2L + seq_len(k)] * tau))
}

# Evaluates the fitted joinpoint function (model scale) at new x values.
jp_eval <- function(coef, tau, xnew, extra_new = NULL) {
  X <- jp_design(xnew, tau, extra_new)
  drop(X %*% coef[seq_len(ncol(X))])
}

# Standard errors of the general changepoint parameters (Lerman 1980;
# Feder 1975). Observations lying exactly on a
# joinpoint are deleted, the remaining data are partitioned by the estimated
# joinpoints, and an unconstrained least squares fit with a separate intercept
# and slope for each segment gives the covariance matrix. `jump` is the
# location of a coding change (jump model) or NULL.
jp_unconstrained <- function(x, y, S, tau, jump = NULL, covars = NULL) {
  k <- length(tau)
  # Joinpoints on observations are copies of observed x values.
  keep <- which(!(x %in% tau))
  xs <- x[keep]
  ys <- y[keep]
  Ss <- jp_subset_cov(S, keep)
  seg <- findInterval(xs, tau) + 1L
  nseg <- k + 1L
  ctr <- mean(xs)
  X <- matrix(0, length(xs), 2L * nseg)
  for (j in seq_len(nseg)) {
    X[seg == j, 2L * j - 1L] <- 1
    X[seg == j, 2L * j] <- xs[seg == j] - ctr
  }
  if (!is.null(jump)) {
    X <- cbind(X, as.numeric(seg == findInterval(jump, tau) + 1L & xs >= jump))
  }
  # Covariates enter every segment with common coefficients.
  if (!is.null(covars)) X <- cbind(X, covars[keep, , drop = FALSE])
  U <- jp_whiten_mat(Ss)
  Xw <- jp_apply_whiten(U, X)
  yw <- drop(jp_apply_whiten(U, ys))
  qx <- qr(Xw)
  p <- qx$rank
  d <- length(xs) - p
  rss <- sum(qr.resid(qx, yw)^2)
  sigma2 <- if (d > 0) rss / d else NA_real_
  V <- matrix(NA_real_, ncol(X), ncol(X))
  if (p > 0 && d > 0) {
    piv <- qx$pivot[seq_len(p)]
    R <- qr.R(qx)[seq_len(p), seq_len(p), drop = FALSE]
    Ri <- backsolve(R, diag(p))
    V[piv, piv] <- sigma2 * Ri %*% t(Ri)
  }
  # Segments with fewer than two retained observations have no estimable
  # slope; their statistics are reported as missing.
  n_seg <- tabulate(seg, nseg)
  bad <- which(n_seg < 2L)
  for (j in bad) {
    V[2L * j, ] <- V[, 2L * j] <- NA_real_
    if (n_seg[j] < 1L) V[2L * j - 1L, ] <- V[, 2L * j - 1L] <- NA_real_
  }
  # Intercepts on the original x scale: a_j = a_j(c) - b_j * ctr. Missing
  # entries are kept out of the product so they do not spread to other
  # segments.
  A <- diag(ncol(X))
  for (j in seq_len(nseg)) A[2L * j - 1L, 2L * j] <- -ctr
  miss <- is.na(diag(V))
  V[is.na(V)] <- 0
  Vo <- A %*% V %*% t(A)
  Vo[miss, ] <- Vo[, miss] <- NA_real_
  list(vcov = Vo, df = d, slope_index = 2L * seq_len(nseg),
       intercept_index = 2L * seq_len(nseg) - 1L)
}

# Assembles a fitted k joinpoint model with its parametric inference.
jp_build_model <- function(x, y, S, tau, log, jump = NULL, covars = NULL) {
  k <- length(tau)
  extra <- cbind(if (!is.null(jump)) cbind(jump = as.numeric(x >= jump)),
                 covars)
  ft <- jp_fit_fixed(x, y, S, tau, extra)
  n <- length(x)
  un <- jp_unconstrained(x, y, S, tau, jump, covars)
  slopes <- jp_slopes(ft$coef, k)
  ints <- jp_intercepts(ft$coef, tau)
  V <- un$vcov
  se_slope <- sqrt(diag(V)[un$slope_index])
  se_int <- sqrt(diag(V)[un$intercept_index])
  bounds <- c(min(x), tau, max(x))
  list(
    k = k, tau = tau, coef = ft$coef, slopes = slopes, intercepts = ints,
    se_slopes = se_slope, se_intercepts = se_int, vcov_gp = V,
    vcov_slopes = V[un$slope_index, un$slope_index, drop = FALSE],
    df = un$df, sse = ft$sse,
    mse = ft$sse / (n - 2 * (k + 1) - length(colnames(extra))),
    fitted = ft$fitted, resid = ft$resid,
    wresid = ft$wresid, n = n, bounds = bounds, jump = jump,
    extra_names = colnames(extra), log = log,
    link = if (log) "log" else "identity"
  )
}

# Continuous joinpoints (Hudson 1966) by exhaustive search. `cand` holds the
# observations and one point inside each gap between them. Inside the gap
# (x_a, x_b) the hinge satisfies (x - tau)_+ = (x - x_b)_+ + (x_b - tau) I(x >= x_b)
# at every observation, so the best joinpoint there follows from a linear
# fit with a hinge and a step at x_b: tau = x_b - step / hinge, admissible
# when it falls inside the gap (otherwise the best placement is at an
# observation, which is examined as well). Joinpoints may lie anywhere in
# the closure of the placement rule: an observation on a joinpoint counts
# toward one of its two segments, so the minimum exists. `cons` must be
# built with `closed = TRUE`. Returns the least squares joinpoints and the
# profile: for
# each joinpoint, the smallest SSE with the joinpoint at each observation or
# inside each gap, and where in the gap. With `fix = list(j, pos, tau)` it
# returns the smallest SSE with joinpoint j fixed at tau inside gap pos.
jp_continuous_search <- function(x, y, S, cand, cons, k, extra = NULL,
                                 fix = NULL, min_end = 2L, min_between = 2L) {
  ux <- sort(unique(x))
  n <- length(x)
  lt <- findInterval(cand, sort(x), left.open = TRUE)
  le <- findInterval(cand, sort(x))
  need0 <- c(min_end, rep(min_between, max(k - 1L, 0L)), min_end)
  # Observations strictly inside each segment must be completed by those on
  # its joinpoints, each given to the left segment only when it needs it.
  admissible <- function(i) {
    short <- need0 - c(lt[i[1L]], if (k > 1L) lt[i[-1L]] - le[i[-k]],
                       n - le[i[k]])
    on <- c(le[i] - lt[i], 0)
    spare <- 0
    for (s in seq_along(short)) {
      if (short[s] <= spare) {
        spare <- on[s]
      } else if (short[s] <= spare + on[s]) {
        spare <- 0
      } else {
        return(FALSE)
      }
    }
    TRUE
  }
  left <- findInterval(cand, ux)
  gap <- !(cand %in% ux)
  right <- ux[left + gap]
  U <- jp_whiten_mat(S)
  Z <- jp_apply_whiten(U, cbind(1, x - mean(x), extra))
  H <- jp_apply_whiten(U, outer(x, right, function(a, b) pmax(a - b, 0)))
  St <- jp_apply_whiten(U, outer(x, right, function(a, b) as.numeric(a >= b)))
  yw <- drop(jp_apply_whiten(U, y))
  p0 <- ncol(Z)
  cf <- jp_configs(cons, k)
  if (!is.null(fix)) {
    cf <- cf[, cf[fix$j, ] == fix$pos, drop = FALSE]
    gap[fix$pos] <- FALSE
    H[, fix$pos] <- jp_apply_whiten(U, pmax(x - fix$tau, 0))
    right[fix$pos] <- fix$tau
  }
  sse <- rep(Inf, ncol(cf))
  taus <- matrix(NA_real_, k, ncol(cf))
  for (c in seq_len(ncol(cf))) {
    i <- cf[, c]
    if (k && !admissible(i)) next
    g <- gap[i]
    X <- cbind(Z, H[, i, drop = FALSE], St[, i[g], drop = FALSE])
    ft <- stats::.lm.fit(X, yw)
    if (ft$rank < ncol(X)) next
    b <- ft$coefficients
    tau <- right[i]
    tau[g] <- tau[g] - b[p0 + k + seq_len(sum(g))] / b[p0 + which(g)]
    if (all(is.finite(tau)) && all(tau[g] > ux[left[i[g]]] & tau[g] < right[i[g]])) {
      sse[c] <- sum(ft$residuals^2)
      taus[, c] <- tau
    }
  }
  if (!is.null(fix)) {
    return(min(c(Inf, sse)))
  }
  if (!any(is.finite(sse))) {
    return(NULL)
  }
  at <- matrix(rep(cand, k), length(cand), k)
  o <- order(sse)
  for (j in seq_len(k)) {
    first <- o[!duplicated(cf[j, o])]
    at[cf[j, first], j] <- taus[j, first]
  }
  list(tau = taus[, o[1L]], profile = jp_config_profile(sse, cf, length(cand)),
       profile_at = at)
}

# Joinpoint intervals of continuous fits: an accepted region ends inside a
# gap, where the profile SSE crosses the cutoff `crit`; the crossing is
# found by root finding in the gap at each end of the accepted candidates.
jp_continuous_limits <- function(tab, m, x, y, S, cand, cons, extra, crit,
                                 ctrl) {
  if (!is.finite(crit)) {
    return(tab)
  }
  ux <- sort(unique(x))
  gap <- !(cand %in% ux)
  for (j in seq_len(m$k)) {
    ok <- which(m$profile[, j] <= crit * (1 + 1e-12))
    for (side in c(-1, 1)) {
      e <- if (side < 0) ok[1L] else ok[length(ok)]
      g <- if (gap[e]) e else e + side
      if (g < 1L || g > length(cand) || !gap[g]) next
      inner <- if (g == e) m$profile_at[e, j] else cand[e]
      outer <- ux[findInterval(cand[g], ux) + (side > 0)]
      f <- function(t) {
        jp_continuous_search(x, y, S, cand, cons, m$k, extra,
                             list(j = j, pos = g, tau = t), ctrl$min_obs_end,
                             ctrl$min_obs_between) - crit
      }
      if (f(inner) > 0) next
      lim <- if (f(outer) <= 0) outer else {
        stats::uniroot(f, sort(c(inner, outer)), tol = 1e-8 * diff(range(x)))$root
      }
      if (side < 0) tab$lower[j] <- min(tab$lower[j], lim)
      else tab$upper[j] <- max(tab$upper[j], lim)
    }
  }
  tab
}
