# Error structures. The covariance of the model scale response is stored as
# `S`: a vector of variances (diagonal case, possibly all ones for constant
# variance) or a full covariance matrix (autocorrelated errors or a supplied
# variance-covariance matrix). Only relative values matter for estimation:
# the residual variance scale is estimated from the data, as in weighted
# least squares.

# AR(1) correlation of observations at positions `pos` of the cohort's time
# sequence: excluded observations keep their place, so the lag between the
# observations on either side of an excluded one is two.
jp_ar1_corr <- function(n, rho, pos = seq_len(n)) {
  rho^abs(outer(pos, pos, "-"))
}

# Builds S for the included observations (at positions `pos` of the
# cohort's time sequence).
jp_covariance <- function(v, rho = 0, vcov = NULL, r = NULL, log = TRUE,
                          pos = seq_along(v)) {
  n <- length(v)
  sd <- if (all(is.na(v))) rep(1, n) else sqrt(v)
  if (!is.null(vcov)) {
    if (!is.matrix(vcov) || any(dim(vcov) != n)) {
      stop("`vcov` must be an n x n matrix matching the included ",
           "observations.", call. = FALSE)
    }
    if (!isSymmetric(unname(vcov)) || any(diag(vcov) <= 0) ||
        all(vcov == 0)) {
      stop("`vcov` must be symmetric with positive diagonal.", call. = FALSE)
    }
    return(if (log) vcov / outer(r, r) else vcov)
  }
  if (rho != 0) {
    return(outer(sd, sd) * jp_ar1_corr(n, rho, pos))
  }
  sd^2
}

# Weight (inverse covariance) representation used by the fitting code.
# The weight matrix W = S^-1 and its Cholesky factor U (U'U = W) of the
# last covariance matrix are cached: a fit factors the same matrix many
# times, at a cost cubic in the number of observations.
jp_cov_cache <- new.env(parent = emptyenv())

jp_cov_factors <- function(S) {
  if (!identical(S, jp_cov_cache$S)) {
    W <- chol2inv(chol(S))
    jp_cov_cache$U <- chol(W)
    jp_cov_cache$W <- W
    jp_cov_cache$S <- S
  }
  jp_cov_cache
}

# TRUE for a covariance matrix whose scale is known (it is not multiplied by
# an estimated residual variance).
jp_known_scale <- function(S) isTRUE(attr(S, "known_scale"))

jp_weights <- function(S) {
  if (is.matrix(S)) jp_cov_factors(S)$W else 1 / S
}

jp_subset_cov <- function(S, idx) {
  if (is.matrix(S)) S[idx, idx, drop = FALSE] else S[idx]
}

# Whitening transform U with U' U = W, so that U e has identity covariance.
jp_whiten_mat <- function(S) {
  if (is.matrix(S)) {
    jp_cov_factors(S)$U
  } else {
    1 / sqrt(S)
  }
}

jp_apply_whiten <- function(U, A) {
  if (is.matrix(U)) U %*% A else A * U
}

jp_unwhiten <- function(U, A) {
  if (is.matrix(U)) backsolve(U, A) else A / U
}

# Lag one autocorrelation estimated from the residuals of the model with the
# maximum number of joinpoints fitted under independence (Kim et al. 2000,
# section 2.3). Residuals are standardized by the heteroscedastic weights
# first. The estimate is truncated to (-0.99, 0.99).
#
# With excluded observations the included ones are `gaps` time steps apart
# (positions in the time sequence), and their correlation is rho^gap, the
# model used to build the covariance. rho is then the maximum of the
# Gaussian AR(1) likelihood of the residuals at these gaps (profiled over
# the innovation variance). When every gap is even, rho and -rho give the
# same covariance of the included observations, and the nonnegative value
# is used.
jp_estimate_rho <- function(wresid, tol = 0, gaps = NULL) {
  e <- wresid - mean(wresid)
  n <- length(e)
  if (n < 3 || sum(e^2) <= tol) {
    return(0)
  }
  if (is.null(gaps) || all(gaps == 1)) {
    rho <- sum(e[-1] * e[-n]) / sum(e^2)
    return(max(min(rho, 0.99), -0.99))
  }
  nll <- function(r) {
    rg <- r^gaps
    v <- 1 - rg^2
    ss <- e[1L]^2 + sum((e[-1L] - rg * e[-n])^2 / v)
    n * log(ss / n) + sum(log(v))
  }
  lower <- if (all(gaps %% 2 == 0)) 0 else -0.99
  stats::optimize(nll, c(lower, 0.99))$minimum
}
