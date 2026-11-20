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
jp_estimate_rho <- function(wresid, tol = 0) {
  e <- wresid - mean(wresid)
  n <- length(e)
  if (n < 3 || sum(e^2) <= tol) {
    return(0)
  }
  rho <- sum(e[-1] * e[-n]) / sum(e^2)
  max(min(rho, 0.99), -0.99)
}
