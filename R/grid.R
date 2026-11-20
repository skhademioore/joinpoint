# Grid search engine helpers.

# Candidate joinpoint locations: every observed x value plus `n_grid` equally
# spaced points between adjacent observed values.
jp_candidates <- function(x, n_grid = 0L) {
  ux <- sort(unique(x))
  if (n_grid <= 0L || length(ux) < 2L) {
    return(ux)
  }
  steps <- seq_len(n_grid) / (n_grid + 1)
  between <- unlist(lapply(seq_len(length(ux) - 1L), function(i) {
    ux[i] + steps * (ux[i + 1L] - ux[i])
  }))
  sort(c(ux, between))
}

# Placement constraints for each candidate: a minimum number of
# observations from a joinpoint to either end of the data and between two
# joinpoints. Observations that coincide with a joinpoint are not counted,
# or, with `closed = TRUE`, counted on both sides (a superset of the
# placements of continuous joinpoints, which jp_continuous_search narrows).
jp_constraints <- function(x, cand, min_end = 2L, min_between = 2L,
                           closed = FALSE) {
  xs <- sort(x)
  n <- length(xs)
  lt <- findInterval(cand, xs, left.open = TRUE)
  le <- findInterval(cand, xs)
  before <- if (closed) le else lt
  after <- n - if (closed) lt else le
  # Observations between candidates i < j: at or beyond i (closed) or
  # beyond i, less those at or beyond j (open) or beyond j.
  from <- n - if (closed) lt else le
  upto <- n - if (closed) le else lt
  m <- length(cand)
  next_min <- vapply(seq_len(m), function(i) {
    if (i == m) {
      return(m)
    }
    j <- (i + 1L):m
    ok <- which(from[i] - upto[j] >= min_between)
    if (length(ok)) j[ok[1L]] - 1L else m
  }, integer(1))
  list(
    first_ok = as.integer(before >= min_end),
    last_ok = as.integer(after >= min_end),
    next_min = as.integer(next_min)
  )
}

# Builds the cross product structures used by the compiled search.
# `base` is the n x p0 matrix of fixed columns (usually intercept and x),
# `xh` the variable used for the hinge columns, `Y` an n x nd matrix of
# responses and `W` either a vector of weights or a full weight matrix
# (inverse covariance).
jp_crossprod <- function(base, xh, cand, Y, W) {
  X <- cbind(base, outer(xh, cand, function(a, b) pmax(a - b, 0)))
  Y <- as.matrix(Y)
  if (is.matrix(W)) {
    WX <- W %*% X
    WY <- W %*% Y
  } else {
    WX <- X * W
    WY <- Y * W
  }
  list(
    G = crossprod(X, WX),
    B = crossprod(X, WY),
    yy = colSums(Y * WY)
  )
}

# Minimum weighted SSE (and its joinpoints) for every number of joinpoints
# 0..kmax, for every column of Y. The fixed columns are an intercept, x and
# the optional `extra` columns (jump indicators, group intercepts).
#
# For numerical stability the responses and the hinge columns are first
# whitened and projected off the fixed columns by a QR decomposition; the
# compiled search then works on these residuals, whose sums of squares are
# the SSEs themselves rather than small differences of large sums (a large
# common trend no longer swamps the differences between candidates). A
# hinge column in the span of the fixed columns (to within rounding) cannot
# enter a model.
jp_grid_search <- function(x, Y, W, cand, cons, kmax, extra = NULL,
                           profile = FALSE) {
  ctr <- mean(x)
  scl <- stats::sd(x)
  if (!is.finite(scl) || scl == 0) scl <- 1
  u <- (x - ctr) / scl
  base <- cbind(1, u, extra)
  Y <- as.matrix(Y)
  if (is.matrix(W)) {
    U <- chol(W)
    wh <- function(A) U %*% A
  } else {
    sw <- sqrt(W)
    wh <- function(A) A * sw
  }
  qb <- qr(wh(base))
  if (qb$rank < ncol(base)) {
    stop("base design columns are linearly dependent", call. = FALSE)
  }
  Yr <- qr.resid(qb, wh(Y))
  yy <- colSums(Yr^2)
  if (kmax <= 0L || !length(cand)) {
    return(list(sse = matrix(yy, ncol = 1L), idx = list(matrix(integer(0),
                                                             ncol(Y), 0L)),
                profile = list(NULL)))
  }
  Hw <- wh(outer(u, (cand - ctr) / scl, function(a, b) pmax(a - b, 0)))
  Hr <- qr.resid(qb, Hw)
  Hr[, colSums(Hr^2) <= 1e-14 * colSums(Hw^2)] <- 0
  jp_grid_cpp(
    crossprod(Hr), crossprod(Yr, Hr), yy, 0L, as.integer(kmax),
    cons$first_ok, cons$last_ok, cons$next_min, profile
  )
}

# Configurations of the exhaustive search at each depth 1..kmax: `nodes`
# counts all admissible placements of the first j joinpoints (the work of
# the search), `models` those that are complete j joinpoint models.
jp_node_counts <- function(cons, kmax) {
  m <- length(cons$first_ok)
  cnt <- as.numeric(cons$first_ok)
  nodes <- models <- numeric(kmax)
  for (j in seq_len(kmax)) {
    if (j > 1L) {
      d <- numeric(m + 1L)
      for (i in which(cons$next_min < m & cnt > 0)) {
        d[cons$next_min[i] + 1L] <- d[cons$next_min[i] + 1L] + cnt[i]
      }
      cnt <- cumsum(d)[seq_len(m)]
    }
    nodes[j] <- sum(cnt)
    models[j] <- sum(cnt[cons$last_ok == 1L])
  }
  list(nodes = nodes, models = models)
}

# Every admissible set of k candidate indices, one set per column.
jp_configs <- function(cons, k) {
  m <- length(cons$first_ok)
  if (k == 0L) {
    return(matrix(integer(0), 0L, 1L))
  }
  cf <- matrix(which(cons$first_ok == 1L), nrow = 1L)
  for (d in seq_len(k - 1L)) {
    nxt <- lapply(cf[d, ], function(i) {
      seq.int(cons$next_min[i] + 1L, length.out = m - cons$next_min[i])
    })
    cf <- rbind(cf[, rep(seq_len(ncol(cf)), lengths(nxt)), drop = FALSE],
                unlist(nxt))
  }
  cf[, cons$last_ok[cf[k, ]] == 1L, drop = FALSE]
}

# Smallest value of `v` (one per configuration) for each joinpoint and
# candidate: the profile used for joinpoint intervals.
jp_config_profile <- function(v, cf, m) {
  prof <- matrix(Inf, m, nrow(cf))
  for (j in seq_len(nrow(cf))) {
    best <- tapply(v, cf[j, ], min)
    prof[as.integer(names(best)), j] <- best
  }
  prof
}
