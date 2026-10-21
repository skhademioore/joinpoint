# Model selection criteria (Kim et al. 2009; Zhang and Siegmund 2007;
# Kim and Kim 2016; Kim et al. 2023).

jp_bic <- function(sse, n, k) log(sse / n) + 2 * (k + 1) * log(n) / n

jp_bic3 <- function(sse, n, k) log(sse / n) + (3 * k + 2) * log(n) / n

# Modified BIC adapted from Zhang and Siegmund (2007) to segmented line
# regression. The determinant uses the whitened hinge design.
jp_mbic <- function(sse, n, k, Xw) {
  ld <- determinant(crossprod(Xw), logarithm = TRUE)$modulus
  jp_bic(sse, n, k) + as.numeric(ld) / n -
    (2 / n) * lgamma((n - k - 3) / 2) - ((k + 3) / n) * log(sse)
}

# Local statistics for adjacent segment pairs used by DDS and WBIC. For the
# i-th joinpoint the observations in (tau_{i-1}, tau_{i+1}] are used, with
# tau_0 = min(x) - 1 and tau_{k+1} = max(x). Returns the partial R^2 and the
# Delta measure for every joinpoint.
jp_local_stats <- function(x, y, S, model, sigma2_max = NA_real_,
                           extra = NULL) {
  k <- model$k
  if (k == 0) {
    return(list(r2 = numeric(0), delta = numeric(0)))
  }
  tau <- model$tau
  ext <- c(min(x) - 1, tau, max(x))
  r2 <- delta <- numeric(k)
  for (i in seq_len(k)) {
    idx <- which(x > ext[i] & x <= ext[i + 2L])
    xi <- x[idx]
    U <- jp_whiten_mat(jp_subset_cov(S, idx))
    X0 <- jp_apply_whiten(U, cbind(1, xi, extra[idx, , drop = FALSE]))
    z <- drop(jp_apply_whiten(U, pmax(xi - tau[i], 0)))
    yy <- drop(jp_apply_whiten(U, y[idx]))
    q0 <- qr(X0)
    rz <- qr.resid(q0, z)
    ry <- qr.resid(q0, yy)
    zz <- sum(rz^2)
    den <- sum(ry^2) * zz
    r2[i] <- if (den > 0) sum(ry * rz)^2 / den else 0
    delta[i] <- model$coef[[2L + i]]^2 * zz / sigma2_max
  }
  list(r2 = r2, delta = delta)
}

# Criteria table for all fitted models.
jp_criteria <- function(x, y, S, models, method, extra = NULL) {
  n <- length(x)
  ks <- vapply(models, `[[`, numeric(1), "k")
  sse <- vapply(models, `[[`, numeric(1), "sse")
  # Exact fits share the numerical zero, so the penalty decides among them.
  tol <- jp_tiny(jp_whiten_mat(S), y)
  ls <- pmax(sse, tol)
  tab <- data.frame(
    joinpoints = ks, sse = sse,
    bic = jp_bic(ls, n, ks), bic3 = jp_bic3(ls, n, ks)
  )
  if (method == "wbic") {
    wt <- vapply(models, function(m) {
      if (m$k == 0) NA_real_ else max(jp_local_stats(x, y, S, m, extra = extra)$r2)
    }, numeric(1))
    tab$weight <- wt
    tab$wbic <- ifelse(is.na(wt), tab$bic, (1 - wt) * tab$bic + wt * tab$bic3)
  }
  if (method == "mbic") {
    U <- jp_whiten_mat(S)
    tab$mbic <- vapply(models, function(m) {
      Xw <- jp_apply_whiten(U, jp_design(x - mean(x), m$tau - mean(x),
                                         extra))
      jp_mbic(max(m$sse, tol), n, m$k, Xw)
    }, numeric(1))
  }
  tab
}

# Data dependent selection (DDS): chooses between the BIC and BIC3 selected
# models using the Delta measure with cutoffs c and d.
jp_select_dds <- function(x, y, S, models, tab, dds_c = 10, dds_d = 200,
                          extra = NULL) {
  kb <- tab$joinpoints[which.min(tab$bic)]
  k3 <- tab$joinpoints[which.min(tab$bic3)]
  ks <- tab$joinpoints
  mmax <- models[[which.max(ks)]]
  sigma2_max <- mmax$mse
  dk <- function(k) {
    if (k == 0) {
      return(0)
    }
    m <- models[[which(ks == k)]]
    min(jp_local_stats(x, y, S, m, sigma2_max, extra)$delta)
  }
  if (kb == k3) {
    return(list(k = kb, k_bic = kb, k_bic3 = k3, delta_bic = NA_real_,
                delta_bic3 = NA_real_))
  }
  db <- dk(kb)
  d3 <- dk(k3)
  dmax <- max(db, d3)
  dmin <- min(db, d3)
  k <- if (dmax <= dds_c || dmax - dmin > dds_d) kb else k3
  list(k = k, k_bic = kb, k_bic3 = k3, delta_bic = db, delta_bic3 = d3)
}

# Minimum APC difference worth detecting: lowers k until all adjacent
# segment APCs differ by at least `psi` percentage points.
jp_madwd <- function(models, k_sel, psi, kmin, log) {
  ks <- vapply(models, `[[`, numeric(1), "k")
  rows <- list()
  k <- k_sel
  repeat {
    m <- models[[which(ks == k)]]
    ch <- if (log) 100 * exp(m$slopes) else m$slopes
    d <- abs(diff(ch))
    pass <- k == 0 || all(d >= psi)
    rows[[length(rows) + 1L]] <- data.frame(
      joinpoints = k,
      min_difference = if (k) min(d) else NA_real_,
      result = if (pass) "Pass" else if (k <= kmin) "Fail (at kmin)" else "Fail"
    )
    if (pass || k <= kmin) break
    k <- k - 1
  }
  list(k = k, table = do.call(rbind, rows))
}
