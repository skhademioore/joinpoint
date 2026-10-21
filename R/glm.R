# Likelihood based joinpoint models for counts: Poisson regression with log
# link and the population as exposure, or binomial regression with logit
# link (Czajkowski et al. 2008). The maximum likelihood joinpoints are found
# by fitting every admissible set of joinpoints.

# Response, prior weights and offset of the generalized linear model. The
# Poisson offset makes exp(eta) the rate per `rate_per`.
jp_glm_data <- function(inc, family, rate_per) {
  if (family == "poisson") {
    list(y = inc$count, w = rep(1, nrow(inc)),
         off = log(inc$population / rate_per))
  } else {
    list(y = inc$count / inc$population, w = inc$population,
         off = rep(0, nrow(inc)))
  }
}

# Maximum likelihood fit with fixed joinpoints; coefficients refer to the
# original x. `ok` is FALSE when the estimates do not exist or were not
# found: no convergence, a singular design, or fitted values of 0 (or 1),
# as under complete separation.
jp_glm_fit <- function(x, gd, family, tau, extra = NULL) {
  ctr <- mean(x)
  X <- jp_design(x - ctr, tau - ctr, extra)
  fam <- jp_glm_family(family)
  ft <- suppressWarnings(stats::glm.fit(X, gd$y, weights = gd$w,
                                        offset = gd$off, family = fam))
  b <- ft$coefficients
  b[is.na(b)] <- 0
  eta <- drop(X %*% b)
  mu_eta <- fam$mu.eta(eta + gd$off)
  mu <- fam$linkinv(eta + gd$off)
  eps <- 10 * .Machine$double.eps
  coef <- b
  coef[1L] <- b[1L] - b[2L] * ctr
  names(coef) <- colnames(X)
  list(coef = coef, eta = eta, W = gd$w * mu_eta^2 / fam$variance(mu),
       z = eta + (gd$y - mu) / mu_eta, deviance = ft$deviance,
       pearson = sum(gd$w * (gd$y - mu)^2 / fam$variance(mu)), X = X,
       ctr = ctr, tau = tau,
       ok = ft$converged && !ft$boundary && ft$rank == ncol(X) &&
         all(mu > eps) && (family == "poisson" || all(mu < 1 - eps)))
}

jp_glm_family <- function(family) {
  if (family == "poisson") stats::poisson() else stats::binomial()
}

# Exhaustive search for k joinpoints: the deviance of every admissible set.
# The profile holds the smallest deviance with each joinpoint at each
# candidate.
jp_glm_search <- function(x, gd, family, cand, cons, k, extra = NULL) {
  cf <- jp_configs(cons, k)
  if (!ncol(cf)) {
    return(NULL)
  }
  fam <- jp_glm_family(family)
  u <- x - mean(x)
  base <- cbind(1, u, extra)
  H <- outer(u, cand - mean(x), function(a, b) pmax(a - b, 0))
  dev <- vapply(seq_len(ncol(cf)), function(i) {
    suppressWarnings(stats::glm.fit(cbind(base, H[, cf[, i], drop = FALSE]),
                                    gd$y, weights = gd$w, offset = gd$off,
                                    family = fam)$deviance)
  }, numeric(1))
  ft <- jp_glm_fit(x, gd, family, cand[cf[, which.min(dev)]], extra)
  ft$profile <- jp_config_profile(dev, cf, length(cand))
  ft
}

# Model object with Wald inference (dispersion `phi`) in the same form as the
# least squares models.
jp_glm_model <- function(ft, x, family, phi, jump = NULL) {
  k <- length(ft$tau)
  tau <- ft$tau
  p <- ncol(ft$X)
  Vc <- if (ft$ok) {
    phi * chol2inv(chol(crossprod(ft$X * sqrt(ft$W))))
  } else {
    matrix(NA_real_, p, p)
  }
  # Segment intercepts (original x) and slopes from the centered coefficients.
  nseg <- k + 1L
  ne <- p - 2L - k
  A <- matrix(0, 2L * nseg + ne, p)
  for (j in seq_len(nseg)) {
    A[2L * j - 1L, 1:2] <- c(1, -ft$ctr)
    A[2L * j, 2L] <- 1
    for (h in seq_len(j - 1L)) {
      A[2L * j - 1L, 2L + h] <- -tau[h]
      A[2L * j, 2L + h] <- 1
    }
  }
  for (e in seq_len(ne)) A[2L * nseg + e, 2L + k + e] <- 1
  V <- A %*% Vc %*% t(A)
  si <- 2L * seq_len(nseg)
  n <- length(x)
  list(
    k = k, tau = tau, coef = ft$coef, slopes = jp_slopes(ft$coef, k),
    intercepts = jp_intercepts(ft$coef, tau),
    se_slopes = sqrt(diag(V)[si]), se_intercepts = sqrt(diag(V)[si - 1L]),
    vcov_gp = V, vcov_slopes = V[si, si, drop = FALSE], df = Inf,
    sse = ft$deviance, mse = phi, fitted = ft$eta,
    resid = ft$z - ft$eta, wresid = (ft$z - ft$eta) * sqrt(ft$W), n = n,
    bounds = c(min(x), tau, max(x)), jump = jump,
    extra_names = colnames(ft$X)[-seq_len(2L + k)], log = TRUE,
    link = if (family == "poisson") "log" else "logit", z = ft$z, W = ft$W,
    profile = ft$profile, converged = ft$ok
  )
}

# Best fitting models for kmin..kmax with a common dispersion: estimated
# from the Pearson statistic of the largest model, or 1.
jp_glm_models <- function(x, gd, family, cand, cons, kmin, kmax, extra,
                          jump, overdispersion) {
  fits <- lapply(kmin:kmax, function(k) {
    jp_glm_search(x, gd, family, cand, cons, k, extra)
  })
  ok <- !vapply(fits, is.null, logical(1))
  fits <- fits[ok]
  ks <- (kmin:kmax)[ok]
  bad <- ks[!vapply(fits, `[[`, logical(1), "ok")]
  if (length(bad)) {
    warning("No maximum likelihood estimate for ", paste(bad, collapse = ", "),
            " joinpoints (no convergence, complete separation or fitted ",
            "values of 0); Wald and joinpoint intervals withheld.",
            call. = FALSE)
  }
  last <- fits[[length(fits)]]
  # Pearson dispersion of the largest model; its degrees of freedom also
  # count the estimated joinpoint locations.
  pearson <- last$pearson / (length(x) - ncol(last$X) - max(ks))
  phi <- if (overdispersion) max(1, pearson) else 1
  models <- lapply(fits, jp_glm_model, x = x, family = family, phi = phi,
                   jump = jump)
  names(models) <- ks
  attr(models, "pearson") <- pearson
  models
}

# Information criteria with the deviance divided by n and by the dispersion
# in place of ln(SSE / n).
jp_glm_criteria <- function(x, models, method, extra) {
  n <- length(x)
  ks <- vapply(models, `[[`, numeric(1), "k")
  dev <- vapply(models, `[[`, numeric(1), "sse")
  sc <- dev / models[[1L]]$mse
  tab <- data.frame(joinpoints = ks, deviance = dev,
                    bic = sc / n + 2 * (ks + 1) * log(n) / n,
                    bic3 = sc / n + (3 * ks + 2) * log(n) / n)
  if (method == "wbic") {
    tab$weight <- vapply(models, function(m) {
      if (m$k == 0) NA_real_ else {
        max(jp_local_stats(x, m$z, 1 / m$W, m, extra = extra)$r2)
      }
    }, numeric(1))
    tab$wbic <- ifelse(is.na(tab$weight), tab$bic,
                       (1 - tab$weight) * tab$bic + tab$weight * tab$bic3)
  }
  tab
}

# Sequential Monte Carlo likelihood ratio test: as the permutation test, with
# data sets simulated from the fitted null model (parametric bootstrap) and
# fitted by the same exhaustive search.
jp_glm_lrt <- function(x, gd, family, models, cand, cons, kmin, kmax, alpha,
                       B, extra, n_exceed = NULL) {
  ks <- vapply(models, `[[`, numeric(1), "k")
  a <- kmin
  b <- kmax
  rows <- list()
  while (a < b) {
    m0 <- models[[which(ks == a)]]
    stat <- m0$sse - models[[which(ks == b)]]$sse
    mu <- jp_glm_family(family)$linkinv(m0$fitted + gd$off)
    mc <- jp_mc_pvalue(stat, function(nb) vapply(seq_len(nb), function(i) {
      gs <- gd
      gs$y <- if (family == "poisson") {
        stats::rpois(length(mu), mu)
      } else {
        stats::rbinom(length(mu), gd$w, mu) / gd$w
      }
      f0 <- jp_glm_search(x, gs, family, cand, cons, a, extra)
      f1 <- jp_glm_search(x, gs, family, cand, cons, b, extra)
      if (is.null(f0) || is.null(f1)) NA_real_ else f0$deviance - f1$deviance
    }, numeric(1)), B, n_exceed, batch = if (is.null(n_exceed)) B else 1L)
    lev <- alpha / (kmax - a)
    reject <- mc[["p"]] < lev
    rows[[length(rows) + 1L]] <- data.frame(
      null = a, alternative = b, statistic = stat, p_value = mc[["p"]],
      alpha = lev, reject = reject, simulations = mc[["n"]]
    )
    if (reject) a <- a + 1 else b <- b - 1
  }
  list(k = a, table = do.call(rbind, rows))
}
