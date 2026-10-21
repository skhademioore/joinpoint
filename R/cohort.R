# Fitting one cohort (one combination of by-variable values).

jp_fit_cohort <- function(d, spec, cohort_index = 1L) {
  ctrl <- spec$control
  inc <- d[d$included, , drop = FALSE]
  x <- inc$x
  y <- inc$y
  n <- length(x)
  log <- spec$log

  # Coding change: jump location halfway between the last old coding
  # observation and the first new coding observation, both included.
  jump_at <- NULL
  if (!is.null(spec$jump)) {
    pos <- match(spec$jump, d$x)
    pos_in <- match(spec$jump, x)
    if (is.na(pos) || pos == nrow(d) || !all(d$included[pos + 0:1]) ||
        pos_in < 3L || pos_in > n - 2L) {
      stop("The last old coding observation and the next observation must ",
           "be included, at least two observations from either end.",
           call. = FALSE)
    }
    jump_at <- (x[pos_in] + x[pos_in + 1L]) / 2
  }
  model_jump <- if (!is.null(jump_at) && is.null(spec$comparability_ratio)) {
    jump_at
  }
  covars <- if (!is.null(spec$covariates)) {
    spec$covariates[d$included, , drop = FALSE]
  }
  extra <- cbind(if (!is.null(model_jump)) cbind(jump = as.numeric(x >= jump_at)),
                 covars)
  glm <- spec$family != "gaussian"

  # Continuous joinpoints: the observations and one point inside each gap.
  cand <- jp_candidates(x, if (spec$continuous) 1L else ctrl$n_grid)
  cons <- jp_constraints(x, cand, ctrl$min_obs_end, ctrl$min_obs_between,
                         closed = spec$continuous)
  lim <- jp_kmax(spec, ctrl, cons, n)
  kmin <- lim$kmin
  kmax <- lim$kmax
  notes <- lim$notes

  V <- if (!is.null(spec$vcov)) spec$vcov[d$included, d$included, drop = FALSE]
  rho <- 0
  if (!identical(spec$autocorrelation, 0)) {
    if (length(unique(round(diff(x), 8))) > 1L) {
      warning("Autocorrelated errors assume equally spaced observations.",
              call. = FALSE)
    }
    rho <- spec$autocorrelation
    if (identical(rho, "estimate")) {
      S0 <- jp_covariance(inc$v, 0, V, inc$r, log)
      g0 <- jp_grid_batch(x, y, S0, cand, cons, kmax, extra)
      tau0 <- if (kmax) cand[g0$idx[[kmax + 1L]][1L, ]] else numeric(0)
      rho <- jp_estimate_rho(jp_fit_fixed(x, y, S0, tau0, extra)$wresid,
                             jp_tiny(jp_whiten_mat(S0), y))
    }
  }
  S <- jp_covariance(inc$v, rho, V, inc$r, log)

  if (glm) {
    gd <- jp_glm_data(inc, spec$family, spec$rate_per)
    models <- jp_glm_models(x, gd, spec$family, cand, cons, kmin, kmax, extra,
                            model_jump, spec$overdispersion)
    if (!spec$overdispersion && attr(models, "pearson") > 1.5) {
      notes <- c(notes, sprintf(paste(
        "the Pearson dispersion of the largest model is %.2f; consider",
        "overdispersion = TRUE."), attr(models, "pearson")))
    }
    ks <- as.numeric(names(models))
    kmax <- max(ks)
    if (spec$type == "percent") for (i in seq_along(models)) models[[i]]$scale <- 100
  } else {
    g <- if (!spec$continuous) {
      jp_grid_search(x, y, jp_weights(S), cand, cons, kmax, extra,
                     profile = TRUE)
    }
    models <- lapply(kmin:kmax, function(k) {
      if (spec$continuous) {
        cs <- jp_continuous_search(x, y, S, cand, cons, k, extra, NULL,
                                   ctrl$min_obs_end, ctrl$min_obs_between)
        if (is.null(cs)) {
          return(NULL)
        }
      } else {
        cs <- list(tau = if (k) cand[g$idx[[k + 1L]][1L, ]] else numeric(0),
                   profile = g$profile[[k + 1L]])
      }
      m <- jp_build_model(x, y, S, cs$tau, log, jump = model_jump,
                          covars = covars)
      m$profile <- cs$profile
      m$profile_at <- cs$profile_at
      m
    })
    names(models) <- kmin:kmax
    models <- Filter(Negate(is.null), models)
    ks <- as.numeric(names(models))
    kmax <- max(ks)
  }

  # Model selection.
  method <- spec$method
  crit <- if (glm) {
    jp_glm_criteria(x, models, method, extra)
  } else {
    jp_criteria(x, y, S, models, method, extra)
  }
  perm <- NULL
  dds <- NULL
  k_sel <- switch(method,
    bic = crit$joinpoints[which.min(crit$bic)],
    bic3 = crit$joinpoints[which.min(crit$bic3)],
    wbic = crit$joinpoints[which.min(crit$wbic)],
    mbic = crit$joinpoints[which.min(crit$mbic)],
    dds = {
      dds <- jp_select_dds(x, y, S, models, crit, ctrl$dds_cutoffs[1L],
                           ctrl$dds_cutoffs[2L], extra)
      dds$k
    },
    permutation = {
      seed <- jp_cohort_seed(ctrl$perm_seed, cohort_index, ctrl$seed_mode)
      perm <- jp_with_seed(seed, jp_permutation_select(
        x, S, models, cand, cons, kmin, kmax, ctrl$alpha, ctrl$n_perm, extra,
        ctrl$n_exceed
      ))
      perm$k
    },
    lrt = {
      seed <- jp_cohort_seed(ctrl$perm_seed, cohort_index, ctrl$seed_mode)
      perm <- jp_with_seed(seed, jp_glm_lrt(
        x, gd, spec$family, models, cand, cons, kmin, kmax, ctrl$alpha,
        ctrl$n_perm, extra, ctrl$n_exceed
      ))
      perm$k
    }
  )
  k_stat <- k_sel
  madwd <- NULL
  if (!is.null(ctrl$madwd) && ctrl$madwd > 0) {
    madwd <- jp_madwd(models, k_sel, ctrl$madwd, kmin, log)
    k_sel <- madwd$k
  }

  # Inference for every model: parametric always; empirical quantile for the
  # final model (or all models).
  ci_method <- spec$ci_method
  eq_for <- if (ci_method == "empirical") {
    if (ctrl$eq_final_only) k_sel else ks
  } else {
    integer(0)
  }
  x_all <- d$x
  ranges <- jp_aapc_ranges(x_all, d$included, ctrl$aapc_full, ctrl$aapc_last,
                           ctrl$aapc_ranges)
  eq_seed <- jp_cohort_seed(ctrl$eq_seed, cohort_index, ctrl$seed_mode)
  for (i in seq_along(models)) {
    m <- models[[i]]
    eq <- NULL
    if (m$k %in% eq_for) {
      eq <- jp_with_seed(eq_seed, jp_eq_resample(x, S, m, cand, cons,
                                                 ctrl$n_resample, extra))
    }
    method_i <- if (is.null(eq)) "parametric" else "empirical"
    m$ci_method <- method_i
    m$apc <- jp_apc_table(m, ctrl$alpha_apc, method_i, eq)
    m$aapc <- jp_aapc_table(m, ranges, ctrl$alpha_aapc, method_i, eq)
    m$joinpoints <- jp_tau_table(m, ctrl$alpha_joinpoint, method_i, eq, cand,
                                 ctrl$tau_interval)
    if (spec$continuous && m$k) {
      m$joinpoints <- jp_continuous_limits(
        m$joinpoints, m, x, y, S, cand, cons, extra,
        jp_tau_crit(m, ctrl$alpha_joinpoint, ctrl$tau_interval), ctrl
      )
    }
    m$jump_table <- jp_jump_table(m, ctrl$alpha_jump)
    m$eq <- eq
    models[[i]] <- m
  }

  f <- structure(list(
    data = d, x = x, y = y, S = S, rho = rho, models = models,
    k_selected = k_sel, k_statistical = k_stat, kmin = kmin, kmax = kmax,
    criteria = crit, permutation = perm, dds = dds, madwd = madwd,
    method = method, ci_method = ci_method, candidates = cand,
    constraints = cons, jump_at = jump_at,
    comparability_ratio = spec$comparability_ratio, notes = notes,
    ranges = ranges, covariates = spec$covariates,
    covariate_info = spec$covariate_info
  ), class = "jp_cohort")
  # Fitted values on the original scale for every observation (excluded
  # points are predicted from the model).
  f$fitted <- lapply(models, function(m) jp_fitted_at(f, m, d$x))
  f
}

# Linear predictor at `xs`: the trend with the jump, plus the effect of the
# centered covariates `Z` (NULL for the trend at their average).
jp_linear <- function(m, xs, Z = NULL) {
  ex <- if (!is.null(m$jump)) cbind(as.numeric(xs >= m$jump))
  eta <- jp_eval(m$coef, m$tau, xs, ex)
  if (is.null(Z)) eta else eta + drop(Z %*% m$coef[colnames(Z)])
}

# Fitted values on the original response scale at `xs`, with the
# comparability ratio back transformation applied.
jp_fitted_at <- function(f, m, xs, Z = NULL) {
  out <- jp_linkinv(m, jp_linear(m, xs, Z))
  cr <- f$comparability_ratio
  if (!is.null(cr)) out[xs < f$jump_at] <- out[xs < f$jump_at] / cr
  out
}

# Number of joinpoints to search. An explicit `kmax` that the placement
# constraints cannot accommodate is handled by the `infeasible` setting; a
# recommended `kmax` is lowered when the exhaustive search would exceed the
# time budget.
jp_kmax <- function(spec, ctrl, cons, n) {
  notes <- character(0)
  rec <- jp_recommended_kmax(n)
  cnt <- jp_node_counts(cons, max(rec, spec$kmax))
  nodes <- c(1, 1 + cumsum(cnt$nodes))
  feasible <- sum(cumprod(cnt$models > 0))
  say <- function(msg) {
    warning("Cohort ", spec$label, ": ", msg, call. = FALSE)
    notes <<- c(notes, msg)
  }
  if (is.null(spec$kmax)) {
    kmax <- min(rec, feasible)
    # Seconds per node of the compiled search (single response, and per
    # response in the batched searches of resampling and permutation), and
    # per model fitted by the likelihood and continuous searches (growing
    # with n), which the likelihood ratio test repeats for every simulated
    # data set. Early stopping is not counted.
    cost <- function(k) {
      nd <- nodes[k + 1L]
      tests <- max(1, k - spec$kmin)
      if (spec$family != "gaussian" || spec$continuous) {
        fits <- 1 + sum(cnt$models[seq_len(k)])
        per <- n * if (spec$continuous) 3e-6 else 1.5e-5
        per * fits * (1 + (spec$method == "lrt") * 2 * ctrl$n_perm * tests)
      } else {
        eq <- (spec$ci_method == "empirical") *
          if (ctrl$eq_final_only) 1 else k - spec$kmin + 1
        nd * 2.5e-7 + nd * 1.5e-8 * (
          eq * ctrl$n_resample +
            (spec$method == "permutation") * ctrl$n_perm * tests
        )
      }
    }
    k0 <- kmax
    while (kmax > spec$kmin && cost(kmax) > ctrl$time_budget) kmax <- kmax - 1
    if (kmax < k0) {
      msg <- sprintf(paste(
        "recommended maximum of %d joinpoints lowered to %d to keep the",
        "exhaustive search within about %s seconds; set kmax to override."
      ), k0, kmax, format(signif(ctrl$time_budget, 2)))
      if (kmax == spec$kmin && spec$method %in% c("permutation", "lrt")) {
        say(paste(msg, "No test is performed; set kmax, or lower n_perm or",
                  "set n_exceed in jp_control()."))
      } else {
        message("Cohort ", spec$label, ": ", msg)
        notes <- c(notes, msg)
      }
    }
  } else {
    kmax <- spec$kmax
    if (kmax > feasible) {
      if (ctrl$infeasible == "skip") {
        stop("too few observations for ", kmax, " joinpoints.", call. = FALSE)
      }
      new <- if (ctrl$infeasible == "recommended") min(rec, feasible) else feasible
      say(sprintf("too few observations for %d joinpoints; using %d.", kmax,
                  new))
      kmax <- new
    }
  }
  kmin <- spec$kmin
  if (kmin > kmax) {
    say(sprintf("minimum number of joinpoints lowered from %d to %d.", kmin,
                kmax))
    kmin <- kmax
  }
  list(kmin = kmin, kmax = kmax, notes = notes)
}

jp_cohort_seed <- function(seed, g, mode) {
  if (is.null(seed)) {
    return(NULL)
  }
  if (identical(mode, "vary")) (seed * g + (g - 1)) %% .Machine$integer.max else seed
}

# Jump value and comparability ratio estimated by the jump model.
jp_jump_table <- function(m, alpha) {
  if (is.null(m$jump)) {
    return(NULL)
  }
  g <- m$coef[["jump"]]
  se <- sqrt(m$vcov_gp[2L * (m$k + 1L) + 1L, 2L * (m$k + 1L) + 1L])
  z <- stats::qnorm(1 - alpha / 2)
  tq <- if (m$df > 0) stats::qt(1 - alpha / 2, m$df) else NA_real_
  cr <- exp(g)
  se_cr <- cr * se
  data.frame(
    statistic = c("Jump value", "Comparability ratio"),
    estimate = c(g, cr),
    se = c(se, se_cr),
    lower = c(g - tq * se, cr - z * se_cr),
    upper = c(g + tq * se, cr + z * se_cr),
    test_statistic = c(g / se, (cr - 1) / se_cr),
    p_value = c(2 * stats::pt(-abs(g / se), m$df),
                2 * stats::pnorm(-abs((cr - 1) / se_cr)))
  )
}

# Inverse link: identity, exp (log-linear and Poisson models) or the
# logistic function scaled to proportions or percents.
jp_linkinv <- function(m, eta) {
  switch(m$link, identity = eta, log = exp(eta),
         logit = (m$scale %||% 1) * stats::plogis(eta))
}
