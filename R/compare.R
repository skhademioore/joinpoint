# Pairwise comparison of joinpoint regression functions (Kim et al. 2004)
# and comparison of AAPCs between independent groups.

# Light fits with fixed joinpoints for every k, used where only SSEs and
# coefficients are needed (comparisons and clustering).
jp_light_models <- function(x, y, S, cand, cons, kmin, kmax, extra = NULL) {
  g <- jp_grid_batch(x, y, S, cand, cons, kmax, extra)
  ks <- (kmin:kmax)[is.finite(g$sse[1L, kmin:kmax + 1L])]
  models <- lapply(ks, function(k) {
    tau <- if (k) cand[g$idx[[k + 1L]][1L, ]] else numeric(0)
    ft <- jp_fit_fixed(x, y, S, tau, extra)
    n <- length(x)
    p <- 2 * (k + 1) + if (is.null(extra)) 0 else ncol(extra)
    c(ft, list(k = k, slopes = jp_slopes(ft$coef, k), n = n,
               mse = ft$sse / (n - p)))
  })
  names(models) <- ks
  models
}

# Selects the number of joinpoints among light models.
jp_light_select <- function(x, y, S, models, method, ctrl, cand, cons,
                            extra = NULL, seed = NULL) {
  if (method == "permutation") {
    ks <- as.numeric(names(models))
    r <- jp_with_seed(seed, jp_permutation_select(
      x, S, models, cand, cons, min(ks), max(ks), ctrl$alpha, ctrl$n_perm,
      extra, ctrl$n_exceed
    ))
    return(r$k)
  }
  crit <- jp_criteria(x, y, S, models, method, extra)
  switch(method,
    dds = jp_select_dds(x, y, S, models, crit, ctrl$dds_cutoffs[1L],
                        ctrl$dds_cutoffs[2L], extra)$k,
    crit$joinpoints[which.min(crit[[method]])]
  )
}

#' Pairwise comparison of joinpoint trends
#'
#' Tests whether the joinpoint regression functions of two cohorts are
#' identical (test of coincidence) or parallel (test of parallelism, equal
#' joinpoints and slopes with different intercepts), using the permutation
#' test of Kim et al. (2004). Every pair of levels of the innermost
#' by-variable is compared within each combination of the other
#' by-variables. The test is carried out with
#' `K = max(K0, K1, K2)` joinpoints, where `K1` and `K2` are the numbers
#' selected for each cohort and `K0` the number selected for the two cohorts
#' fitted together under the null hypothesis. For each pair whose test
#' rejects the null hypothesis (after any p value adjustment) the AAPCs of
#' the two cohorts of a log-linear fit are compared by a z test with
#' parametric variances, also
#' when the fit used empirical quantile intervals; these follow-up p values
#' are not adjusted.
#'
#' @param object A `joinpoint` object with at least one by-variable,
#'   fitted with uncorrelated errors.
#' @param type `"parallel"` or `"coincidence"`.
#' @param alpha Significance level of the test.
#' @param n_perm Number of permuted data sets.
#' @param seed Random number seed.
#' @param p_adjust Method of [stats::p.adjust()] applied to the p values of
#'   the tests of all pairs, for example `"holm"` when many pairs are
#'   compared; the unadjusted p values are kept in `p_unadjusted`.
#' @return An object of class `jp_comparison`: a data frame with one row per
#'   pair and the attributes `aapc` (AAPC comparisons), `fits` (null model
#'   fits) and `p_adjust`.
#' @references Kim HJ, Fay MP, Yu B, Barrett MJ, Feuer EJ (2004).
#'   Comparability of segmented line regression models. Biometrics
#'   60:1005-1014.
#' @export
#' @examples
#' d <- subset(us_mortality, cause %in% c("Stroke", "Heart Disease"))
#' fit <- joinpoint(d, year, rate, by = "cause", kmax = 2,
#'                  ci_method = "parametric")
#' cmp <- jp_compare(fit, n_perm = 199)
#' cmp
jp_compare <- function(object, type = c("parallel", "coincidence"),
                       alpha = 0.05, n_perm = 4499, seed = 7160,
                       p_adjust = "none") {
  stopifnot(inherits(object, "joinpoint"))
  type <- match.arg(type)
  s <- object$settings
  if (is.null(object$by)) {
    stop("Pairwise comparison needs at least one by-variable.", call. = FALSE)
  }
  if (any(vapply(object$fits, function(f) is.matrix(f$S), logical(1)))) {
    stop("Pairwise comparison requires uncorrelated errors.", call. = FALSE)
  }
  if (!is.null(s$jump) || s$family != "gaussian" || !is.null(s$covariates) ||
      s$continuous) {
    stop("Pairwise comparison is available for least squares fits with grid ",
         "joinpoints and without jumps or covariates.", call. = FALSE)
  }
  ctrl <- s$control
  keys <- object$cohorts
  inner <- utils::tail(object$by, 1L)
  outer <- setdiff(object$by, inner)
  strata <- if (length(outer)) {
    interaction(keys[outer], drop = TRUE, lex.order = TRUE)
  } else {
    factor(rep("All", nrow(keys)))
  }
  rows <- list()
  nulls <- list()
  for (st in levels(strata)) {
    ids <- which(strata == st & !vapply(object$fits, is.null, logical(1)))
    if (length(ids) < 2L) next
    for (pr in utils::combn(ids, 2L, simplify = FALSE)) {
      r <- jp_compare_pair(object, pr[1L], pr[2L], type, alpha, n_perm,
                           seed, ctrl)
      lab <- c(group1 = as.character(keys[[inner]][pr[1L]]),
               group2 = as.character(keys[[inner]][pr[2L]]))
      key <- keys[pr[1L], outer, drop = FALSE]
      rows[[length(rows) + 1L]] <- cbind(key, as.data.frame(as.list(lab)),
                                         r$test)
      nulls[[length(nulls) + 1L]] <- c(r$null, list(
        cohorts = pr,
        outer = if (length(outer)) paste(unlist(lapply(key, as.character)),
                                         collapse = ", ")
      ))
    }
  }
  if (!length(rows)) {
    stop("No pair of fitted cohorts to compare.", call. = FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (p_adjust != "none") {
    out$p_unadjusted <- out$p_value
    out$p_value <- stats::p.adjust(out$p_value, p_adjust)
  }
  rej <- out$p_value < alpha
  out$result <- if (type == "parallel") {
    ifelse(rej, "Not parallel", "Parallel")
  } else {
    ifelse(rej, "Not coincident", "Coincident")
  }
  aapcs <- lapply(which(rej & s$log), function(i) {
    pr <- nulls[[i]]$cohorts
    a <- jp_compare_aapc(object$fits[[pr[1L]]], object$fits[[pr[2L]]],
                         ctrl$alpha_aapc)
    if (!is.null(a)) {
      cbind(out[rep(i, nrow(a)), c(outer, "group1", "group2"), drop = FALSE], a)
    }
  })
  aapc <- do.call(rbind, aapcs)
  if (!is.null(aapc)) rownames(aapc) <- NULL
  structure(out, class = c("jp_comparison", "data.frame"), type = type,
            aapc = aapc, fits = nulls, p_adjust = p_adjust, object = object)
}

jp_compare_pair <- function(object, i1, i2, type, alpha, n_perm, seed,
                            ctrl) {
  f1 <- object$fits[[i1]]
  f2 <- object$fits[[i2]]
  if (!isTRUE(all.equal(f1$x, f2$x))) {
    stop("Compared cohorts must have the same independent variable values.",
         call. = FALSE)
  }
  x <- c(f1$x, f2$x)
  y <- c(f1$y, f2$y)
  S <- c(f1$S, f2$S)
  grp <- rep(0:1, c(length(f1$x), length(f2$x)))
  extra <- if (type == "parallel") cbind(group = grp)
  cand <- f1$candidates
  cons <- f1$constraints
  method <- object$settings$method
  kmax <- min(f1$kmax, f2$kmax)

  kmin <- min(f1$kmin, f2$kmin)
  null_models <- jp_light_models(x, y, S, cand, cons, kmin, kmax, extra)
  k0 <- jp_light_select(x, y, S, null_models, method, ctrl, cand, cons, extra,
                        seed)
  K <- min(max(k0, f1$k_selected, f2$k_selected),
           max(as.numeric(names(null_models))))
  sep_fit <- function(yy, xx, SS) {
    jp_grid_batch(xx, yy, SS, cand, cons, K)$sse[, K + 1L]
  }
  m0 <- null_models[[as.character(K)]]
  rss0 <- m0$sse
  rss1 <- sep_fit(f1$y, f1$x, f1$S) + sep_fit(f2$y, f2$x, f2$S)
  n <- length(y)
  d1 <- if (type == "parallel") 2 * K + 1 else 2 * K + 2
  d2 <- n - 2 * (2 * K + 2)
  U <- jp_whiten_mat(S)
  tol <- jp_tiny(U, y)
  fstat <- function(r0, r1) (jp_ratio(r0, r1, tol) - 1) * d2 / d1
  stat <- fstat(rss0, rss1)
  ew <- drop(jp_apply_whiten(U, m0$resid))
  pstat <- jp_with_seed(seed, {
    E <- vapply(seq_len(n_perm), function(b) ew[sample.int(n)], numeric(n))
    Ys <- m0$fitted + jp_unwhiten(U, E)
    r0 <- jp_grid_batch(x, Ys, S, cand, cons, K, extra)$sse[, K + 1L]
    r1 <- sep_fit(Ys[grp == 0, , drop = FALSE], f1$x, f1$S) +
      sep_fit(Ys[grp == 1, , drop = FALSE], f2$x, f2$S)
    fstat(r0, r1)
  })
  p <- (sum(pstat >= stat * (1 - 1e-10)) + 1) / (n_perm + 1)
  test <- data.frame(
    test = type, joinpoints = K, k_null = k0, k_group1 = f1$k_selected,
    k_group2 = f2$k_selected, statistic = stat, df1 = d1, df2 = d2,
    p_value = p, alpha = alpha
  )
  list(test = test,
       null = list(model = null_models[[as.character(k0)]], x = x, grp = grp))
}

# AAPC comparison of two independent groups over the same ranges.
jp_compare_aapc <- function(f1, f2, alpha) {
  a1 <- jp_aapc_table(jp_pick(f1, "final"), f1$ranges, alpha, "parametric")
  a2 <- jp_aapc_table(jp_pick(f2, "final"), f2$ranges, alpha, "parametric")
  if (is.null(a1) || is.null(a2)) {
    return(NULL)
  }
  common <- intersect(a1$range, a2$range)
  a1 <- a1[match(common, a1$range), ]
  a2 <- a2[match(common, a2$range), ]
  z <- stats::qnorm(1 - alpha / 2)
  d <- a1$estimate - a2$estimate
  se_d <- 100 * sqrt(exp(2 * a1$log_estimate) * a1$se_log^2 +
                       exp(2 * a2$log_estimate) * a2$se_log^2)
  zstat <- (a1$log_estimate - a2$log_estimate) /
    sqrt(a1$se_log^2 + a2$se_log^2)
  data.frame(
    range = common, aapc1 = a1$estimate, aapc2 = a2$estimate,
    difference = d, lower = d - z * se_d, upper = d + z * se_d,
    statistic = zstat, p_value = 2 * stats::pnorm(-abs(zstat)),
    significant = 2 * stats::pnorm(-abs(zstat)) < alpha,
    ci_method = "parametric"
  )
}

#' @export
print.jp_comparison <- function(x, ...) {
  cat("Pairwise comparison: test of", attr(x, "type"))
  adj <- attr(x, "p_adjust")
  if (!is.null(adj) && adj != "none") {
    cat(";", adj, "adjusted p values")
  }
  cat("\n\n")
  print(jp_plain(x), digits = 4, row.names = FALSE)
  if (!is.null(attr(x, "aapc"))) {
    cat("\nAAPC comparisons (group 1 minus group 2):\n")
    print(attr(x, "aapc"), digits = 4, row.names = FALSE)
  }
  invisible(x)
}

#' Plot a pairwise comparison
#'
#' Observed values of the compared cohorts with their separately fitted
#' models (solid) and the model fitted under the null hypothesis of
#' parallelism or coincidence (dashed), one panel per pair.
#'
#' @param object A `jp_comparison` object.
#' @param log_y Use a logarithmic y axis.
#' @param ... Ignored.
#' @return A ggplot object.
#' @examples
#' d <- subset(us_mortality, cause %in% c("Stroke", "Heart Disease"))
#' fit <- joinpoint(d, year, rate, by = "cause", kmax = 2,
#'                  ci_method = "parametric")
#' ggplot2::autoplot(jp_compare(fit, n_perm = 99))
#' @export
autoplot.jp_comparison <- function(object, log_y = FALSE, ...) {
  fit <- attr(object, "object")
  nf <- attr(object, "fits")
  pd <- lapply(seq_along(nf), function(i) {
    pr <- nf[[i]]$cohorts
    panel <- paste0(if (!is.null(nf[[i]]$outer)) paste0(nf[[i]]$outer, ": "),
                    object$group1[i], " vs ", object$group2[i], "\n",
                    object$result[i], ", p ", jp_pformat(object$p_value[i]))
    one <- function(g) {
      f <- fit$fits[[g]]
      cv <- jp_curve(f, jp_pick(f, "final"), band = FALSE)
      list(
        pts = data.frame(panel = panel, cohort = jp_cohort_label(fit, g),
                         x = f$data$x, y = f$data$r),
        crv = data.frame(panel = panel, cohort = jp_cohort_label(fit, g),
                         x = cv$x, fit = cv$fit, piece = cv$piece)
      )
    }
    a <- one(pr[1L])
    b <- one(pr[2L])
    m <- nf[[i]]$model
    nul <- data.frame(panel = panel,
                      cohort = rep(c(a$pts$cohort[1L], b$pts$cohort[1L]),
                                   table(nf[[i]]$grp)),
                      x = nf[[i]]$x,
                      fit = if (fit$settings$log) exp(m$fitted) else m$fitted)
    list(pts = rbind(a$pts, b$pts), crv = rbind(a$crv, b$crv), nul = nul)
  })
  pts <- do.call(rbind, lapply(pd, `[[`, "pts"))
  crv <- do.call(rbind, lapply(pd, `[[`, "crv"))
  nul <- do.call(rbind, lapply(pd, `[[`, "nul"))
  g <- ggplot2::ggplot(mapping = ggplot2::aes(x = .data$x,
                                              color = .data$cohort)) +
    ggplot2::geom_point(data = pts, ggplot2::aes(y = .data$y), size = 1.4,
                        alpha = 0.55) +
    ggplot2::geom_line(data = crv, ggplot2::aes(
      y = .data$fit, group = interaction(.data$cohort, .data$piece)
    ), linewidth = 0.9) +
    ggplot2::geom_line(data = nul, ggplot2::aes(y = .data$fit,
                                                group = .data$cohort),
                       linewidth = 0.6, linetype = "22") +
    scale_color_joinpoint(name = NULL) +
    ggplot2::facet_wrap(~panel, scales = "free_y") +
    ggplot2::labs(x = jp_title(fit$settings$x), y = jp_ylab(fit$settings),
                  caption = "Solid: separate fits; dashed: fit under the null hypothesis.") +
    theme_joinpoint()
  if (log_y) g <- g + ggplot2::scale_y_log10()
  g
}

# The comparison table without its class and attributes.
jp_plain <- function(x) {
  structure(x, class = "data.frame", type = NULL, aapc = NULL, fits = NULL,
            object = NULL)
}
