#' Control settings for joinpoint fitting
#'
#' Settings of the grid search, model selection, confidence intervals,
#' AAPC ranges and random number generation.
#'
#' @param min_obs_end Minimum number of observations from a joinpoint to
#'   either end of the data, excluding the joinpoint itself when it falls on
#'   an observation (1 to 9, default 2).
#' @param min_obs_between Minimum number of observations between two
#'   joinpoints, excluding joinpoints that fall on observations (0 to 9,
#'   default 2).
#' @param n_grid Number of grid points placed between adjacent observed x
#'   values in the grid search (0 to 9, default 0).
#' @param alpha Overall significance level of the sequential permutation test.
#' @param n_perm Number of randomly permuted data sets for the permutation
#'   test (default 4499).
#' @param n_exceed Optional early stopping of the permutation and likelihood
#'   ratio tests (Besag and Clifford 1991): simulation stops once
#'   `n_exceed` simulated statistics reach the observed one, at simulation
#'   `L`, and the p value is `n_exceed / L`. `NULL` runs all `n_perm`
#'   simulations. Values of 10 to 20 save most of the time of tests that do
#'   not reject.
#' @param perm_seed Random number seed for the permutation test (`NULL` for
#'   the current random number stream).
#' @param n_resample Number of resampled data sets for the empirical quantile
#'   confidence intervals (default 5001).
#' @param eq_seed Random number seed for the empirical quantile method
#'   (`NULL` for the current random number stream).
#' @param eq_final_only Use the empirical quantile method only for the final
#'   selected model; other models get parametric intervals.
#' @param seed_mode `"constant"` uses the same seed for every cohort;
#'   `"vary"` uses `seed * g + (g - 1)` (modulo `2^31 - 1`) for cohort `g`.
#' @param alpha_apc,alpha_aapc,alpha_joinpoint,alpha_jump Significance levels
#'   for the APC, AAPC, joinpoint location and jump or comparability ratio
#'   intervals.
#' @param aapc_full Compute the AAPC over the full range of the data.
#' @param aapc_last Integer vector (up to two values) of "Last X observations"
#'   AAPC ranges.
#' @param aapc_ranges List of up to seven numeric `c(start, end)` ranges
#'   anywhere within the data range.
#' @param madwd Minimum APC difference worth detecting, in percentage points
#'   (`NULL` to disable).
#' @param dds_cutoffs Cutoffs `c < d` of the data dependent selection.
#' @param infeasible What to do when a cohort has too few observations for
#'   the requested `kmax` under the placement constraints: `"reduce"` lowers `kmax` to the largest feasible
#'   value, `"recommended"` uses the recommended maximum for the cohort, and
#'   `"skip"` leaves the cohort out. A warning is given in every case.
#' @param tau_interval Parametric interval of joinpoint locations:
#'   `"lerman"` (Lerman 1980) or the more conservative `"literal"` variant
#'   `Q (1 + k / (n - p)) F`, which reproduces some published intervals.
#' @param time_budget Approximate computing time, in seconds, allowed for
#'   the exhaustive search when `kmax` is left to its recommended value. The
#'   recommended maximum is lowered, with a message, when the search (with
#'   its permutation tests or resamples) would take longer; the budget is an
#'   estimate, not a deadline, and `Inf` removes it. An explicit `kmax` is
#'   always honored.
#' @return A list of class `jp_control`.
#' @references Besag J, Clifford P (1991). Sequential Monte Carlo p-values.
#'   Biometrika 78:301-304.
#' @export
#' @examples
#' jp_control(n_grid = 1, aapc_last = c(5, 10))
jp_control <- function(min_obs_end = 2, min_obs_between = 2, n_grid = 0,
                       alpha = 0.05, n_perm = 4499, n_exceed = NULL,
                       perm_seed = 7160,
                       n_resample = 5001, eq_seed = 10000,
                       eq_final_only = TRUE,
                       seed_mode = c("constant", "vary"),
                       alpha_apc = 0.05, alpha_aapc = 0.05,
                       alpha_joinpoint = 0.05, alpha_jump = 0.05,
                       aapc_full = TRUE, aapc_last = NULL,
                       aapc_ranges = NULL, madwd = NULL,
                       dds_cutoffs = c(10, 200),
                       infeasible = c("reduce", "recommended", "skip"),
                       tau_interval = c("lerman", "literal"),
                       time_budget = 60) {
  seed_mode <- match.arg(seed_mode)
  infeasible <- match.arg(infeasible)
  tau_interval <- match.arg(tau_interval)
  whole <- function(v, lo, hi) {
    length(v) == 1L && is.finite(v) && v == round(v) && v >= lo && v <= hi
  }
  prob <- function(v) length(v) == 1L && is.finite(v) && v > 0 && v < 1
  seed <- function(v) is.null(v) || whole(v, -.Machine$integer.max, .Machine$integer.max)
  flag <- function(v) isTRUE(v) || isFALSE(v)
  ok <- c(
    min_obs_end = whole(min_obs_end, 1, 9),
    min_obs_between = whole(min_obs_between, 0, 9),
    n_grid = whole(n_grid, 0, 9),
    n_perm = whole(n_perm, 1, Inf), n_resample = whole(n_resample, 10, Inf),
    n_exceed = is.null(n_exceed) || whole(n_exceed, 1, Inf),
    perm_seed = seed(perm_seed), eq_seed = seed(eq_seed),
    eq_final_only = flag(eq_final_only), aapc_full = flag(aapc_full),
    alpha = prob(alpha), alpha_apc = prob(alpha_apc),
    alpha_aapc = prob(alpha_aapc), alpha_joinpoint = prob(alpha_joinpoint),
    alpha_jump = prob(alpha_jump),
    aapc_last = length(aapc_last) <= 2 &&
      all(vapply(aapc_last, whole, logical(1), lo = 2, hi = Inf)),
    aapc_ranges = length(aapc_ranges) <= 7 &&
      all(vapply(aapc_ranges, function(r) {
        length(r) == 2L && all(is.finite(r)) && r[1L] < r[2L]
      }, logical(1))),
    madwd = is.null(madwd) ||
      (length(madwd) == 1L && is.finite(madwd) && madwd >= 0),
    dds_cutoffs = length(dds_cutoffs) == 2L && all(is.finite(dds_cutoffs)) &&
      dds_cutoffs[1L] > 0 && dds_cutoffs[1L] < dds_cutoffs[2L],
    time_budget = is.numeric(time_budget) && length(time_budget) == 1L &&
      !is.na(time_budget) && time_budget > 0
  )
  if (!all(ok)) {
    stop("Invalid control settings: ", paste(names(ok)[!ok], collapse = ", "),
         ".", call. = FALSE)
  }
  structure(as.list(environment())[c(
    "min_obs_end", "min_obs_between", "n_grid", "alpha", "n_perm",
    "n_exceed", "perm_seed", "n_resample", "eq_seed", "eq_final_only", "seed_mode",
    "alpha_apc", "alpha_aapc", "alpha_joinpoint", "alpha_jump", "aapc_full",
    "aapc_last", "aapc_ranges", "madwd", "dds_cutoffs", "infeasible",
    "tau_interval", "time_budget"
  )], class = "jp_control")
}

#' Fit joinpoint regression models
#'
#' Fits joinpoint (segmented line) regression models to one or more trend
#' series (Kim et al. 2000). For each cohort (combination of `by` values)
#' the best fitting model with `kmin` to `kmax` joinpoints is found by an
#' exhaustive grid search, the number of joinpoints is selected with the
#' chosen criterion, and the annual percent change (APC), average annual
#' percent change (AAPC) and joinpoint locations are estimated with
#' confidence intervals.
#'
#' @param data A data frame.
#' @param x Independent variable (column name, bare or quoted), usually the
#'   calendar year. Values must be unique within each cohort.
#' @param y Dependent variable (a rate, proportion, count or other value).
#'   May be omitted when `count` and `population` are given.
#' @param by Optional character vector of by-variables defining cohorts.
#' @param se Optional column of standard errors of `y` (on the original
#'   scale).
#' @param count,population Optional columns of counts and populations (or
#'   numerators and denominators) used to calculate crude rates, proportions
#'   or percents and their standard errors.
#' @param type Type of dependent variable: `"rate"` (a provided crude or age
#'   adjusted rate), `"crude_rate"`, `"proportion"`, `"percent"`, `"count"`
#'   or `"other"`.
#' @param rate_per Rate multiplier for calculated crude rates.
#' @param family `"gaussian"` (default) fits by (weighted or generalized)
#'   least squares. `"poisson"` fits a Poisson regression of `count` with
#'   log link and `population` as exposure, and `"binomial"` a logistic
#'   regression of `count` events out of `population` trials, both by
#'   maximum likelihood (Czajkowski et al. 2008); for these the segment
#'   changes are percent changes in the rate or in the odds.
#' @param log Fit the log-linear model `ln(y) = xb` (`TRUE`, default) or the
#'   linear model `y = xb` (`FALSE`).
#' @param variance Heteroscedastic errors option: `"constant"`, `"se"`
#'   (standard errors provided or calculated) or `"poisson"`. Defaults to
#'   `"se"` when standard errors are available and `"constant"` otherwise.
#' @param autocorrelation First order autocorrelation of the errors: a number
#'   in (-1, 1), or `"estimate"` to estimate it from the data.
#' @param vcov Optional variance-covariance matrix of the responses (single
#'   cohort only).
#' @param kmin,kmax Minimum and maximum number of joinpoints. `kmax = NULL`
#'   uses the recommended maximum for the number of observations.
#' @param method Model selection method: `"wbic"` (weighted BIC, default),
#'   `"bic"`, `"bic3"`, `"mbic"` (modified BIC), `"dds"` (data dependent
#'   selection) or `"permutation"` for least squares fits; `"wbic"`,
#'   `"bic"`, `"bic3"` or `"lrt"` (sequential Monte Carlo likelihood ratio
#'   test) for Poisson and binomial fits.
#' @param ci_method Confidence interval method for APCs, AAPCs and joinpoint
#'   locations: `"empirical"` (empirical quantile, the default for least
#'   squares fits with grid joinpoints) or `"parametric"` (the only method for
#'   Poisson and binomial fits and for continuous joinpoints). Parametric
#'   intervals condition on the estimated joinpoints, except the interval of
#'   a joinpoint location; for Poisson and binomial fits they are Wald
#'   intervals and the joinpoint interval is based on the profile deviance.
#' @param covariates Optional one-sided formula of covariates adjusted for in
#'   every segment, for example `~ factor(month)` for seasonality. The
#'   covariates are centered, so the fitted trend refers to their average.
#' @param continuous Estimate joinpoints anywhere between observations
#'   (Hudson 1966) instead of only at the grid points: every placement of
#'   joinpoints at observations or inside the gaps between them is examined,
#'   so the least squares solution is exact. Placements are those allowed by
#'   the placement rules and their limits, in which an observation on a
#'   joinpoint counts toward one of its two segments. The limits of the
#'   parametric joinpoint intervals are located inside the gaps. Available
#'   for least squares fits with parametric intervals and selection by a
#'   criterion; `n_grid` is not used.
#' @param overdispersion For Poisson and binomial fits, estimate a
#'   dispersion parameter (quasi-likelihood) instead of fixing it at 1
#'   (maximum likelihood): the Pearson statistic of the largest model divided
#'   by the number of observations less its coefficients and joinpoints, and
#'   at least 1. It scales the Wald covariances and the deviances in the
#'   criteria and is reported by [model_selection()].
#' @param exclude Values of `x` to exclude from model fitting. Excluded
#'   points are kept in the data and graphs.
#' @param shift Constant added to `x` before fitting.
#' @param jump For coding changes: the last value of `x` under the old
#'   coding. Fits the jump model, or the comparability ratio model when
#'   `comparability_ratio` is given.
#' @param comparability_ratio,cr_variance Comparability ratio (new coding
#'   over old coding) and its variance for the comparability ratio model.
#' @param interval Interval label of `x` used in output, e.g. `"Annual"` or
#'   `"Monthly"`.
#' @param control Settings from [jp_control()].
#' @return An object of class `joinpoint`.
#' @references
#' Kim HJ, Fay MP, Feuer EJ, Midthune DN (2000). Permutation tests for
#' joinpoint regression with applications to cancer rates. Statistics in
#' Medicine 19:335-351.
#'
#' Clegg LX, Hankey BF, Tiwari R, Feuer EJ, Edwards BK (2009). Estimating
#' average annual percent change in trend analysis. Statistics in Medicine
#' 28:3670-3682.
#'
#' Kim HJ, Luo J, Chen HS, et al. (2017). Improved confidence interval for
#' average annual percent change in trend analysis. Statistics in Medicine
#' 36:3059-3074.
#'
#' Kim HJ, Chen HS, Midthune D, et al. (2023). Data-driven choice of a model
#' selection method in joinpoint regression. Journal of Applied Statistics
#' 50:1992-2013.
#'
#' Hudson DJ (1966). Fitting segmented curves whose join points have to be
#' estimated. Journal of the American Statistical Association 61:1097-1129.
#'
#' Czajkowski M, Gill R, Rempala G (2008). Model selection in logistic
#' joinpoint regression with applications to analyzing cohort mortality
#' patterns. Statistics in Medicine 27:1508-1526.
#' @export
#' @examples
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, x = year, y = rate, kmax = 3, ci_method = "parametric")
#' fit
#' apc(fit)
joinpoint <- function(data, x, y = NULL, by = NULL, se = NULL, count = NULL,
                      population = NULL,
                      type = c("rate", "crude_rate", "proportion", "percent",
                               "count", "other"),
                      rate_per = 100000,
                      family = c("gaussian", "poisson", "binomial"),
                      log = TRUE, variance = NULL, autocorrelation = 0,
                      vcov = NULL, kmin = 0, kmax = NULL,
                      method = c("wbic", "bic", "bic3", "mbic", "dds",
                                 "permutation", "lrt"),
                      ci_method = NULL, covariates = NULL,
                      continuous = FALSE, overdispersion = FALSE,
                      exclude = NULL, shift = 0, jump = NULL,
                      comparability_ratio = NULL, cr_variance = 0,
                      interval = "Annual", control = jp_control()) {
  cl <- match.call()
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  xn <- jp_col(substitute(x), data, "x")
  yn <- jp_col(substitute(y), data, "y")
  sen <- jp_col(substitute(se), data, "se")
  cn <- jp_col(substitute(count), data, "count")
  pn <- jp_col(substitute(population), data, "population")
  byn <- if (is.null(by)) NULL else jp_col(by, data, "by")
  type <- match.arg(type)
  family <- match.arg(family)
  method <- match.arg(method)
  glm <- family != "gaussian"
  ci_method <- match.arg(ci_method %||%
                           if (glm || continuous) "parametric" else "empirical",
                         c("empirical", "parametric"))
  if (is.null(yn) && !is.null(cn) && !is.null(pn) && type == "rate") {
    type <- "crude_rate"
  }
  if (glm) {
    if (is.null(cn) || is.null(pn) || !is.null(yn) || !log ||
        !is.null(vcov) || !identical(autocorrelation, 0) ||
        ci_method != "parametric" || continuous ||
        !(method %in% c("wbic", "bic", "bic3", "lrt")) ||
        !is.null(comparability_ratio)) {
      stop("Poisson and binomial fits need `count` and `population` (and ",
           "no `y`), the log or logit scale, independent errors, parametric ",
           "intervals, grid joinpoints, a method among \"wbic\", \"bic\", ",
           "\"bic3\" and \"lrt\", and no comparability ratio.", call. = FALSE)
    }
    if (method == "lrt" && overdispersion) {
      stop("The likelihood ratio test simulates Poisson or binomial counts ",
           "and cannot allow for overdispersion; use a criterion.",
           call. = FALSE)
    }
    type <- if (family == "poisson") "crude_rate" else if (type == "percent") {
      "percent"
    } else {
      "proportion"
    }
    variance <- "glm"
  } else if (method == "lrt") {
    stop("`method = \"lrt\"` is for Poisson and binomial fits; use ",
         "\"permutation\" for least squares fits.", call. = FALSE)
  }
  if (continuous && (ci_method != "parametric" || method == "permutation")) {
    stop("Continuous joinpoints are available with parametric intervals and ",
         "selection by a criterion.", call. = FALSE)
  }
  if (is.null(variance)) {
    variance <- if (!is.null(sen) ||
                    (type %in% c("crude_rate", "proportion", "percent") &&
                     is.null(yn) && !is.null(cn))) "se" else "constant"
  }
  variance <- match.arg(variance, c("constant", "se", "poisson", "glm"))
  if (!(identical(autocorrelation, "estimate") ||
        (is.numeric(autocorrelation) && abs(autocorrelation) < 1))) {
    stop("`autocorrelation` must be \"estimate\" or a number in (-1, 1).",
         call. = FALSE)
  }
  if ((identical(autocorrelation, "estimate") || autocorrelation != 0) &&
      variance == "poisson") {
    stop("Autocorrelated errors are available with constant variance or ",
         "standard errors only.", call. = FALSE)
  }
  if (!is.null(jump) && !log) {
    stop("Jump and comparability ratio models require `log = TRUE`.",
         call. = FALSE)
  }
  if (!is.null(comparability_ratio)) {
    if (is.null(jump)) {
      stop("A comparability ratio needs `jump`, the last x of the old coding.",
           call. = FALSE)
    }
    if (comparability_ratio < 0.01 || comparability_ratio > 100 ||
        cr_variance < 0 || cr_variance > 100) {
      stop("`comparability_ratio` must lie in [0.01, 100] and `cr_variance` ",
           "in [0, 100].", call. = FALSE)
    }
    if (!is.null(vcov)) {
      stop("The comparability ratio model cannot be combined with `vcov`.",
           call. = FALSE)
    }
  }
  if (!inherits(control, "jp_control")) {
    stop("`control` must be created by jp_control().", call. = FALSE)
  }
  flag <- function(v) isTRUE(v) || isFALSE(v)
  if (!flag(log) || !flag(continuous) || !flag(overdispersion)) {
    stop("`log`, `continuous` and `overdispersion` must be TRUE or FALSE.",
         call. = FALSE)
  }
  nonneg <- function(v) {
    is.numeric(v) && length(v) == 1L && is.finite(v) && v == round(v) && v >= 0
  }
  if (!nonneg(kmin) || !(is.null(kmax) || nonneg(kmax))) {
    stop("`kmin` and `kmax` must be nonnegative whole numbers.", call. = FALSE)
  }
  if (!is.null(kmax) && kmin > kmax) {
    stop("`kmin` must not exceed `kmax`.", call. = FALSE)
  }
  if (!(length(rate_per) == 1L && is.finite(rate_per) && rate_per > 0) ||
      !(length(shift) == 1L && is.finite(shift))) {
    stop("`rate_per` must be positive and `shift` finite.", call. = FALSE)
  }
  if (method == "wbic" && control$min_obs_end == 1 &&
      control$min_obs_between == 0 && control$n_grid == 0) {
    stop("WBIC cannot be used with the grid settings (1, 0, 0).",
         call. = FALSE)
  }
  if (method == "permutation" && control$n_perm < 4499) {
    message("Fewer than 4499 permutations give less consistent p values.")
  }

  if (!is.null(byn) && anyNA(data[byn])) {
    warning(sum(!stats::complete.cases(data[byn])), " rows with missing ",
            "by-variable values are left out.", call. = FALSE)
  }
  groups <- if (is.null(byn)) {
    list(seq_len(nrow(data)))
  } else {
    split(seq_len(nrow(data)), data[byn], drop = TRUE, lex.order = TRUE)
  }
  if (!is.null(vcov) && length(groups) > 1L) {
    stop("A variance-covariance matrix can only be used with one cohort.",
         call. = FALSE)
  }
  keys <- if (is.null(byn)) {
    data.frame(cohort = "All")[, 0, drop = FALSE]
  } else {
    do.call(rbind, lapply(groups, function(i) data[i[1L], byn, drop = FALSE]))
  }
  rownames(keys) <- NULL

  settings <- list(
    x = xn, y = yn, se = sen, count = cn, population = pn, type = type,
    rate_per = rate_per, family = family, log = log, variance = variance,
    covariates = covariates, continuous = continuous,
    overdispersion = overdispersion,
    autocorrelation = autocorrelation, vcov = vcov, kmin = kmin,
    kmax = kmax, method = method, ci_method = ci_method, exclude = exclude,
    shift = shift, jump = jump, comparability_ratio = comparability_ratio,
    cr_variance = cr_variance, interval = interval, control = control
  )
  get <- function(nm, i) if (is.null(nm)) NULL else data[[nm]][i]
  # A cohort whose data cannot be analyzed is skipped with a warning; the call fails only when no cohort can be fitted.
  errs <- rep(NA_character_, length(groups))
  fits <- lapply(seq_along(groups), function(g) {
    i <- groups[[g]]
    tryCatch({
      d <- jp_prepare(
        x = get(xn, i), r = get(yn, i), se = get(sen, i),
        count = get(cn, i), population = get(pn, i), type = type,
        variance = variance, log = log, rate_per = rate_per,
        exclude = exclude, shift = shift
      )
      spec <- settings
      spec$label <- jp_cohort_label(list(by = byn, cohorts = keys), g)
      if (!is.null(covariates)) {
        mf <- stats::model.frame(covariates, data[i, , drop = FALSE],
                                 na.action = stats::na.pass)
        Z0 <- stats::model.matrix(stats::terms(mf), mf)
        Z <- Z0[d$row, colnames(Z0) != "(Intercept)", drop = FALSE]
        if (anyNA(Z)) stop("Missing covariate values.", call. = FALSE)
        ctr <- colMeans(Z[d$included, , drop = FALSE])
        spec$covariates <- sweep(Z, 2L, ctr)
        # Everything needed to rebuild the design for new data: the terms
        # with their data dependent transformations, factor levels,
        # contrasts and centering.
        spec$covariate_info <- list(
          terms = stats::terms(mf), xlev = stats::.getXlevels(stats::terms(mf), mf),
          contrasts = attr(Z0, "contrasts"), center = ctr
        )
      }
      if (!is.null(vcov)) spec$vcov <- vcov[d$row, d$row]
      if (!is.null(jump)) {
        spec$jump <- jump + shift
        if (!is.null(comparability_ratio)) {
          d <- jp_apply_cr(d, spec$jump, comparability_ratio, cr_variance,
                           variance, log)
        }
      }
      jp_fit_cohort(d, spec, g)
    }, error = function(e) {
      errs[g] <<- conditionMessage(e)
      NULL
    })
  })
  if (all(!is.na(errs))) stop(unique(errs), call. = FALSE)
  for (g in which(!is.na(errs))) {
    warning("Cohort ", jp_cohort_label(list(by = byn, cohorts = keys), g),
            " skipped: ", errs[g], call. = FALSE)
  }
  skipped <- if (any(!is.na(errs))) {
    cbind(keys[!is.na(errs), , drop = FALSE], reason = errs[!is.na(errs)])
  }
  structure(list(call = cl, fits = fits, cohorts = keys, by = byn,
                 settings = settings, skipped = skipped), class = "joinpoint")
}

# Comparability ratio model: observations before the coding change are
# multiplied by the ratio and their variances updated by the delta method.
jp_apply_cr <- function(d, last_old, cr, cr_var, variance, log) {
  pre <- d$x <= last_old
  r_adj <- d$r
  r_adj[pre] <- cr * d$r[pre]
  if (variance != "constant") {
    v_orig <- if (log) d$v * d$r^2 else d$v
    v_adj <- v_orig
    v_adj[pre] <- cr^2 * v_orig[pre] + d$r[pre]^2 * cr_var
    d$v <- if (log) v_adj / r_adj^2 else v_adj
  }
  d$y <- if (log) base::log(r_adj) else r_adj
  d
}
