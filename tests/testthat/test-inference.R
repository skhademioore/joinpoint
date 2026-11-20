L <- function(d, ...) {
  joinpoint(d, x, y, log = FALSE, type = "other", ci_method = "parametric",
            ...)
}
one <- function(f) f$fits[[1]]

test_that("an exactly zero trend is not significant", {
  for (lg in c(TRUE, FALSE)) {
    f <- joinpoint(data.frame(x = 1:12, y = if (lg) 1 else 0), x, y, log = lg,
                   kmax = 0, control = jp_control(n_resample = 50))
    expect_equal(c(apc(f)$p_value, aapc(f)$p_value), c(1, 1))
    expect_false(apc(f)$significant)
  }
  expect_equal(joinpoint:::jp_eq_pvalue(c(-1, 0, 0, 1)), 1)
  expect_equal(joinpoint:::jp_eq_pvalue(c(1, 2, 3, -1)), 0.5)
})

test_that("empirical resampling keeps weighted trends unbiased", {
  x <- 1:12
  s <- rep(c(1, 8), each = 6)
  z <- qr.resid(qr(cbind(1, x) / s), rep(1, 12))
  d <- data.frame(x, s, y = s * z)
  f <- joinpoint(d, x, y, se = s, type = "other", log = FALSE, kmax = 0)
  sl <- one(f)$models[[1]]$eq$slopes
  expect_lt(abs(mean(sl)), 3 * stats::sd(sl) / sqrt(length(sl)))
})

test_that("Poisson variance is computed for supplied crude rates", {
  d <- data.frame(x = 1:9, n = 1:9, p = 100, y = 1000 * (1:9))
  f <- joinpoint(d, x, y, count = n, population = p, type = "crude_rate",
                 variance = "poisson", kmax = 0, ci_method = "parametric")
  expect_equal(one(f)$S, 1 / (1:9))
  expect_error(joinpoint(transform(d, s = 1), x, y, se = s, count = n,
                         population = p, type = "crude_rate",
                         variance = "poisson", kmax = 0), "Poisson")
})

test_that("the comparability ratio uncertainty is shared or per observation", {
  # Old coding values are half the new coding ones (CR = 2).
  x <- 1:12
  d <- data.frame(x, y = exp(1 + .1 * x + .01 * rep(c(-1, 1), 6)) /
                    ifelse(x <= 6, 2, 1), se = .1)
  fit <- function(v, ...) {
    joinpoint(d, x, y, kmax = 0, jump = 6, comparability_ratio = 2,
              cr_variance = v, ci_method = "parametric", ...)
  }
  expect_gt(sqrt(vcov(fit(.4))[2, 2]), sqrt(vcov(fit(0))[2, 2]))
  S <- one(fit(.4, se = se))$S
  expect_equal(S[1, 2], .4 / 4)
  expect_equal(S[1, 12], 0)
  # The NCI form: Var(log(CR r)) = Var(log r) + cr_variance / CR^2 for each
  # old coding observation, no covariance, and nothing with constant variance.
  S <- one(fit(.4, se = se, cr_covariance = "nci"))$S
  r <- d$y
  expect_equal(S[1:6], (d$se[1:6] / r[1:6])^2 + .4 / 4)
  expect_equal(S[7:12], (d$se[7:12] / r[7:12])^2)
  expect_equal(vcov(fit(.4, cr_covariance = "nci")), vcov(fit(0)))
})

test_that("the jump covariance uses the global step", {
  x <- 1:20
  d <- data.frame(x, y = exp(.03 * x + .4 * pmax(x - 14, 0) + .3 * (x > 5) +
                               .05 * sin(x)))
  f <- joinpoint(d, x, y, kmin = 1, kmax = 1, jump = 5,
                 ci_method = "parametric")
  t <- joinpoints(f)$estimate
  g <- as.numeric(x > t)
  X <- cbind(1 - g, (1 - g) * x, g, g * x, as.numeric(x > 5))
  ref <- lm(log(y) ~ X - 1, data = d, subset = x != t)
  expect_equal(vcov(f)[3, 3], vcov(ref)[3, 3], ignore_attr = TRUE)
})

test_that("coefficients and their covariance share names and order", {
  tb <- subset(us_mortality, cause == "Tuberculosis")
  f <- joinpoint(tb, year, rate, kmin = 1, kmax = 1, ci_method = "parametric")
  expect_equal(names(coef(f)), rownames(vcov(f)))
  est <- model_estimates(f)
  expect_equal(unname(coef(f)[c("Intercept 2", "Slope 2")]),
               est$estimate[match(c("Intercept 2", "Slope 2"), est$parameter)])
  expect_named(coef(f, type = "hinge"), c("intercept", "slope", "delta1"))
})

test_that("likelihood fits detect a missing maximum, not small means", {
  d <- data.frame(x = 1:10, n = 10, p = 10)
  expect_warning(f <- joinpoint(d, x, count = n, population = p,
                                family = "binomial", kmax = 0),
                 "No maximum likelihood")
  expect_false(one(f)$models[[1]]$converged)
  d <- data.frame(x = 1:12, n = c(rep(0, 10), 1, 89), p = 1e5)
  f <- suppressWarnings(joinpoint(d, x, count = n, population = p,
                                  family = "poisson", kmax = 0))
  expect_true(one(f)$models[[1]]$converged)
  ref <- suppressWarnings(glm(n ~ x + offset(log(p)), family = poisson(),
                              data = d))
  expect_equal(model_estimates(f)$se, unname(sqrt(diag(vcov(ref)))),
               tolerance = 1e-5)
})

test_that("a rank deficient likelihood design does not win the search", {
  x <- 1:20
  d <- data.frame(x, z = pmax(x - 3, 0), cnt = 100, pop = 1000)
  f <- suppressWarnings(joinpoint(d, x, count = cnt, population = pop,
                                  family = "poisson", kmin = 1, kmax = 1,
                                  method = "bic", covariates = ~ z))
  m <- one(f)$models[[1]]
  X <- cbind(1, x, pmax(x - m$tau, 0), d$z)
  expect_equal(qr(X)$rank, 4L)
  expect_true(isTRUE(m$converged))
})

test_that("an inestimable dispersion is refused, not used", {
  d <- data.frame(x = 1:5, cnt = c(32, 35, 30, 25, 22), pop = 100,
                  z = c(0, 1, 0, 1, 0))
  res <- tryCatch(suppressWarnings(
    joinpoint(d, x, count = cnt, population = pop, family = "poisson",
              kmax = 1, method = "bic", overdispersion = TRUE,
              covariates = ~ z)), error = identity)
  if (inherits(res, "error")) {
    expect_match(conditionMessage(res), "dispersion|degrees of freedom")
  } else {
    disp <- vapply(one(res)$models, `[[`, numeric(1), "mse")
    expect_true(all(is.finite(disp) & disp > 0))
  }
})

test_that("covariance rows follow the data rows", {
  d <- data.frame(x = 1:10, y = 1:10 + sin(1:10), g = c(NA, rep("A", 9)))
  expect_warning(f <- L(d, by = "g", vcov = diag(1:10), kmax = 0),
                 "missing by-variable")
  expect_equal(diag(one(f)$S), 2:10)
  expect_error(L(d[-1, ], vcov = diag(10), kmax = 0), "one row and column")
  expect_error(L(d[-1, ], vcov = diag(9), autocorrelation = .5, kmax = 0),
               "autocorrelation")
})

test_that("AR(1) lags keep the place of excluded observations", {
  f <- L(data.frame(x = 1:10, y = 1:10 + sin(1:10)), kmax = 0, exclude = 5,
         autocorrelation = .5)
  expect_equal(one(f)$S[4, 5], .25)
})

test_that("excluded responses are not validated or transformed", {
  f <- joinpoint(data.frame(x = 1:6, y = c(1, 2, 0, 3, 4, 5)), x, y,
                 exclude = 3, kmax = 0, ci_method = "parametric")
  expect_equal(nobs(f), 5)
  expect_error(joinpoint(data.frame(x = 1:10, y = 1:10), x, y,
                         variance = "glm", kmax = 0))
})

test_that("the known comparability ratio variance is not rescaled", {
  x <- 1:20
  d <- data.frame(x, y = exp(2 + .03 * x + .2 * pmax(x - 10, 0) +
                               .001 * sin(x)) / ifelse(x <= 10, 2, 1))
  d$s <- .1 * d$y
  f <- joinpoint(d, x, y, se = s, jump = 10, comparability_ratio = 2,
                 cr_variance = .4, kmin = 1, kmax = 1,
                 ci_method = "parametric")
  expect_equal(joinpoints(f)$estimate, 10)
  # Every old coding observation shares the log ratio variance .4 / 2^2.
  expect_gte(vcov(f)["Intercept 1", "Intercept 1"], .1)
  # An exact fit leaves no residual variance to combine it with.
  d$y <- exp(2 + .03 * x) / ifelse(x <= 10, 2, 1)
  expect_error(joinpoint(d, x, y, jump = 10, comparability_ratio = 2,
                         cr_variance = .4, kmax = 0, ci_method = "parametric"),
               "degenerate")
})

test_that("estimated autocorrelation respects the gaps of excluded rows", {
  set.seed(20261121)
  x <- 1:300
  e <- as.numeric(stats::arima.sim(list(ar = .8), n = 300, sd = sqrt(1 - .64)))
  d <- data.frame(x, y = 100 + .001 * x + e)
  ex <- seq(1, 300, by = 2)
  fit <- function(dd, ...) {
    one(L(dd, kmax = 0, method = "bic", autocorrelation = "estimate", ...))$S
  }
  S1 <- fit(d, exclude = ex)
  S2 <- fit(d[!d$x %in% ex, ])
  adj <- function(S) S[1, 2] / sqrt(S[1, 1] * S[2, 2])
  expect_lt(abs(adj(S1) - adj(S2)), .05)
})
