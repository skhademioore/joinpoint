test_that("Poisson and binomial joinpoints are maximum likelihood fits", {
  set.seed(11)
  x <- 1:30
  pop <- rep(1e5, 30)
  eta <- log(50 / 1e5) + 0.03 * x - 0.08 * pmax(x - 15, 0) -
    0.03 * pmax(x - 23, 0)
  d <- data.frame(x, cnt = rpois(30, exp(eta) * pop), pop)
  f <- joinpoint(d, x, count = cnt, population = pop, family = "poisson",
                 kmax = 2, method = "bic")
  best <- Inf
  for (a in 3:28) for (b in 3:28) if (b - a >= 3) {
    X <- cbind(1, x, pmax(x - a, 0), pmax(x - b, 0))
    dv <- glm.fit(X, d$cnt, offset = log(pop), family = poisson())$deviance
    if (dv < best) {
      best <- dv
      bt <- c(a, b)
    }
  }
  m2 <- f$fits[[1]]$models[["2"]]
  expect_equal(m2$tau, bt)
  expect_equal(m2$sse, best, tolerance = 1e-8)
  a <- apc(f)
  expect_true(all(a$lower < a$estimate & a$estimate < a$upper))

  p <- plogis(-2 + 0.05 * x - 0.1 * pmax(x - 18, 0))
  db <- data.frame(x, ev = rbinom(30, 500, p), n = 500)
  fb <- joinpoint(db, x, count = ev, population = n, family = "binomial",
                  kmax = 1, kmin = 1, method = "bic")
  expect_equal(fb$fits[[1]]$models[["1"]]$tau, 19)
  expect_match(names(jp_table(fb))[3], "odds")
  expect_s3_class(autoplot(fb), "ggplot")
  expect_error(joinpoint(db, x, count = ev, population = n,
                         family = "binomial", ci_method = "empirical"))
})

test_that("Monte Carlo likelihood ratio test selects the number of joinpoints", {
  set.seed(3)
  x <- 1:25
  d <- data.frame(x, cnt = rpois(25, exp(4 + 0.06 * x - 0.15 * pmax(x - 12, 0))),
                  pop = 1e5)
  f <- joinpoint(d, x, count = cnt, population = pop, family = "poisson",
                 kmax = 2, method = "lrt", control = jp_control(n_perm = 99))
  tests <- model_selection(f)$tests
  expect_true(nrow(tests) >= 1)
  expect_equal(f$fits[[1]]$k_selected, 1)
})

test_that("continuous joinpoints improve on the grid", {
  set.seed(2)
  x <- 1:25
  y <- exp(2 + 0.04 * x - 0.09 * pmax(x - 12.4, 0) + rnorm(25, sd = 0.01))
  d <- data.frame(x, y)
  fc <- joinpoint(d, x, y, kmin = 1, kmax = 1, continuous = TRUE,
                  ci_method = "parametric")
  fg <- joinpoint(d, x, y, kmin = 1, kmax = 1, ci_method = "parametric",
                  control = jp_control(n_grid = 9))
  mc <- fc$fits[[1]]$models[["1"]]
  expect_lte(mc$sse, fg$fits[[1]]$models[["1"]]$sse)
  expect_equal(mc$tau, 12.4, tolerance = 0.2)
  expect_false(mc$tau %in% x)
})

test_that("covariates adjust the trend", {
  set.seed(5)
  t <- 1:72
  mo <- rep(1:12, 6)
  y <- exp(1 + 0.01 * t - 0.03 * pmax(t - 40, 0) +
             0.3 * sin(2 * pi * mo / 12) + rnorm(72, sd = 0.01))
  d <- data.frame(t, mo, y)
  f <- joinpoint(d, t, y, kmin = 1, kmax = 1, covariates = ~ factor(mo),
                 ci_method = "parametric")
  expect_equal(joinpoints(f)$estimate, 40)
  expect_equal(apc(f)$estimate, 100 * expm1(c(0.01, -0.02)), tolerance = 0.05)
  expect_gt(model_selection(f)$criteria$r_squared, 0.99)
})

test_that("the likelihood search finds the maximum where reweighting failed", {
  d <- data.frame(x = 1:15, pop = 100,
                  cnt = c(93, 4, 6, 49, 6, 23, 22, 9, 2, 11, 90, 4, 197, 3, 0))
  f <- joinpoint(d, x, count = cnt, population = pop, family = "poisson",
                 kmin = 1, kmax = 1, method = "bic")
  m <- f$fits[[1]]$models[["1"]]
  dev <- sapply(3:13, function(a) {
    glm.fit(cbind(1, d$x, pmax(d$x - a, 0)), d$cnt, offset = log(d$pop),
            family = poisson())$deviance
  })
  expect_equal(m$tau, 13)
  expect_equal(m$sse, min(dev), tolerance = 1e-8)
  # The joinpoint interval holds the candidates within the chi-square
  # cutoff of the profile deviance.
  j <- joinpoints(f)
  ok <- (3:13)[dev - min(dev) <= qchisq(0.95, 1)]
  expect_equal(c(j$lower, j$upper), range(ok))
})

test_that("likelihood fits without estimates withhold inference", {
  d <- data.frame(x = 1:20, ev = rep(c(0, 100), each = 10), n = 100)
  expect_warning(f <- joinpoint(d, x, count = ev, population = n,
                                family = "binomial", kmax = 0, method = "bic"),
                 "No maximum likelihood estimate")
  expect_true(all(is.na(model_estimates(f)$se)))
  expect_error(joinpoint(d, x, count = ev, population = n,
                         family = "binomial", method = "lrt",
                         overdispersion = TRUE), "overdispersion")
  expect_error(joinpoint(transform(d, y = 0.5), x, y, count = ev,
                         population = n, family = "binomial"), "no `y`")
  expect_error(joinpoint(transform(d, ev = ev + 0.5), x, count = ev,
                         population = n, family = "binomial"), "whole-number")
})

test_that("continuous joinpoints are the global least squares solution", {
  for (s in 1:5) {
    set.seed(s)
    x <- 1:20
    y <- 20 + 0.3 * x + 1.5 * pmax(x - 6.4, 0) - 2.2 * pmax(x - 13.7, 0) +
      rnorm(20)
    d <- data.frame(x, y)
    fc <- joinpoint(d, x, y, log = FALSE, kmin = 2, kmax = 2,
                    continuous = TRUE, method = "bic")
    fg <- joinpoint(d, x, y, log = FALSE, kmin = 2, kmax = 2, method = "bic",
                    ci_method = "parametric", control = jp_control(n_grid = 9))
    expect_lte(fc$fits[[1]]$models[["2"]]$sse,
               fg$fits[[1]]$models[["2"]]$sse + 1e-8)
  }
  y <- 1 + 0.5 * x + 2 * pmax(x - 10.4, 0)
  f <- joinpoint(data.frame(x, y), x, y, log = FALSE, kmin = 1, kmax = 1,
                 continuous = TRUE, method = "bic")
  expect_equal(unlist(joinpoints(f)[c("estimate", "lower", "upper")]),
               c(estimate = 10.4, lower = 10.4, upper = 10.4))
  expect_error(joinpoint(data.frame(x, y), x, y, log = FALSE,
                         continuous = TRUE, ci_method = "empirical"),
               "parametric intervals")
})

test_that("covariates enter predictions, residuals and degrees of freedom", {
  set.seed(5)
  d <- data.frame(x = 1:20, z = rnorm(20))
  d$y <- 1 + 0.2 * d$x + 1.5 * d$z + rnorm(20, sd = 0.1)
  f <- joinpoint(d, x, y, type = "other", log = FALSE, kmax = 0,
                 covariates = ~ z, ci_method = "parametric")
  l <- lm(y ~ x + z, d)
  expect_equal(fitted(f), unname(fitted(l)))
  expect_equal(residuals(f), unname(resid(l)))
  expect_equal(model_selection(f)$criteria$df, 17)
  est <- model_estimates(f)
  expect_equal(est$se[est$parameter == "z"],
               summary(l)$coefficients["z", 2])
  p <- predict(f, newdata = data.frame(x = 5, z = c(0, 1)))
  expect_equal(diff(p$fit), unname(coef(l)["z"]))
  expect_error(predict(f, newdata = data.frame(x = 5)), "z")
})

test_that("covariates named like the model's coefficients are refused", {
  set.seed(5)
  d <- data.frame(x = 1:20, slope = rnorm(20))
  d$y <- 1 + 0.2 * d$x + 1.5 * d$slope + rnorm(20, sd = 0.1)
  expect_error(joinpoint(d, x, y, type = "other", log = FALSE, kmax = 0,
                         covariates = ~ slope, ci_method = "parametric"),
               "rename slope")
  names(d)[2L] <- "delta1"
  expect_error(joinpoint(d, x, y, type = "other", log = FALSE, kmax = 1,
                         covariates = ~ delta1, ci_method = "parametric"),
               "rename delta1")
})

test_that("continuous joinpoints beat every admissible fine grid placement", {
  # Weighted two-joinpoint series on which a local refinement failed; `fine`
  # holds an admissible placement found on a grid of ten points per gap.
  cases <- readRDS(test_path("fixtures", "continuous_cases.rds"))
  for (d in cases) {
    cand <- joinpoint:::jp_candidates(d$x, 1)
    cons <- joinpoint:::jp_constraints(d$x, cand, closed = TRUE)
    cs <- joinpoint:::jp_continuous_search(d$x, d$y, d$S, cand, cons, 2)
    sse <- function(tau) joinpoint:::jp_fit_fixed(d$x, d$y, d$S, tau)$sse
    expect_lte(sse(cs$tau), sse(d$fine) + 1e-8)
  }
})

test_that("continuous joinpoint intervals end where the profile crosses the cutoff", {
  set.seed(4)
  x <- 1:20
  y <- 1 + 0.5 * x + 2 * pmax(x - 10.4, 0) + rnorm(20, sd = 0.3)
  f <- joinpoint(data.frame(x, y), x, y, log = FALSE, kmin = 1, kmax = 1,
                 continuous = TRUE, method = "bic")
  m <- f$fits[[1]]$models[["1"]]
  crit <- joinpoint:::jp_tau_crit(m, 0.05)
  ts <- seq(9, 12, by = 0.001)
  q <- sapply(ts, function(t) sum(resid(lm(y ~ x + pmax(x - t, 0)))^2))
  expect_equal(unlist(joinpoints(f)[c("lower", "upper")]),
               c(lower = min(ts[q <= crit]), upper = max(ts[q <= crit])),
               tolerance = 1e-3)
})
