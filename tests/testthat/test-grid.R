test_that("grid search reproduces the published tutorial fits", {
  fits <- utils::read.csv(test_path("fixtures", "tutorial_model_fits.csv"))
  key <- c(Accidents = "accidents", Cancer = "cancer",
           `Heart Disease` = "heart", `Influenza and Pneumonia` = "pi",
           Stroke = "stroke", Tuberculosis = "tb")
  expected_jp <- list(
    accidents = list(1904, c(1906, 1920), c(1906, 1921, 1967)),
    cancer = list(1928, c(1928, 1993), c(1928, 1976, 1991)),
    heart = list(1954, c(1943, 1962), c(1920, 1937, 1961)),
    pi = list(1955, c(1918, 1949), c(1914, 1918, 1948)),
    stroke = list(1967, c(1947, 1959), c(1924, 1939, 1960)),
    tb = list(1951, c(1932, 1958), c(1918, 1922, 1957))
  )
  # The published fitted values have two decimals. Eight cells are
  # transcription errors: each breaks the piecewise linear shape of its own
  # column (it differs from the mean of its neighbors in the same segment by
  # 0.1 to 2). Eight more differ from the package by 0.008 to 0.036 and are
  # inconsistent with their neighbors at the second decimal, or are the
  # last value of a column. Every other cell agrees to display precision.
  typos <- list(accidents.3 = 1923, cancer.1 = 1974, cancer.2 = 1927,
                cancer.3 = 1970, heart.0 = 1917, pi.1 = 1918,
                stroke.0 = c(1978, 1998))
  near <- list(accidents.0 = 1907, accidents.3 = 1968, heart.1 = 1998,
               heart.2 = c(1960, 1986), pi.0 = 1936, stroke.1 = 1966,
               stroke.2 = 1988)
  for (cause in names(key)) {
    d <- us_mortality[us_mortality$cause == cause, ]
    fit <- joinpoint(d, year, rate, log = FALSE, kmin = 0, kmax = 3,
                     method = "bic", ci_method = "parametric")
    for (k in 0:3) {
      m <- fit$fits[[1]]$models[[as.character(k)]]
      if (k) expect_equal(m$tau, expected_jp[[key[[cause]]]][[k]])
      modeled <- fit$fits[[1]]$fitted[[as.character(k)]]
      ref <- fits[[paste0(key[[cause]], ".", k)]]
      diff <- abs(modeled - ref)
      col <- paste0(key[[cause]], ".", k)
      typo <- d$year %in% typos[[col]]
      off <- d$year %in% near[[col]]
      expect_true(all(diff[off] > 0.0051 & diff[off] < 0.04))
      expect_true(all(diff[!typo & !off] <= 0.0051))
    }
  }
})

test_that("compiled search matches brute force weighted least squares", {
  set.seed(1)
  x <- 1:20
  y <- c(1:8, 8 - 0.5 * (1:12)) + rnorm(20, sd = 0.3)
  w <- runif(20, 0.5, 2)
  cand <- joinpoint:::jp_candidates(x, 1)
  cons <- joinpoint:::jp_constraints(x, cand, 2, 2)
  g <- joinpoint:::jp_grid_search(x, y, w, cand, cons, 2)
  sse <- function(tau) {
    X <- cbind(1, x, outer(x, tau, function(a, b) pmax(a - b, 0)))
    sum(w * lm.wfit(X, y, w)$residuals^2)
  }
  ok <- function(tau) {
    nb <- sum(x < tau[1]); na <- sum(x > tau[length(tau)])
    btw <- if (length(tau) > 1) sum(x > tau[1] & x < tau[2]) else Inf
    nb >= 2 && na >= 2 && btw >= 2
  }
  best1 <- min(vapply(cand[vapply(cand, ok, TRUE)], sse, 1))
  pairs <- utils::combn(cand, 2, simplify = FALSE)
  best2 <- min(vapply(pairs[vapply(pairs, ok, TRUE)], sse, 1))
  expect_equal(g$sse[1, 2], best1, tolerance = 1e-8)
  expect_equal(g$sse[1, 3], best2, tolerance = 1e-8)
  expect_equal(g$sse[1, 1], sse(numeric(0)), tolerance = 1e-8)
})

test_that("placement constraints count observations off the joinpoints", {
  x <- 1990:2014
  cand <- joinpoint:::jp_candidates(x, 0)
  cons <- joinpoint:::jp_constraints(x, cand, 2, 2)
  expect_equal(cand[which(cons$first_ok == 1)[1]], 1992)
  i2000 <- which(cand == 2000)
  expect_equal(cand[cons$next_min[i2000] + 1], 2003)
  cand3 <- joinpoint:::jp_candidates(x, 3)
  expect_true(all(c(1995.25, 2000.5, 2010.75) %in% cand3))
})

test_that("permutation test reproduces the published tutorial selections", {
  skip_on_cran()
  # van Doren (2024), linear models with 4499 permutations: 1 joinpoint for
  # influenza and pneumonia, 3 for every other cause.
  fit <- joinpoint(us_mortality, year, rate, by = "cause", log = FALSE,
                   kmax = 3, method = "permutation", ci_method = "parametric")
  k <- model_selection(fit)$selected$joinpoints
  expect_equal(k, c(3, 3, 3, 1, 3, 3))
})
