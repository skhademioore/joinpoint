test_that("AAPC weights slopes by segment length (Clegg et al. 2009)", {
  bounds <- c(1975, 1988, 1992, 1995, 2003)
  b <- c(0.026, 0.151, -0.113, 0.007)
  aapc <- function(a, e) {
    100 * expm1(sum(joinpoint:::jp_range_weights(bounds, a, e) * b))
  }
  expect_equal(round(aapc(1994, 2003), 1), -0.6)
  expect_equal(round(aapc(1999, 2003), 1), 0.7)
  expect_equal(round(aapc(1975, 2003), 1), 2.4)
})

test_that("empirical quantile p value is twice the smaller tail share", {
  expect_equal(joinpoint:::jp_eq_pvalue(c(-1, 1, 2, 3)), 0.5)
  expect_equal(joinpoint:::jp_eq_pvalue(1:10), 0)
})

test_that("AAPC of a range inside one segment equals the segment APC", {
  tb <- subset(us_mortality, cause == "Tuberculosis")
  fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric",
                   control = jp_control(aapc_ranges = list(c(1960, 1970))))
  a <- apc(fit)
  r <- aapc(fit)
  inside <- r[r$start == 1960, ]
  seg <- a[a$start <= 1960 & a$end >= 1970, ]
  expect_true(inside$within_segment)
  expect_equal(inside$estimate, seg$estimate)
  expect_equal(c(inside$lower, inside$upper), c(seg$lower, seg$upper))
})

test_that("published joinpoint analysis of prepregnancy BMI is reproduced", {
  mb <- unique(ontario_bmi[c("year", "births", "mean_bmi", "sd_bmi")])
  mb$se <- mb$sd_bmi / sqrt(mb$births)
  fit <- joinpoint(mb, year, mean_bmi, se = se, log = FALSE, kmax = 2,
                   method = "bic", ci_method = "parametric")
  # Khademioore et al. (2026), Table 4: BIC 3.73, 2.93 and 3.33; one
  # joinpoint at 2016; slope change 0.109 (p = 0.010).
  crit <- model_selection(fit)$criteria
  expect_equal(round(crit$bic, 2), c(3.73, 2.93, 3.33))
  expect_equal(joinpoints(fit)$estimate, 2016)
  est <- model_estimates(fit)
  change <- est[est$parameter == "Slope 2 - Slope 1", ]
  expect_equal(round(change$estimate, 3), 0.109)
  expect_equal(round(change$p_value, 2), 0.01)
  expect_equal(round(apc(fit)$estimate, 3), c(0.016, 0.125))
})

test_that("excluded points shift last-N AAPC ranges", {
  x <- 2000:2022
  r <- jp_aapc_ranges <- joinpoint:::jp_aapc_ranges
  inc <- x != 2020
  expect_equal(unname(r(x, inc, FALSE, 5)[[1]]), c(2018, 2022))
  inc <- x != 2022
  expect_equal(unname(r(x, inc, FALSE, 5)[[1]]), c(2017, 2021))
})

test_that("an AAPC inside an estimable segment keeps its interval", {
  x <- 1:12
  set.seed(2)
  y <- exp(1 + 0.05 * x - 0.1 * pmax(x - 4, 0) + 0.2 * pmax(x - 5, 0) +
             rnorm(12, sd = 0.01))
  f <- joinpoint(data.frame(x, y), x, y, kmin = 2, kmax = 2,
                 ci_method = "parametric",
                 control = jp_control(min_obs_between = 0))
  a <- aapc(f, ranges = list(c(6, 12)))
  expect_false(anyNA(a[c("lower", "upper", "p_value")]))
  expect_equal(a$ci_method, "parametric")
  expect_error(aapc(f, ranges = list(c(10, 3))), "aapc_ranges")
  expect_equal(aapc(f, ranges = list(c(-100, 100)))$estimate,
               aapc(f)$estimate)
})

test_that("published segment tests and joinpoint interval are reproduced", {
  mb <- unique(ontario_bmi[c("year", "births", "mean_bmi", "sd_bmi")])
  mb$se <- mb$sd_bmi / sqrt(mb$births)
  args <- list(mb, quote(year), quote(mean_bmi), se = quote(se), log = FALSE,
               kmax = 2, method = "bic", ci_method = "parametric")
  fit <- do.call(joinpoint, args)
  # Khademioore et al. (2026), Table 4: segment p values 0.49 and < 0.001;
  # joinpoint interval 2014-2019, reproduced by the conservative variant.
  a <- apc(fit)
  expect_equal(round(a$p_value[1], 2), 0.49)
  expect_lt(a$p_value[2], 0.001)
  expect_equal(unlist(joinpoints(fit)[c("lower", "upper")]),
               c(lower = 2014, upper = 2017))
  lit <- do.call(joinpoint, c(args, list(
    control = jp_control(tau_interval = "literal"))))
  expect_equal(unlist(joinpoints(lit)[c("lower", "upper")]),
               c(lower = 2014, upper = 2019))
})
