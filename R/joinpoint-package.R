#' joinpoint: Trend Analysis by Joinpoint Regression
#'
#' Joinpoint regression for trend analysis of rates, proportions, counts
#' and other series (Kim et al. 2000). The main entry point is [joinpoint()]; results are summarized with [apc()],
#' [aapc()], [model_selection()] and plotted with [autoplot.joinpoint()].
#'
#' @keywords internal
#' @useDynLib joinpoint, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom ggplot2 .data
#' @importFrom stats coef fitted residuals predict vcov nobs
"_PACKAGE"

#' @export
ggplot2::autoplot
