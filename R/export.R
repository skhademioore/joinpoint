#' Export joinpoint results
#'
#' Writes the result tables of a joinpoint fit to delimited text files, one
#' file per table: data and modeled values, APCs, AAPCs, joinpoint
#' locations, model estimates, model selection (criteria, permutation tests,
#' data dependent selection, MADWD and the final selected models with their
#' search limits), coefficient covariances, jump estimates, pairwise and
#' AAPC comparisons, clustering tables, skipped cohorts, the centered
#' covariate design, a supplied covariance matrix and every setting of the
#' run. `saveRDS()` keeps the complete fitted object.
#'
#' @param object A `joinpoint` or `jp_cluster` object.
#' @param dir Output directory (created if needed).
#' @param prefix File name prefix.
#' @param model `"all"` (every fitted model) or `"final"`.
#' @param comparison Optional result of [jp_compare()] to export with the
#'   fit.
#' @param sep Field separator.
#' @return The paths of the written files, invisibly.
#' @examples
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
#' dir <- file.path(tempdir(), "joinpoint-example")
#' files <- jp_export(fit, dir)
#' basename(files)
#' unlink(dir, recursive = TRUE)
#' @export
jp_export <- function(object, dir, prefix = "joinpoint",
                      model = c("all", "final"), comparison = NULL,
                      sep = ",") {
  model <- match.arg(model)
  cluster <- NULL
  if (inherits(object, "jp_cluster")) {
    cluster <- object
    object <- object$fit
  }
  stopifnot(inherits(object, "joinpoint"))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ms <- model_selection(object)
  s <- object$settings
  flat <- c(list(by = object$by), s[setdiff(names(s), c("control", "vcov"))],
            s$control,
            list(package_version = as.character(utils::packageVersion("joinpoint"))))
  settings <- data.frame(
    setting = names(flat),
    value = vapply(flat, function(v) {
      if (inherits(v, "formula")) {
        paste(deparse(v), collapse = " ")
      } else if (is.null(v)) "" else paste(format(unlist(v)), collapse = " ")
    }, character(1))
  )
  covariance <- jp_collect(object, function(f, m) {
    ms <- if (is.null(m)) f$models else list(m)
    do.call(rbind, lapply(ms, function(mm) {
      V <- jp_vcov(mm)
      cbind(jp_model_cols(f, mm),
            data.frame(row = rownames(V)[row(V)], column = colnames(V)[col(V)],
                       covariance = as.vector(V)))
    }))
  }, model)
  tables <- list(
    Data = as.data.frame(object, model = model),
    APCs = apc(object, model),
    AAPCs = aapc(object, model),
    Joinpoints = joinpoints(object, model),
    Model.Estimates = model_estimates(object, model),
    Coefficient.Covariance = covariance,
    Model.Selection = ms$criteria,
    Permutation.Tests = ms$tests,
    Selection.Results = ms$dds,
    MADWD = ms$madwd,
    Final.Selected.Model = ms$selected,
    Jump.Model = jump_estimates(object, model),
    Comparison = if (!is.null(comparison)) jp_plain(comparison),
    AAPC.Comparison = attr(comparison, "aapc"),
    Final.Cluster = cluster$clusters,
    Cluster.Summary = cluster$selection,
    Cluster.MADWD = cluster$madwd,
    Skipped.Cohorts = object$skipped,
    Covariate.Design = jp_collect(object, function(f, m) {
      if (!is.null(f$covariates)) {
        data.frame(x = f$data$x, f$covariates, check.names = FALSE)
      }
    }),
    Supplied.Covariance = if (!is.null(s$vcov)) as.data.frame(s$vcov),
    Run.Summary = settings
  )
  tables <- tables[!vapply(tables, is.null, logical(1))]
  paths <- file.path(dir, paste0(prefix, ".", names(tables), ".txt"))
  for (i in seq_along(tables)) {
    utils::write.table(tables[[i]], paths[i], sep = sep, row.names = FALSE,
                       quote = TRUE)
  }
  invisible(paths)
}
