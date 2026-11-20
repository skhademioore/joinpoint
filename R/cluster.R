# Multi-group clustering of contiguous age groups.

# Pre-analysis count compression: adjacent age groups are combined until the
# count of every combined group exceeds `min_count` at every time point.
# Forward compression merges sparse groups into the preceding group and
# backward compression into the following group.
jp_compress <- function(C, min_count, direction) {
  A <- ncol(C)
  ord <- if (direction == "forward") seq_len(A) else rev(seq_len(A))
  blocks <- list()
  cur <- integer(0)
  for (a in ord) {
    small <- any(C[, a] <= min_count)
    cur_small <- length(cur) && any(rowSums(C[, cur, drop = FALSE]) <= min_count)
    if (!length(cur) || small || cur_small) {
      cur <- c(cur, a)
    } else {
      blocks[[length(blocks) + 1L]] <- cur
      cur <- a
    }
  }
  blocks[[length(blocks) + 1L]] <- cur
  nb <- length(blocks)
  if (nb > 1L && any(rowSums(C[, blocks[[nb]], drop = FALSE]) <= min_count)) {
    blocks[[nb - 1L]] <- c(blocks[[nb - 1L]], blocks[[nb]])
    blocks[[nb]] <- NULL
  }
  blocks <- lapply(blocks, sort)
  if (direction == "backward") blocks <- rev(blocks)
  blocks
}

# Parallel joinpoint model for a block of compressed age groups: log crude
# (or delay adjusted) rates with group specific intercepts and common
# joinpoints and slopes, weighted by the delta method variance of the log
# rate under Poisson counts.
jp_parallel_fit <- function(x, DC, D2C, P, cols, method, ctrl, kmin, kmax,
                            cand, cons) {
  nt <- length(x)
  xs <- rep(x, length(cols))
  y <- log(as.vector(DC[, cols] / P[, cols]))
  S <- as.vector(D2C[, cols] / DC[, cols]^2)
  extra <- if (length(cols) > 1L) {
    g <- factor(rep(seq_along(cols), each = nt))
    stats::model.matrix(~g)[, -1L, drop = FALSE]
  }
  models <- jp_light_models(xs, y, S, cand, cons, kmin, kmax, extra)
  k <- jp_light_select(xs, y, S, models, method, ctrl, cand, cons, extra,
                       ctrl$perm_seed)
  list(models = models, k = k)
}

#' Multi-group clustering of age specific trends
#'
#' Clusters contiguous age groups whose trends are as parallel as possible
#' and fits a joinpoint model to the age adjusted rates of each cluster:
#' pre-analysis count compression, parallel joinpoint fits
#' for every block of contiguous age groups, the partition into `m`
#' clusters with minimum total residual sum of squares, selection of `m` by
#' BIC, and optional clustering MADWD.
#'
#' Every contiguous block of age groups is searched for joinpoints in the
#' first stage. With `kmax = NULL` the recommended maximum is lowered, with
#' a message, when that search would exceed the time budget of [jp_control()]
#' for all blocks together; an explicit `kmax` is kept, with a message giving
#' the expected time when it exceeds the budget.
#'
#' @inheritParams age_adjust
#' @param data A data frame with one row per age group and time point (and
#'   cohort).
#' @param min_count Age groups are combined until counts exceed this value
#'   at every time point.
#' @param compression `"forward"` or `"backward"` count compression.
#' @param standardize Age groups of the direct standardization of the final
#'   cluster rates: `"original"` (default) standardizes each cluster with its
#'   original age groups, their own rates and standard weights, as
#'   [age_adjust()] does; `"compressed"` uses the compressed age groups, as
#'   the NCI Joinpoint program does. Compression pools the counts and
#'   populations of the merged age groups, so when their population shares
#'   change over time the compressed rates can show a trend that no original
#'   age-specific rate has. Compression always serves the choice of the
#'   clusters.
#' @param method Method for the number of joinpoints within clustering:
#'   `"wbic"`, `"bic"`, `"bic3"` or `"permutation"`.
#' @param max_clusters Maximum number of clusters (default: number of
#'   compressed age groups).
#' @param madwd Clustering minimum APC difference worth detecting, in
#'   percentage points (`NULL` to disable).
#' @param kmin,kmax Minimum and maximum number of joinpoints.
#' @param final_method Model selection method for the final fits of the
#'   cluster age adjusted rates (see [joinpoint()]).
#' @param control Settings from [jp_control()], used for the clustering
#'   fits and the final fits.
#' @param ... Further arguments to [joinpoint()] for the final cluster fits
#'   (for example `ci_method`).
#' @return An object of class `jp_cluster` with the compressed age groups,
#'   the cluster selection table, the final clusters, the MADWD steps, the
#'   cluster rates, `fit`, the [joinpoint()] fit of the cluster age
#'   adjusted rates, and `settings`, the settings of the clustering
#'   (standardization, minimum count, maximum clusters, joinpoint limits and
#'   selection methods).
#' @examples
#' set.seed(1)
#' ages <- c("40-49", "50-59", "60-69", "70-79")
#' d <- expand.grid(year = 2000:2019, age = ages, stringsAsFactors = FALSE)
#' d$pop <- 1e5
#' slope <- ifelse(d$age %in% ages[1:2], 0.03, -0.02)
#' d$count <- rpois(nrow(d), 200 * exp(slope * (d$year - 2000)))
#' cl <- jp_cluster(d, year, count, pop, age, standard = setNames(rep(1, 4), ages),
#'                  kmax = 1, ci_method = "parametric")
#' cl$clusters
#' @export
jp_cluster <- function(data, x, count, population, age, standard, by = NULL,
                       delay = NULL, min_count = 0,
                       compression = c("forward", "backward"),
                       method = c("wbic", "bic", "bic3", "permutation"),
                       max_clusters = NULL, madwd = NULL, kmin = 0,
                       kmax = NULL, final_method = "wbic",
                       rate_per = 100000, control = jp_control(),
                       standardize = c("original", "compressed"), ...) {
  compression <- match.arg(compression)
  method <- match.arg(method)
  standardize <- match.arg(standardize)
  whole <- function(v, lo, name, null_ok = TRUE) {
    if (is.null(v) && null_ok) return(invisible())
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v != round(v) ||
        v < lo) {
      stop("`", name, "` must be a whole number of at least ", lo, ".",
           call. = FALSE)
    }
  }
  whole(max_clusters, 1, "max_clusters")
  whole(kmin, 0, "kmin", FALSE)
  whole(kmax, 0, "kmax")
  if (!is.numeric(min_count) || length(min_count) != 1L ||
      is.na(min_count) || min_count < 0) {
    stop("`min_count` must be a nonnegative number.", call. = FALSE)
  }
  if (!is.null(madwd) && (!is.numeric(madwd) || length(madwd) != 1L ||
                          is.na(madwd) || madwd < 0)) {
    stop("`madwd` must be NULL or a nonnegative number.", call. = FALSE)
  }
  xn <- jp_col(substitute(x), data, "x")
  cn <- jp_col(substitute(count), data, "count")
  pn <- jp_col(substitute(population), data, "population")
  an <- jp_col(substitute(age), data, "age")
  dn <- jp_col(substitute(delay), data, "delay")
  data$.std <- jp_std_rows(substitute(standard), standard, data, an)
  data$.delay <- if (is.null(dn)) 1 else data[[dn]]
  jp_check_age_data(data, xn, an, by, data[[cn]], data[[pn]], data$.std,
                    data$.delay, rate_per)
  # Age order: factor levels, numeric order, or the leading number of each
  # label ("0-4", "5-9", "85+").
  a <- data[[an]]
  ages <- if (is.factor(a)) {
    levels(droplevels(a))
  } else if (is.numeric(a)) {
    as.character(sort(unique(a)))
  } else {
    u <- unique(as.character(a))
    lead <- suppressWarnings(as.numeric(sub("^\\D*(\\d+(\\.\\d+)?).*$", "\\1", u)))
    if (anyNA(lead) || anyDuplicated(lead)) {
      stop("The order of the age groups is ambiguous; give `age` as numbers ",
           "or as a factor with levels in age order.", call. = FALSE)
    }
    u[order(lead)]
  }
  cohorts <- if (is.null(by)) list(seq_len(nrow(data))) else {
    unname(jp_split_keys(data[by]))
  }
  res <- lapply(cohorts, function(i) {
    d <- data[i, , drop = FALSE]
    x <- sort(unique(d[[xn]]))
    pos <- cbind(match(d[[xn]], x), match(as.character(d[[an]]), ages))
    M <- function(v) {
      out <- matrix(NA_real_, length(x), length(ages))
      out[pos] <- v
      out
    }
    C <- M(d[[cn]])
    P <- M(d[[pn]])
    D <- M(d$.delay)
    Sd <- M(d$.std)
    if (anyNA(C)) {
      stop("Every age group of the data needs a row in every cohort.",
           call. = FALSE)
    }
    jp_cluster_one(x, C, P, D, Sd[1L, ], ages, min_count, compression, method,
                   max_clusters, madwd, kmin, kmax, control, rate_per,
                   standardize)
  })
  keys <- if (is.null(by)) NULL else {
    do.call(rbind, lapply(cohorts, function(i) data[i[1L], by, drop = FALSE]))
  }
  tab <- function(name) {
    out <- do.call(rbind, lapply(seq_along(res), function(g) {
      t <- res[[g]][[name]]
      if (is.null(t)) NULL else if (is.null(keys)) t else {
        cbind(keys[rep(g, nrow(t)), , drop = FALSE], t)
      }
    }))
    if (!is.null(out)) rownames(out) <- NULL
    out
  }
  rates <- tab("rates")
  fit <- joinpoint(rates, x = "x", y = "rate", se = "se",
                   by = c(by, "cluster"), kmin = kmin, kmax = kmax,
                   method = final_method, control = control, ...)
  fit$settings$x <- xn
  structure(list(
    compression = tab("compression"), selection = tab("selection"),
    clusters = tab("clusters"), madwd = tab("madwd"), rates = rates,
    fit = fit,
    settings = list(standardize = standardize, min_count = min_count,
                    max_clusters = max_clusters, kmin = kmin, kmax = kmax,
                    method = method, final_method = final_method,
                    rate_per = rate_per)
  ), class = "jp_cluster")
}

jp_cluster_one <- function(x, C, P, D, std, ages, min_count, compression,
                           method, max_clusters, madwd, kmin, kmax, ctrl,
                           rate_per, standardize = "original") {
  blocks <- jp_compress(C, min_count, compression)
  if (length(blocks) == 1L && any(rowSums(C) <= min_count)) {
    stop("Counts are too small even after combining all age groups.",
         call. = FALSE)
  }
  lab <- vapply(blocks, function(b) {
    if (length(b) == 1L) ages[b] else paste0(ages[b[1L]], " to ", ages[max(b)])
  }, character(1))
  # Compressed groups: summed counts, delay adjusted counts (D C), their
  # variance terms (D^2 C), populations and standard populations.
  comb <- function(A) {
    matrix(vapply(blocks, function(b) rowSums(A[, b, drop = FALSE]),
                  numeric(length(x))), length(x))
  }
  Cc <- comb(C)
  DCc <- comb(D * C)
  D2Cc <- comb(D^2 * C)
  Pc <- comb(P)
  Sc <- vapply(blocks, function(b) sum(std[b]), numeric(1))
  G <- length(blocks)
  nt <- length(x)
  cand <- jp_candidates(x, ctrl$n_grid)
  cons <- jp_constraints(x, cand, ctrl$min_obs_end, ctrl$min_obs_between)
  kmax <- jp_cluster_kmax(kmax, kmin, nt, G, cons, method, ctrl)

  fits <- list()
  key <- function(i, j) paste(i, j)
  for (i in seq_len(G)) {
    for (j in i:G) {
      fits[[key(i, j)]] <- jp_parallel_fit(x, DCc, D2Cc, Pc, i:j, method,
                                           ctrl, kmin, kmax, cand, cons)
    }
  }
  best <- function(i, j) {
    k <- max(vapply(c(key(i, j), vapply(i:j, function(l) key(l, l), "")),
                    function(kk) fits[[kk]]$k, numeric(1)))
    ms <- fits[[key(i, j)]]$models
    k <- min(k, max(as.numeric(names(ms))))
    ms[[as.character(k)]]
  }
  rss <- matrix(Inf, G, G)
  kb <- matrix(NA_real_, G, G)
  for (i in seq_len(G)) {
    for (j in i:G) {
      m <- best(i, j)
      rss[i, j] <- m$sse
      kb[i, j] <- m$k
    }
  }
  M <- min(G, if (is.null(max_clusters)) G else max_clusters)
  # Best partition into m contiguous clusters by dynamic programming.
  Q <- matrix(Inf, M, G)
  arg <- matrix(NA_integer_, M, G)
  Q[1L, ] <- rss[1L, ]
  for (m in seq_len(M)[-1L]) {
    for (j in m:G) {
      v <- Q[m - 1L, (m - 1L):(j - 1L)] + rss[m:j, j]
      Q[m, j] <- min(v)
      arg[m, j] <- (m - 1L) + which.min(v) - 1L
    }
  }
  partition <- function(m) {
    cuts <- integer(0)
    j <- G
    for (mm in m:2) {
      i <- arg[mm, j]
      cuts <- c(i, cuts)
      j <- i
    }
    starts <- c(1L, cuts + 1L)
    ends <- c(cuts, G)
    Map(function(s, e) s:e, starts, ends)
  }
  blk_lab <- function(b) {
    a <- c(blocks[[min(b)]][1L], max(blocks[[max(b)]]))
    if (a[1L] == a[2L]) ages[a[1L]] else paste(ages[a[1L]], ages[a[2L]], sep = " to ")
  }
  N <- nt * G
  # Residual sums of squares below the rounding error of the weighted log
  # rates are numerically zero: partitions that all fit exactly are tied,
  # and the penalty chooses the fewest clusters instead of rounding error.
  yw <- log(DCc / Pc) / sqrt(D2Cc / DCc^2)
  yw <- yw[is.finite(yw)]
  qtiny <- max(100 * N * (.Machine$double.eps * max(abs(yw), 0))^2,
               .Machine$double.xmin)
  sel <- do.call(rbind, lapply(seq_len(M), function(m) {
    pt <- if (m == 1L) list(seq_len(G)) else partition(m)
    S <- sum(vapply(pt, function(b) kb[min(b), max(b)], numeric(1)))
    p <- m - 1 + 2 * S + m + G
    data.frame(clusters = m, Q = Q[m, G], joinpoints = S, parameters = p,
               bic = log(max(Q[m, G], qtiny) / N) + p * log(N) / N,
               partition = paste(vapply(pt, blk_lab, ""), collapse = " | "))
  }))
  m_best <- sel$clusters[which.min(sel$bic)]
  sel$selected <- sel$clusters == m_best

  mad <- NULL
  if (!is.null(madwd) && madwd > 0) {
    rows <- list()
    m <- m_best
    repeat {
      pt <- if (m == 1L) list(seq_len(G)) else partition(m)
      dist <- if (m > 1L) {
        vapply(seq_len(m - 1L), function(c) {
          jp_apc_distance(best(min(pt[[c]]), max(pt[[c]])),
                          best(min(pt[[c + 1L]]), max(pt[[c + 1L]])),
                          min(x), max(x))
        }, numeric(1))
      } else {
        numeric(0)
      }
      pass <- m == 1L || all(dist >= madwd)
      rows[[length(rows) + 1L]] <- data.frame(
        clusters = m, min_distance = if (length(dist)) min(dist) else NA,
        result = if (pass) "Pass" else "Fail"
      )
      if (pass) break
      m <- m - 1L
    }
    mad <- do.call(rbind, rows)
    m_best <- m
  }

  final <- if (m_best == 1L) list(seq_len(G)) else partition(m_best)
  cl_lab <- vapply(final, blk_lab, "")
  # The age adjusted rate of each cluster: direct standardization of its
  # original age groups (their own counts, populations, delay factors and
  # standard weights, as age_adjust() computes it) or of its compressed
  # age groups.
  rates <- do.call(rbind, lapply(seq_along(final), function(c) {
    if (standardize == "original") {
      a <- unlist(blocks[final[[c]]])
      w <- std[a] / sum(std[a])
      Pa <- P[, a, drop = FALSE]
      dc <- D[, a, drop = FALSE] * C[, a, drop = FALSE]
      d2c <- D[, a, drop = FALSE]^2 * C[, a, drop = FALSE]
    } else {
      b <- final[[c]]
      w <- Sc[b] / sum(Sc[b])
      Pa <- Pc[, b, drop = FALSE]
      dc <- DCc[, b, drop = FALSE]
      d2c <- D2Cc[, b, drop = FALSE]
    }
    data.frame(
      cluster = cl_lab[c], x = x,
      rate = rate_per * drop((dc / Pa) %*% w),
      se = rate_per * sqrt(drop((d2c / Pa^2) %*% w^2))
    )
  }))
  rates$cluster <- factor(rates$cluster, cl_lab)
  list(
    compression = data.frame(
      age_group = lab,
      min_count = apply(Cc, 2L, min), mean_count = colMeans(Cc),
      max_count = apply(Cc, 2L, max)
    ),
    selection = sel, madwd = mad, rates = rates,
    clusters = data.frame(
      cluster = rep(cl_lab, lengths(final)),
      age_group = lab[unlist(final)]
    )
  )
}

# Number of joinpoints searched for every contiguous block of the first
# stage, under one time budget for all G (G + 1) / 2 blocks together.
jp_cluster_kmax <- function(kmax, kmin, nt, G, cons, method, ctrl) {
  rec <- jp_recommended_kmax(nt)
  cnt <- jp_node_counts(cons, max(rec, kmax, kmin))
  nodes <- c(1, 1 + cumsum(cnt$nodes))
  feasible <- min(sum(cumprod(cnt$models > 0)), max(0, floor((nt - 3) / 2)))
  blocks <- G * (G + 1) / 2
  cost <- function(k) {
    blocks * nodes[k + 1L] * (2.5e-7 + (method == "permutation") * 1.5e-8 *
                                ctrl$n_perm * max(1, k - kmin))
  }
  secs <- function(k) format(signif(cost(k), 2))
  if (is.null(kmax)) {
    k0 <- kmax <- min(rec, feasible)
    while (kmax > kmin && cost(kmax) > ctrl$time_budget) kmax <- kmax - 1
    if (kmax < k0) {
      message(sprintf(paste(
        "Clustering: recommended maximum of %d joinpoints lowered to %d to",
        "keep the search of the %d age group blocks within about %s",
        "seconds; set kmax to override."
      ), k0, kmax, blocks, format(signif(ctrl$time_budget, 2))))
    }
  } else if (kmax <= length(cnt$nodes) && cost(kmax) > ctrl$time_budget) {
    message(sprintf(paste(
      "Clustering: searching %d joinpoints in each of %d age group blocks",
      "is expected to take about %s seconds."
    ), kmax, blocks, secs(kmax)))
  }
  if (kmin > kmax) {
    stop("`kmin` exceeds the largest feasible number of joinpoints (", kmax,
         ").", call. = FALSE)
  }
  kmax
}

# Duration weighted mean absolute APC difference of two fits over [a, b].
jp_apc_distance <- function(m1, m2, a, b) {
  br <- sort(unique(c(a, m1$tau, m2$tau, b)))
  mids <- (utils::head(br, -1L) + utils::tail(br, -1L)) / 2
  len <- diff(br)
  apc_at <- function(m, t) 100 * expm1(m$slopes[findInterval(t, m$tau) + 1L])
  sum(len * abs(apc_at(m1, mids) - apc_at(m2, mids))) / (b - a)
}

#' @export
print.jp_cluster <- function(x, ...) {
  cat("Multi-group clustering of contiguous age groups\n\n")
  cat("Compressed age groups:\n")
  print(x$compression, digits = 4, row.names = FALSE)
  cat("\nCluster selection (BIC):\n")
  print(x$selection, digits = 4, row.names = FALSE)
  if (!is.null(x$madwd)) {
    cat("\nClustering MADWD:\n")
    print(x$madwd, digits = 4, row.names = FALSE)
  }
  cat("\nFinal clusters:\n")
  print(x$clusters, row.names = FALSE)
  cat("\nJoinpoint fits of the cluster age adjusted rates:\n")
  print(x$fit)
  invisible(x)
}
