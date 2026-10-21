# ggplot2 graphics.

jp_pal <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
            "#4a3aa7", "#e34948")
jp_ink <- c(primary = "#0b0b0b", secondary = "#52514e", muted = "#898781",
            grid = "#e1e0d9", axis = "#c3c2b7")
jp_trend_cols <- c(Increasing = "#d03b3b", Decreasing = "#2a78d6",
                   Stable = "#898781", `Not tested` = "#c3c2b7")

#' Theme and color scales for joinpoint graphics
#'
#' `theme_joinpoint()` is a clean publication theme with hairline horizontal
#' gridlines and recessive axes. `scale_color_joinpoint()` and
#' `scale_fill_joinpoint()` use an eight color qualitative palette whose
#' order keeps adjacent colors distinguishable under color vision
#' deficiency.
#'
#' @param base_size Base font size in points.
#' @param base_family Base font family.
#' @param ... Passed to [ggplot2::discrete_scale()].
#' @param aesthetics Aesthetics the color scale applies to.
#' @return A ggplot2 theme or scale.
#' @examples
#' library(ggplot2)
#' ggplot(us_mortality, aes(year, rate, color = cause)) +
#'   geom_line() +
#'   scale_color_joinpoint() +
#'   theme_joinpoint()
#' @export
theme_joinpoint <- function(base_size = 11, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = jp_ink[["grid"]],
                                                 linewidth = 0.3),
      axis.line.x = ggplot2::element_line(color = jp_ink[["axis"]],
                                          linewidth = 0.4),
      axis.ticks.x = ggplot2::element_line(color = jp_ink[["axis"]],
                                           linewidth = 0.4),
      axis.ticks.length = ggplot2::unit(3, "pt"),
      axis.text = ggplot2::element_text(color = jp_ink[["secondary"]]),
      axis.title = ggplot2::element_text(color = jp_ink[["primary"]]),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0,
                                         color = jp_ink[["primary"]]),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.title.position = "plot",
      plot.caption = ggplot2::element_text(color = jp_ink[["muted"]],
                                           hjust = 0, size = base_size * 0.8),
      plot.caption.position = "plot",
      legend.position = "bottom",
      legend.title = ggplot2::element_text(color = jp_ink[["secondary"]]),
      panel.spacing = ggplot2::unit(1.2, "lines")
    )
}

#' @rdname theme_joinpoint
#' @export
scale_color_joinpoint <- function(..., aesthetics = "colour") {
  ggplot2::discrete_scale(aesthetics, palette = jp_pal_fun, ...)
}

#' @rdname theme_joinpoint
#' @export
scale_fill_joinpoint <- function(...) {
  ggplot2::discrete_scale("fill", palette = jp_pal_fun, ...)
}

jp_pal_fun <- function(n) {
  if (n > length(jp_pal)) {
    stop("At most ", length(jp_pal), " colors; facet instead.", call. = FALSE)
  }
  jp_pal[seq_len(n)]
}

# Fitted curve on a fine grid, split into pieces at joinpoints and at a
# jump; `segment` is the joinpoint segment of each piece.
jp_curve <- function(f, m, band = TRUE, n = 120L) {
  inc <- f$data$x[f$data$included]
  brk <- sort(unique(c(min(inc), m$tau, f$jump_at, max(inc))))
  cv <- do.call(rbind, lapply(seq_len(length(brk) - 1L), function(i) {
    xs <- seq(brk[i], brk[i + 1L], length.out = n)
    # The curve before a jump ends just short of it.
    if (identical(brk[i + 1L], f$jump_at)) xs[n] <- xs[n] - 1e-9 * diff(range(inc))
    data.frame(x = xs, piece = i,
               segment = findInterval((brk[i] + brk[i + 1L]) / 2, m$tau) + 1L)
  }))
  cv$fit <- jp_fitted_at(f, m, cv$x)
  cv$lower <- cv$upper <- NA_real_
  if (band && !is.null(m$eq)) {
    q <- jp_eq_band(f, m, cv$x)
    cv$lower <- q[, 1L]
    cv$upper <- q[, 2L]
  }
  cv
}

# Pointwise percentile band of the fitted curve from the empirical quantile
# resamples.
jp_eq_band <- function(f, m, xs, alpha = 0.05) {
  eq <- m$eq
  ex <- if (!is.null(m$jump)) cbind(as.numeric(xs >= m$jump))
  vals <- vapply(seq_len(nrow(eq$coef)), function(b) {
    jp_eval(eq$coef[b, ], if (m$k) eq$tau[b, ] else numeric(0), xs, ex)
  }, numeric(length(xs)))
  q <- t(apply(vals, 1L, stats::quantile, c(alpha / 2, 1 - alpha / 2),
               na.rm = TRUE, names = FALSE))
  q <- jp_linkinv(m, q)
  cr <- f$comparability_ratio
  if (!is.null(cr)) q[xs < f$jump_at, ] <- q[xs < f$jump_at, ] / cr
  q
}

# Fixed decimals; adding zero turns a rounded -0 into 0.
jp_fmt <- function(v, digits = 1) {
  formatC(round(v, digits) + 0, digits = digits, format = "f")
}

# Direction of each segment: significant increase or decrease, stable, or
# not tested when no p value is available.
jp_trend <- function(estimate, p, alpha) {
  factor(ifelse(is.na(p), "Not tested",
                ifelse(p >= alpha, "Stable",
                       ifelse(estimate > 0, "Increasing", "Decreasing"))),
         names(jp_trend_cols))
}

# Assembles plotting data for the selected cohorts and model.
jp_plot_data <- function(object, model, cohorts = NULL, digits = 1,
                         band = FALSE) {
  idx <- seq_along(object$fits)
  if (!is.null(cohorts)) idx <- idx[cohorts]
  idx <- idx[!vapply(object$fits[idx], is.null, logical(1))]
  labs <- vapply(idx, function(g) jp_cohort_label(object, g), character(1))
  alpha <- object$settings$control$alpha_apc
  out <- lapply(seq_along(idx), function(i) {
    f <- object$fits[[idx[i]]]
    m <- jp_pick(f, model)
    a <- m$apc
    mid <- (a$start + a$end) / 2
    list(
      points = data.frame(cohort = labs[i], x = f$data$x, y = f$data$r,
                          se = f$data$se, included = f$data$included),
      curve = data.frame(cohort = labs[i], jp_curve(f, m, band)),
      joinpoints = if (m$k) {
        data.frame(cohort = labs[i], x = m$tau, y = jp_fitted_at(f, m, m$tau),
                   lower = m$joinpoints$lower, upper = m$joinpoints$upper)
      },
      segments = data.frame(
        cohort = labs[i], segment = a$segment, start = a$start, end = a$end,
        x = mid, y = jp_fitted_at(f, m, mid), estimate = a$estimate,
        lower = a$lower, upper = a$upper,
        trend = jp_trend(a$estimate, a$p_value, alpha),
        label = paste0(jp_fmt(a$estimate, digits),
                       ifelse(a$significant, "*", ""))
      ),
      jump = if (!is.null(f$jump_at)) {
        cr <- f$comparability_ratio %||% exp(m$coef[["jump"]])
        data.frame(cohort = labs[i], x = f$jump_at,
                   label = paste0("CR ", formatC(cr, digits = 3, format = "f")))
      },
      aapc = if (!is.null(m$aapc)) {
        data.frame(cohort = labs[i], m$aapc[1L, , drop = FALSE])
      }
    )
  })
  part <- function(nm) {
    d <- do.call(rbind, lapply(out, `[[`, nm))
    if (!is.null(d)) d$cohort <- factor(d$cohort, unique(labs))
    d
  }
  stats::setNames(lapply(names(out[[1L]]), part), names(out[[1L]]))
}

#' Plot a joinpoint fit
#'
#' Publication quality graph of observed values and the fitted joinpoint
#' model: observed points (excluded points hollow), the fitted curve,
#' joinpoints, and the APC of each segment as a direct label (`*` marks
#' APCs significantly different from zero). A coding change is marked by a
#' dotted line with its comparability ratio.
#'
#' @param object A `joinpoint` object.
#' @param model `"final"` or a number of joinpoints.
#' @param color_by `"cohort"` colors each cohort; `"trend"` colors each
#'   segment by its direction: significantly increasing, significantly
#'   decreasing, stable, or not tested.
#' @param layout `"facet"` draws one panel per cohort; `"overlay"` draws all
#'   cohorts in one panel.
#' @param cohorts Optional integer indices of cohorts to draw.
#' @param labels Segment labels: `"apc"` (APC, or slope for linear models),
#'   `"ci"` (APC with its confidence interval) or `"none"`. The default is
#'   `"apc"`, or `"none"` when several cohorts are overlaid.
#' @param aapc Add the full range AAPC to the title (one cohort), to each
#'   panel title (facets) or to each legend entry (overlaid cohorts).
#' @param band Draw the pointwise 95% band of the fitted curve from the
#'   empirical quantile resamples, when available.
#' @param joinpoint_ci Shade the confidence interval of each joinpoint
#'   location.
#' @param error_bars Draw 95% confidence intervals of the observed values
#'   (observed value plus or minus 1.96 standard errors) when standard
#'   errors are available.
#' @param log_y Use a logarithmic y axis.
#' @param digits Decimals of the segment labels (default 1 for APCs and 3
#'   for slopes of linear models).
#' @param ncol Number of facet columns.
#' @param ... Ignored.
#' @return A ggplot object.
#' @importFrom ggplot2 autoplot
#' @export
#' @examples
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, year, rate, kmax = 3, ci_method = "parametric")
#' ggplot2::autoplot(fit)
#' ggplot2::autoplot(fit, color_by = "trend", log_y = TRUE)
autoplot.joinpoint <- function(object, model = "final",
                               color_by = c("cohort", "trend"),
                               layout = c("facet", "overlay"),
                               cohorts = NULL, labels = NULL, aapc = FALSE,
                               band = TRUE, joinpoint_ci = FALSE,
                               error_bars = FALSE, log_y = FALSE,
                               digits = NULL, ncol = NULL, ...) {
  color_by <- match.arg(color_by)
  layout <- match.arg(layout)
  s <- object$settings
  nm <- jp_stat_names(s$interval, s$log, s$family)
  if (is.null(digits)) digits <- if (s$log) 1L else 3L
  pd <- jp_plot_data(object, model, cohorts, digits, band)
  overlay <- layout == "overlay"
  single <- nlevels(pd$points$cohort) == 1L
  colored <- overlay && !single
  if (overlay && color_by == "trend" && !single) {
    stop("`color_by = \"trend\"` needs `layout = \"facet\"`.", call. = FALSE)
  }
  labels <- match.arg(labels %||% if (colored) "none" else "apc",
                      c("apc", "ci", "none"))
  if (labels == "ci") {
    sg <- pd$segments
    pd$segments$label <- paste0(sg$label, "\n(", jp_fmt(sg$lower, digits),
                                ", ", jp_fmt(sg$upper, digits), ")")
  }
  title <- NULL
  aapc_cap <- NULL
  if (aapc && !is.null(pd$aapc)) {
    a <- pd$aapc
    # A range shared by every cohort is named once, in the caption.
    shared <- length(unique(paste(a$start, a$end))) == 1L && !single
    txt <- paste0(nm$aapc_abbr, if (!shared) paste0(" ", a$start, "-", a$end),
                  ": ", jp_fmt(a$estimate, digits),
                  ifelse(a$significant, "*", ""), " (",
                  jp_fmt(a$lower, digits), ", ", jp_fmt(a$upper, digits), ")")
    if (shared) {
      aapc_cap <- paste0(nm$aapc_abbr, ": ", a$start[1L], "-", a$end[1L],
                         "; * significant at alpha = ", s$control$alpha_aapc,
                         ".")
    }
    if (single) {
      title <- txt
    } else {
      lv <- levels(pd$points$cohort)
      new <- stats::setNames(lv, lv)
      new[as.character(a$cohort)] <- paste0(a$cohort, "\n", txt)
      for (nmd in c("points", "curve", "joinpoints", "segments", "jump")) {
        if (!is.null(pd[[nmd]])) {
          pd[[nmd]]$cohort <- factor(new[as.character(pd[[nmd]]$cohort)],
                                     unname(new))
        }
      }
    }
  }
  # Every layer maps color to `key`: the cohort when cohorts are overlaid,
  # otherwise a role (observed, fitted, joinpoint) or the trend direction.
  key <- function(d, role) {
    d$key <- if (colored) d$cohort else role
    d
  }
  curve <- if (color_by == "trend") {
    cv <- pd$curve
    cv$key <- pd$segments$trend[match(paste(cv$cohort, cv$segment),
                                      paste(pd$segments$cohort,
                                            pd$segments$segment))]
    cv
  } else {
    key(pd$curve, "Fitted")
  }
  pts <- key(pd$points, "Observed")

  g <- ggplot2::ggplot(mapping = ggplot2::aes(x = .data$x, y = .data$y,
                                              color = .data$key))
  if (joinpoint_ci && !is.null(pd$joinpoints)) {
    g <- g + ggplot2::geom_rect(
      data = pd$joinpoints[!is.na(pd$joinpoints$lower), ], inherit.aes = FALSE,
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
      ymin = -Inf, ymax = Inf, fill = jp_ink[["grid"]], alpha = 0.6
    )
  }
  if (!is.null(pd$jump)) {
    g <- g + ggplot2::geom_vline(data = pd$jump,
                                 ggplot2::aes(xintercept = .data$x),
                                 linetype = "dotted", color = jp_ink[["muted"]]) +
      ggplot2::geom_text(data = pd$jump,
                         ggplot2::aes(x = .data$x, label = .data$label),
                         y = Inf, vjust = 1.4, hjust = -0.08, size = 3,
                         color = jp_ink[["secondary"]], inherit.aes = FALSE)
  }
  if (band && any(!is.na(curve$lower))) {
    g <- g + ggplot2::geom_ribbon(
      data = curve, ggplot2::aes(ymin = .data$lower, ymax = .data$upper,
                                 y = NULL, fill = .data$key,
                                 group = interaction(.data$cohort, .data$piece)),
      color = NA, alpha = 0.15, show.legend = FALSE
    )
  }
  if (error_bars && any(!is.na(pts$se))) {
    g <- g + ggplot2::geom_errorbar(
      data = pts[!is.na(pts$se), ],
      ggplot2::aes(ymin = .data$y - 1.96 * .data$se,
                   ymax = .data$y + 1.96 * .data$se),
      width = 0, linewidth = 0.4, alpha = 0.7, show.legend = FALSE
    )
  }
  g <- g +
    ggplot2::geom_point(data = pts[pts$included, ], shape = 16, size = 1.6,
                        alpha = 0.55, show.legend = FALSE) +
    ggplot2::geom_point(data = pts[!pts$included, ], shape = 1, size = 1.9,
                        stroke = 0.6, show.legend = FALSE) +
    ggplot2::geom_line(data = curve, ggplot2::aes(y = .data$fit,
                                                  group = interaction(.data$cohort, .data$piece)),
                       linewidth = 0.9, lineend = "round")
  pad <- if (labels == "none") 0.05 else if (labels == "ci") 0.2 else 0.12
  if (labels != "none") {
    # Labels are placed panel by panel, jointly for overlaid cohorts, with
    # text sizes judged for a figure about 7.5 by 5 inches; the panels span
    # the data ranges plus the axis expansions.
    tf <- if (log_y) log10 else identity
    np <- if (overlay || single) 1L else nlevels(pd$points$cohort)
    nc <- if (np == 1L) 1L else ncol %||% ceiling(sqrt(np))
    panel <- function(d) if (overlay) rep(1L, nrow(d)) else d$cohort
    sg <- do.call(rbind, lapply(split(pd$segments, panel(pd$segments),
                                      drop = TRUE), function(d) {
      on <- panel(pts) == panel(d)[1L]
      cv <- curve[panel(curve) == panel(d)[1L] & !is.na(curve$fit), ]
      jp_place_labels(d, range(pts$x[on]),
                      range(tf(c(pts$y[on], cv$fit)), na.rm = TRUE), tf,
                      cv$x, tf(cv$fit), cw = 1.1 * 0.065 * nc / 6.9,
                      lh = (1 + 2 * pad) * 0.14 * ceiling(np / nc) / 3.8)
    }))
    g <- g + ggplot2::geom_label(
      data = sg, ggplot2::aes(y = .data$label_y, label = .data$label,
                              hjust = .data$hjust),
      color = jp_ink[["secondary"]], size = 3,
      fill = grDevices::adjustcolor("white", 0.8), border.colour = NA,
      label.padding = ggplot2::unit(1.5, "pt")
    )
  }
  # Joinpoint markers are drawn last so that no label hides them.
  if (!is.null(pd$joinpoints)) {
    g <- g + ggplot2::geom_point(data = key(pd$joinpoints, "Joinpoint"),
                                 shape = 21, size = 2.2, stroke = 0.8,
                                 fill = "white", show.legend = FALSE)
  }
  roles <- c(Observed = jp_ink[["secondary"]], Fitted = jp_pal[1L],
             Joinpoint = jp_ink[["primary"]])
  g <- g + if (colored) {
    list(scale_color_joinpoint(name = NULL, aesthetics = c("colour", "fill")),
         if (aapc) ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2)))
  } else if (color_by == "trend") {
    ggplot2::scale_color_manual(
      values = c(jp_trend_cols, roles), name = NULL,
      breaks = intersect(names(jp_trend_cols), as.character(curve$key)),
      aesthetics = c("colour", "fill")
    )
  } else {
    ggplot2::scale_color_manual(values = roles, guide = "none",
                                aesthetics = c("colour", "fill"))
  }
  if (!overlay && !single) {
    g <- g + ggplot2::facet_wrap(~cohort, scales = "free_y", ncol = ncol)
  }
  g <- g + if (log_y) {
    ggplot2::scale_y_log10(expand = ggplot2::expansion(mult = pad))
  } else {
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = pad))
  }
  cap <- c(
    if (labels != "none") paste0("Labels: ", nm$apc, " of each segment; * ",
                                 "significant at alpha = ",
                                 s$control$alpha_apc, "."),
    if (band && any(!is.na(curve$lower))) "Band: 95% empirical quantile interval.",
    if (error_bars) "Error bars: 95% confidence intervals of observed values.",
    aapc_cap
  )
  g + ggplot2::labs(x = jp_title(s$x), y = jp_ylab(s), title = title,
                    caption = if (length(cap)) paste(cap, collapse = "\n")) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_joinpoint()
}

# Positions of the segment labels of one panel, in panel fractions: each
# label goes next to the middle of its segment on the side facing the middle
# of the panel, else on the other side, else stacked further out, at the
# first position that overlaps neither a label already placed nor a fitted
# curve (`cx`, `cy`); a two-line label that fits nowhere is cut to its
# first line, and a label that still fits nowhere takes the position of
# least overlap. `cw` and `lh` are the width of a character and the height
# of a line as fractions of the panel.
jp_place_labels <- function(d, xr, yr, tf, cx, cy, cw, lh) {
  d <- d[order(d$x), ]
  X <- (d$x - xr[1L]) / diff(xr)
  Y <- (tf(d$y) - yr[1L]) / diff(yr)
  CX <- (cx - xr[1L]) / diff(xr)
  CY <- (cy - yr[1L]) / diff(yr)
  d$hjust <- ifelse(X < 0.08, 0, ifelse(X > 0.92, 1, 0.5))
  placed <- matrix(numeric(0), 0L, 4L)
  mid <- numeric(nrow(d))
  for (i in seq_len(nrow(d))) {
    lines <- strsplit(d$label[i], "\n")[[1L]]
    spot <- function(nl) {
      w <- (max(nchar(lines[seq_len(nl)])) + 1) * cw
      h <- nl * lh
      sides <- if (Y[i] < 0.5) c(1, -1) else c(-1, 1)
      # Candidates: nearest first, on either side, aligned as given or
      # shifted half a label to the left or right.
      cand <- expand.grid(shift = c(0, -0.5, 0.5), side = sides, level = 0:5)
      cand$hj <- pmin(1, pmax(0, d$hjust[i] - cand$shift))
      cand$left <- X[i] - cand$hj * w
      cand$lo <- ifelse(cand$side > 0, Y[i] + 0.02 + cand$level * h,
                        Y[i] - 0.02 - (cand$level + 1) * h)
      inside <- cand$lo > -0.18 & cand$lo + h < 1.18 &
        cand$left > -0.02 & cand$left + w < 1.02
      if (any(inside)) cand <- cand[inside, ]
      # Overlap with placed labels (area) and with the curves (points).
      clash <- mapply(function(l, lo) {
        sum(pmax(0, pmin(l + w, placed[, 2L]) - pmax(l, placed[, 1L])) *
              pmax(0, pmin(lo + h, placed[, 4L]) - pmax(lo, placed[, 3L]))) +
          1e-3 * sum(CX > l & CX < l + w & CY > lo & CY < lo + h)
      }, cand$left, cand$lo)
      b <- which.min(clash)
      list(nl = nl, clash = clash[b], hj = cand$hj[b],
           box = c(cand$left[b], cand$left[b] + w, cand$lo[b], h))
    }
    sp <- spot(length(lines))
    if (sp$clash > 0 && length(lines) > 1L) {
      one <- spot(1L)
      if (one$clash < sp$clash) sp <- one
    }
    box <- sp$box
    d$hjust[i] <- sp$hj
    d$label[i] <- paste(lines[seq_len(sp$nl)], collapse = "\n")
    placed <- rbind(placed, c(box[1:3], box[3L] + box[4L]))
    mid[i] <- box[3L] + box[4L] / 2
  }
  y <- yr[1L] + mid * diff(yr)
  d$label_y <- if (identical(tf, log10)) 10^y else y
  d
}

#' @rdname autoplot.joinpoint
#' @param x A `joinpoint` object.
#' @export
plot.joinpoint <- function(x, ...) {
  print(autoplot.joinpoint(x, ...))
  invisible(x)
}

#' Plot all fitted models
#'
#' The best fitting model for every number of joinpoints of one cohort,
#' either as small multiples (the selection statistic in each panel title)
#' or overlaid in one panel, with the final selected model highlighted.
#'
#' @param object A `joinpoint` object.
#' @param cohort Index of the cohort to draw.
#' @param overlay Draw all models in one panel.
#' @param ncol Number of columns of the small multiples.
#' @return A ggplot object.
#' @examples
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
#' plot_models(fit)
#' @export
plot_models <- function(object, cohort = 1L, overlay = FALSE, ncol = NULL) {
  f <- object$fits[[cohort]]
  if (is.null(f)) stop("Cohort was not fitted.", call. = FALSE)
  crit <- f$criteria
  stat <- switch(f$method, wbic = "wbic", mbic = "mbic", bic3 = "bic3",
                 "bic")
  ks <- as.numeric(names(f$models))
  titles <- vapply(ks, function(k) {
    sprintf("%d joinpoint%s%s\n%s = %s", k, if (k == 1) "" else "s",
            if (k == f$k_selected) " (selected)" else "", toupper(stat),
            formatC(crit[[stat]][crit$joinpoints == k], digits = 3,
                    format = "f"))
  }, character(1))
  if (overlay) titles <- sub("\n.*", "", titles)
  crv <- do.call(rbind, lapply(seq_along(ks), function(i) {
    cv <- jp_curve(f, f$models[[i]], band = FALSE)
    data.frame(cv[c("x", "piece", "fit")], model = titles[i],
               final = ks[i] == f$k_selected)
  }))
  crv$model <- factor(crv$model, titles)
  pts <- f$data
  g <- ggplot2::ggplot(crv, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_point(data = pts, ggplot2::aes(y = .data$r,
                                                 shape = .data$included),
                        color = jp_ink[["secondary"]], size = 1.3,
                        alpha = 0.6, show.legend = FALSE) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1))
  g <- if (overlay) {
    g + ggplot2::geom_line(ggplot2::aes(
      y = .data$fit, color = .data$model, linetype = .data$final,
      group = interaction(.data$model, .data$piece)
    ), linewidth = 0.8) +
      ggplot2::scale_linetype_manual(values = c(`TRUE` = "solid",
                                                `FALSE` = "22"),
                                     guide = "none") +
      scale_color_joinpoint(name = NULL)
  } else {
    g + ggplot2::geom_line(ggplot2::aes(y = .data$fit, color = .data$final,
                                        group = .data$piece),
                           linewidth = 0.8, show.legend = FALSE) +
      ggplot2::scale_color_manual(values = c(`TRUE` = jp_pal[1L],
                                             `FALSE` = jp_ink[["muted"]])) +
      ggplot2::facet_wrap(~model, ncol = ncol)
  }
  # Poor fits of small models can run far outside the data; the axis follows
  # the observations.
  g + ggplot2::coord_cartesian(ylim = range(pts$r)) +
    ggplot2::labs(x = jp_title(object$settings$x),
                  y = jp_ylab(object$settings),
                  title = jp_cohort_label(object, cohort)) +
    theme_joinpoint()
}

#' Plot segment APCs and AAPCs
#'
#' `plot_apc()` draws, for each cohort, the APC of every segment as a step
#' over the years it covers, with its confidence interval as a shaded band
#' and a zero reference line, so changes in trend read directly off the
#' time axis. `plot_aapc()` draws a forest plot of AAPCs with confidence
#' intervals by cohort and range.
#'
#' @param object A `joinpoint` object.
#' @param model `"final"` or a number of joinpoints.
#' @param ncol Number of facet columns.
#' @return A ggplot object.
#' @examples
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
#' plot_apc(fit)
#' plot_aapc(fit)
#' @export
plot_apc <- function(object, model = "final", ncol = NULL) {
  pd <- jp_plot_data(object, model)
  sg <- pd$segments
  nm <- jp_stat_names(object$settings$interval, object$settings$log,
                       object$settings$family)
  ylab <- if (object$settings$log) paste0(nm$apc, " (%)") else nm$apc
  g <- ggplot2::ggplot(sg) +
    ggplot2::geom_hline(yintercept = 0, color = jp_ink[["axis"]],
                        linewidth = 0.4) +
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$start, xmax = .data$end,
                                    ymin = .data$lower, ymax = .data$upper,
                                    fill = .data$trend), alpha = 0.18) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$start, xend = .data$end,
                                       y = .data$estimate,
                                       yend = .data$estimate,
                                       color = .data$trend),
                          linewidth = 1, lineend = "round") +
    ggplot2::geom_text(ggplot2::aes(x = .data$x, y = .data$upper,
                                    label = .data$label),
                       vjust = -0.6, size = 3, color = jp_ink[["secondary"]]) +
    ggplot2::scale_color_manual(values = jp_trend_cols, name = NULL,
                                drop = TRUE) +
    ggplot2::scale_fill_manual(values = jp_trend_cols, name = NULL,
                               drop = TRUE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = jp_title(object$settings$x), y = ylab) +
    theme_joinpoint()
  if (nlevels(sg$cohort) > 1L) {
    g <- g + ggplot2::facet_wrap(~cohort, scales = "free_y", ncol = ncol)
  }
  g
}

#' @rdname plot_apc
#' @export
plot_aapc <- function(object, model = "final") {
  a <- aapc(object, model)
  if (is.null(a)) {
    stop("No AAPC available.", call. = FALSE)
  }
  a$cohort <- if (is.null(object$by)) {
    "All"
  } else {
    do.call(paste, c(a[object$by], sep = ", "))
  }
  a$cohort <- factor(a$cohort, rev(unique(a$cohort)))
  a$range <- factor(a$range, unique(a$range))
  a$trend <- jp_trend(a$estimate, a$p_value, object$settings$control$alpha_aapc)
  # Each label sits above its interval, aligned inward near the panel edges.
  lo <- stats::ave(a$lower, a$range, FUN = min)
  hi <- stats::ave(a$upper, a$range, FUN = max)
  rel <- (a$estimate - lo) / pmax(hi - lo, 1e-12)
  a$hjust <- ifelse(rel < 0.3, 0, ifelse(rel > 0.7, 1, 0.5))
  nm <- jp_stat_names(object$settings$interval, object$settings$log,
                       object$settings$family)
  ggplot2::ggplot(a, ggplot2::aes(x = .data$estimate, y = .data$cohort)) +
    ggplot2::geom_vline(xintercept = 0, color = jp_ink[["axis"]],
                        linewidth = 0.4) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower,
                                        xmax = .data$upper,
                                        color = .data$trend),
                           width = 0, linewidth = 0.7, orientation = "y") +
    ggplot2::geom_point(ggplot2::aes(color = .data$trend), size = 2.4) +
    ggplot2::geom_text(ggplot2::aes(hjust = .data$hjust,
                                    label = paste0(jp_fmt(.data$estimate),
                                                   " (", jp_fmt(.data$lower),
                                                   ", ", jp_fmt(.data$upper),
                                                   ")")),
                       nudge_y = 0.3, size = 2.8,
                       color = jp_ink[["secondary"]]) +
    ggplot2::scale_color_manual(values = jp_trend_cols, name = NULL,
                                drop = TRUE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.1)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::facet_wrap(~range, scales = "free_x") +
    ggplot2::labs(x = paste0(nm$aapc, if (object$settings$log) " (%)"),
                  y = NULL) +
    theme_joinpoint() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(
                     color = jp_ink[["grid"]], linewidth = 0.3))
}

#' Plot model selection
#'
#' Selection statistic against the number of joinpoints for each cohort.
#' For the permutation test the sequence of tests is drawn instead: the
#' p value of each test against its significance level.
#'
#' @param object A `joinpoint` object.
#' @return A ggplot object.
#' @examples
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
#' plot_selection(fit)
#' @export
plot_selection <- function(object) {
  ms <- model_selection(object)
  method <- object$settings$method
  coh <- function(d) {
    if (is.null(object$by)) "All" else do.call(paste, c(d[object$by], sep = ", "))
  }
  if (method %in% c("permutation", "lrt")) {
    t <- ms$tests
    if (is.null(t)) stop("No tests were performed.", call. = FALSE)
    what <- if (method == "lrt") "Monte Carlo likelihood ratio" else "permutation"
    t$cohort <- coh(t)
    t$test <- factor(paste0(t$null, " vs ", t$alternative))
    return(
      ggplot2::ggplot(t, ggplot2::aes(x = .data$test, y = .data$p_value)) +
        ggplot2::geom_col(ggplot2::aes(fill = .data$reject), width = 0.6) +
        ggplot2::geom_point(ggplot2::aes(y = .data$alpha), shape = 95,
                            size = 8, color = jp_ink[["primary"]]) +
        ggplot2::scale_fill_manual(values = c(`TRUE` = jp_pal[1L],
                                              `FALSE` = jp_ink[["muted"]]),
                                   labels = c(`TRUE` = "Null rejected",
                                              `FALSE` = "Not rejected"),
                                   name = NULL) +
        ggplot2::facet_wrap(~cohort) +
        ggplot2::labs(x = "Test (joinpoints under H0 vs H1)", y = "p value",
                      caption = paste0("Bars: ", what, " p value; ticks: significance level of the test.")) +
        theme_joinpoint()
    )
  }
  cr <- ms$criteria
  cr$cohort <- coh(cr)
  stats <- intersect(c("bic", "bic3", "wbic", "mbic"), names(cr))
  long <- do.call(rbind, lapply(stats, function(s) {
    data.frame(cohort = cr$cohort, joinpoints = cr$joinpoints,
               criterion = toupper(s), value = cr[[s]],
               selected = cr$selected)
  }))
  long$criterion <- factor(long$criterion, toupper(stats))
  sel <- long[long$selected & long$criterion == toupper(switch(
    method, dds = "bic", permutation = "bic", method)), ]
  ggplot2::ggplot(long, ggplot2::aes(x = .data$joinpoints, y = .data$value,
                                     color = .data$criterion)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_point(data = sel, shape = 21, size = 4, stroke = 1,
                        color = jp_ink[["primary"]], fill = NA) +
    scale_color_joinpoint(name = NULL) +
    ggplot2::scale_x_continuous(breaks = function(l) {
      seq(ceiling(l[1L]), floor(l[2L]))
    }) +
    ggplot2::facet_wrap(~cohort, scales = "free_y") +
    ggplot2::labs(x = "Number of joinpoints", y = "Criterion",
                  caption = "Circled: final selected model.") +
    theme_joinpoint()
}

#' Residual diagnostics
#'
#' Standardized residuals of the selected model (weighted by the error
#' structure, so constant variance is expected) against the independent
#' variable, as a normal quantile plot, or as an autocorrelation function
#' to judge whether an autocorrelated errors model is warranted.
#'
#' @param object A `joinpoint` object.
#' @param type `"time"`, `"qq"` or `"acf"`.
#' @param model `"final"` or a number of joinpoints.
#' @return A ggplot object.
#' @examples
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
#' plot_residuals(fit)
#' plot_residuals(fit, type = "acf")
#' @export
plot_residuals <- function(object, type = c("time", "qq", "acf"),
                           model = "final") {
  type <- match.arg(type)
  d <- do.call(rbind, lapply(seq_along(object$fits), function(g) {
    f <- object$fits[[g]]
    if (is.null(f)) {
      return(NULL)
    }
    m <- jp_pick(f, model)
    if (sum(m$wresid^2) <= jp_tiny(jp_whiten_mat(f$S), f$y)) {
      stop("Cohort ", jp_cohort_label(object, g), " is fitted exactly; its ",
           "residuals are numerically zero.", call. = FALSE)
    }
    r <- m$wresid / sqrt(sum(m$wresid^2) /
                           max(1, m$n - 2 * (m$k + 1) - length(m$extra_names)))
    data.frame(cohort = jp_cohort_label(object, g), x = f$x, resid = r)
  }))
  d$cohort <- factor(d$cohort, unique(d$cohort))
  facet <- if (nlevels(d$cohort) > 1L) ggplot2::facet_wrap(~cohort)
  ref <- ggplot2::geom_hline(yintercept = 0, color = jp_ink[["axis"]],
                             linewidth = 0.4)
  g <- switch(type,
    time = ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$resid)) +
      ref +
      ggplot2::geom_segment(ggplot2::aes(xend = .data$x, yend = 0),
                            color = jp_ink[["muted"]], linewidth = 0.3) +
      ggplot2::geom_point(color = jp_pal[1L], size = 1.4) +
      ggplot2::labs(x = jp_title(object$settings$x), y = "Standardized residual"),
    qq = ggplot2::ggplot(d, ggplot2::aes(sample = .data$resid)) +
      ggplot2::geom_qq_line(color = jp_ink[["muted"]], linewidth = 0.5) +
      ggplot2::geom_qq(color = jp_pal[1L], size = 1.4) +
      ggplot2::labs(x = "Normal quantile", y = "Standardized residual"),
    acf = {
      a <- do.call(rbind, lapply(split(d, d$cohort), function(dd) {
        ac <- stats::acf(dd$resid, plot = FALSE,
                         lag.max = min(10L, nrow(dd) - 2L))
        data.frame(cohort = dd$cohort[1L], lag = drop(ac$lag)[-1L],
                   acf = drop(ac$acf)[-1L], bound = 1.96 / sqrt(nrow(dd)))
      }))
      ggplot2::ggplot(a, ggplot2::aes(x = .data$lag, y = .data$acf)) + ref +
        ggplot2::scale_x_continuous(breaks = function(l) {
          seq(ceiling(l[1L]), floor(l[2L]))
        }) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = .data$bound),
                            linetype = "dashed", color = jp_ink[["muted"]]) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = -.data$bound),
                            linetype = "dashed", color = jp_ink[["muted"]]) +
        ggplot2::geom_segment(ggplot2::aes(xend = .data$lag, yend = 0),
                              color = jp_pal[1L], linewidth = 0.9) +
        ggplot2::labs(x = "Lag", y = "Autocorrelation",
                      caption = "Dashed: approximate 95% limits under no autocorrelation.")
    }
  )
  g + facet + theme_joinpoint()
}

# Axis titles: column names with spaces and a capital first letter.
jp_title <- function(v) {
  v <- gsub("_", " ", v)
  paste0(toupper(substr(v, 1L, 1L)), substring(v, 2L))
}

# Default y axis title from the response specification.
jp_ylab <- function(s) {
  if (!is.null(s$y)) {
    return(jp_title(s$y))
  }
  switch(s$type,
    crude_rate = paste("Rate per", format(s$rate_per, big.mark = ",",
                                          scientific = FALSE)),
    proportion = "Proportion", percent = "Percent", count = "Count", "Value"
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a
