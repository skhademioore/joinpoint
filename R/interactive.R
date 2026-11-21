# Interactive graphics: the ggplot2 figures of the package drawn with the
# interactive layers of ggiraph (tooltips, linked highlighting, zoom) and
# shown as htmlwidgets.

#' Interactive joinpoint graphics
#'
#' Draws a joinpoint figure as an interactive graphic: an htmlwidget made
#' with the `ggiraph` package that looks like the static figure and adds
#' tooltips, linked highlighting and zoom. Hovering a fitted segment shows
#' its period and APC with confidence interval and p value, and highlights
#' the segment with its label and band while dimming the rest; hovering an
#' observation shows its value, confidence interval and fitted trend; hovering
#' a joinpoint shows its location and confidence interval. The toolbar zooms
#' and saves the graphic as a PNG image.
#'
#' The widget works in the RStudio viewer, R Markdown and Quarto HTML
#' documents, pkgdown sites and Shiny applications (with
#' `ggiraph::girafeOutput()` and `ggiraph::renderGirafe()`). [jp_save()]
#' writes it to a standalone HTML file. Documents rendered to Word or PDF
#' need the static figures ([autoplot()][autoplot.joinpoint] and the other
#' plot functions).
#'
#' @param x A `joinpoint` object, or a ggplot made by one of the plot
#'   functions of the package with `interactive = TRUE` (and possibly
#'   modified, for example with further labels or theme settings).
#' @param plot The figure drawn from a `joinpoint` object: `"trend"`
#'   ([autoplot()][autoplot.joinpoint]), `"apc"` ([plot_apc()]), `"aapc"`
#'   ([plot_aapc()]) or `"selection"` ([plot_selection()]).
#' @param ... Further arguments to the plot function, for example
#'   `color_by = "trend"` or `log_y = TRUE`.
#' @param width,height Size of the graphic in inches; the widget scales to
#'   the width of its container, keeping this aspect ratio.
#' @param zoom Maximum zoom factor (`1` for no zoom).
#' @param embed_fonts Embed the font files in the widget, so that it looks
#'   exactly the same on every computer. By default the text uses fonts
#'   found on the system (Arial, Helvetica or DejaVu Sans, for example),
#'   which keeps the widget small: about 0.2 MB instead of several MB.
#' @return An htmlwidget of class `girafe`.
#' @examplesIf requireNamespace("ggiraph", quietly = TRUE)
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, year, rate, kmax = 3, ci_method = "parametric")
#' jp_interactive(fit, log_y = TRUE)
#' jp_interactive(fit, plot = "apc")
#'
#' # A modified figure.
#' g <- ggplot2::autoplot(fit, interactive = TRUE) +
#'   ggplot2::labs(title = "Tuberculosis mortality")
#' jp_interactive(g)
#' @export
jp_interactive <- function(x, plot = c("trend", "apc", "aapc", "selection"),
                           ..., width = 7.5, height = 5, zoom = 5,
                           embed_fonts = FALSE) {
  jp_need("ggiraph", "interactive graphics")
  if (inherits(x, "joinpoint")) {
    plot <- match.arg(plot)
    x <- switch(plot,
      trend = autoplot.joinpoint(x, ..., interactive = TRUE),
      apc = plot_apc(x, ..., interactive = TRUE),
      aapc = plot_aapc(x, ..., interactive = TRUE),
      selection = plot_selection(x, ..., interactive = TRUE)
    )
  } else if (!inherits(x, "ggplot")) {
    stop("`x` must be a joinpoint fit or a ggplot.", call. = FALSE)
  }
  if (!is.numeric(zoom) || length(zoom) != 1L || is.na(zoom) || zoom < 1) {
    stop("`zoom` must be a number of at least 1.", call. = FALSE)
  }
  jp_need("gdtools", "interactive graphics")
  # gdtools may set this option while it builds a font set.
  old <- options(GFONTS_DOWNLOAD_SLEEPTIME = getOption("GFONTS_DOWNLOAD_SLEEPTIME"))
  on.exit(options(old), add = TRUE)
  fs <- if (isTRUE(embed_fonts)) gdtools::font_set_auto() else jp_system_fonts()
  w <- ggiraph::girafe(
    ggobj = x, width_svg = width, height_svg = height, font_set = fs,
    canvas_id = jp_canvas_id(),
    options = list(
      ggiraph::opts_tooltip(css = jp_tooltip_css, opacity = 1,
                            use_fill = FALSE, offx = 12, offy = 12),
      ggiraph::opts_hover(css = "stroke-width:2.5px;opacity:1;"),
      ggiraph::opts_hover_inv(css = "opacity:0.3;"),
      ggiraph::opts_selection(type = "none"),
      ggiraph::opts_zoom(min = 1, max = zoom),
      ggiraph::opts_toolbar(saveaspng = TRUE, position = "topright"),
      ggiraph::opts_sizing(rescale = TRUE)
    )
  )
  if (!isTRUE(embed_fonts)) {
    fonts <- vapply(fs$dependencies, function(d) d$name, character(1))
    w$dependencies <- Filter(function(d) !d$name %in% fonts, w$dependencies)
  }
  w
}

# Fonts installed on the computer, named in the graphic without copies of
# their files (nothing is written to the font cache of gdtools); the viewer
# draws the text with its own fonts of these names or their generic
# family. Without any of the usual system fonts, the font set of gdtools
# is used.
jp_system_fonts <- function() {
  pick <- function(v) {
    ok <- vapply(v, function(f) isTRUE(gdtools::font_family_exists(f)),
                 logical(1))
    if (any(ok)) v[ok][1L]
  }
  fam <- list(sans = pick(c("Arial", "Helvetica", "DejaVu Sans")),
              serif = pick(c("Times New Roman", "Times", "DejaVu Serif")),
              mono = pick(c("Courier New", "Courier", "DejaVu Sans Mono")))
  if (is.null(fam$sans) || is.null(fam$serif) || is.null(fam$mono)) {
    return(gdtools::font_set_auto())
  }
  fam$symbol <- pick(c("Symbol", "Apple Symbols")) %||% fam$sans
  structure(c(fam, list(dependencies = list(), dsvg_fonts = fam,
                        sources = lapply(fam, function(f) "system"))),
            class = "font_set")
}

# Identifier of the SVG canvas of a widget, unique within the session and
# made without the random number generator (ggiraph would otherwise draw
# it and change the user's random number stream).
jp_state <- new.env(parent = emptyenv())
jp_canvas_id <- function() {
  jp_state$n <- (jp_state$n %||% 0L) + 1L
  paste0("jp", Sys.getpid(), "_", sprintf("%.0f", as.numeric(Sys.time()) * 1e3),
         "_", jp_state$n)
}

jp_tooltip_css <- paste(
  "background:#ffffff;color:#0b0b0b;border:1px solid #c3c2b7;",
  "border-radius:4px;padding:6px 8px;font-family:sans-serif;",
  "font-size:12px;line-height:1.35;box-shadow:0 1px 4px rgba(0,0,0,0.15);"
)

# Stops with an informative message when a suggested package is missing.
jp_need <- function(pkg, what) {
  ok <- tryCatch(loadNamespace(pkg), error = function(e) e)
  if (inherits(ok, "error")) {
    stop("The ", pkg, " package is needed for ", what, " but could not be ",
         "loaded: ", conditionMessage(ok), call. = FALSE)
  }
}

# The ggplot2 geom `name`, or its interactive ggiraph version.
jp_geom <- function(name, interactive) {
  if (interactive) {
    getExportedValue("ggiraph", paste0(name, "_interactive"))
  } else {
    getExportedValue("ggplot2", name)
  }
}

# Mapping `m` plus the tooltip and data_id aesthetics of interactive layers.
jp_ia_aes <- function(m, interactive) {
  if (!interactive) {
    return(m)
  }
  utils::modifyList(m, ggplot2::aes(tooltip = .data$tooltip,
                                    data_id = .data$data_id))
}

jp_html <- function(v) {
  v <- gsub("&", "&amp;", as.character(v), fixed = TRUE)
  v <- gsub("<", "&lt;", v, fixed = TRUE)
  gsub(">", "&gt;", v, fixed = TRUE)
}

jp_num <- function(v) {
  out <- formatC(signif(v, 6), digits = 6, format = "fg", drop0trailing = TRUE)
  out[is.na(v)] <- "NA"
  trimws(out)
}

jp_tip_p <- function(p) {
  ifelse(is.na(p), "", paste0("<br>p ", jp_pformat(p)))
}

jp_tip_ci <- function(lower, upper, digits, alpha = 0.05) {
  ifelse(is.na(lower) | is.na(upper), "",
         paste0(" (", jp_level(alpha), " CI ", jp_fmt(lower, digits), ", ",
                jp_fmt(upper, digits), ")"))
}

# Tooltips and data ids of the trend figure: segments (shared by the curve,
# its band and its label), observations and joinpoints.
jp_trend_tips <- function(pd, nm, digits, s) {
  head <- function(cohort) paste0("<b>", jp_html(cohort), "</b><br>")
  xl <- jp_html(jp_title(s$x))
  sg <- pd$segments
  sg$data_id <- paste0(sg$cohort, "|segment|", sg$segment)
  sg$tooltip <- paste0(
    head(sg$cohort), "Segment ", sg$segment, ": ", jp_num(sg$start), "-",
    jp_num(sg$end), "<br>", jp_html(nm$apc_abbr), " ", jp_fmt(sg$estimate, digits),
    ifelse(sg$trend %in% c("Increasing", "Decreasing"), "*", ""),
    jp_tip_ci(sg$lower, sg$upper, digits, s$control$alpha_apc),
    jp_tip_p(sg$p_value), "<br>",
    as.character(sg$trend)
  )
  pd$segments <- sg
  key <- match(paste(pd$curve$cohort, pd$curve$segment),
               paste(sg$cohort, sg$segment))
  pd$curve$data_id <- sg$data_id[key]
  pd$curve$tooltip <- sg$tooltip[key]
  pt <- pd$points
  pt$data_id <- paste0(pt$cohort, "|x|", pt$x)
  pt$tooltip <- paste0(
    head(pt$cohort), xl, " ", jp_num(pt$x), "<br>Observed ", jp_num(pt$y),
    ifelse(is.na(pt$lo), "",
           paste0(" (95% CI ", jp_num(pt$lo), ", ", jp_num(pt$hi), ")")),
    "<br>Fitted trend ", jp_num(pt$fit),
    ifelse(pt$included, "", "<br><i>Excluded from the fit</i>")
  )
  pd$points <- pt
  if (!is.null(pd$joinpoints)) {
    jp <- pd$joinpoints
    jp$data_id <- paste0(jp$cohort, "|joinpoint|", jp$x)
    jp$tooltip <- paste0(
      head(jp$cohort), "Joinpoint ", jp_num(jp$x),
      ifelse(is.na(jp$lower), "",
             paste0(" (", jp_level(s$control$alpha_joinpoint), " CI ",
                    jp_num(jp$lower), ", ", jp_num(jp$upper), ")")),
      "<br>Fitted ", jp_num(jp$y)
    )
    pd$joinpoints <- jp
  }
  pd
}
