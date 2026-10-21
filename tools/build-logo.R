#!/usr/bin/env Rscript

# Rebuild the package hex logo, man/figures/logo.png and logo.svg.
# Run from the repository root: Rscript tools/build-logo.R [output-directory]
#
# Motif: a trend that falls steeply between two joinpoints and then levels
# off, drawn in white on the ink field with blue rings at the joinpoints and
# gray observations scattered behind it. The
# wordmark is drawn as outlines of the bundled Roboto Medium font
# (tools/Roboto-Medium.ttf, Apache License 2.0, tools/Roboto-LICENSE.txt), so
# the build depends on nothing outside the repository and the SVG embeds no
# font.

for (package in c("grid", "systemfonts", "ragg", "svglite")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package '", package, "' is required to rebuild the logo.",
         call. = FALSE)
  }
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.", call. = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) {
  args[[1L]]
} else {
  file.path(dirname(script_dir), "man", "figures")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

font_file <- file.path(script_dir, "Roboto-Medium.ttf")
systemfonts::register_font("joinpoint-logo", plain = font_file)

width_in <- 1.6
height_in <- 1.85
color <- list(blue = "#2A78D6", ink = "#0B2540", gray = "#7F8FA6",
              white = "#FFFFFF")

# Glyph outlines of the wordmark, shaped once and reused by both devices.
shaped <- systemfonts::shape_string("joinpoint", family = "joinpoint-logo",
                                    size = 1000)
letters <- systemfonts::glyph_outline(shaped$shape$index, path = font_file,
                                      size = 1000, tolerance = 0.08)
letters$x <- letters$x + shaped$shape$x_offset[letters$glyph]
letters$y <- letters$y + shaped$shape$y_offset[letters$glyph]
letters$id <- as.integer(interaction(letters$glyph, letters$contour,
                                     drop = TRUE))

draw_logo <- function() {
  grid::grid.newpage()
  # Hexagon with a vertex at the top.
  angle <- (90 + (0:5) * 60) * pi / 180
  grid::grid.polygon(
    x = grid::unit(width_in / 2 + 0.885 * cos(angle), "inches"),
    y = grid::unit(height_in / 2 + 0.885 * sin(angle), "inches"),
    gp = grid::gpar(fill = color$ink, col = color$blue, lwd = 2.5,
                    linejoin = "mitre")
  )
  x <- c(0.18, 0.40, 0.64, 0.82)
  y <- c(0.77, 0.70, 0.48, 0.48)
  # Observations around the trend, with fixed offsets so the build is
  # reproducible.
  px <- seq(0.19, 0.81, length.out = 13)
  py <- stats::approx(x, y, px)$y +
    c(0.04, -0.036, 0.03, 0.045, -0.034, 0.04, -0.042, 0.032, -0.038, 0.04,
      -0.03, 0.036, -0.04)
  grid::grid.circle(px, py, r = grid::unit(0.024, "inches"),
                    gp = grid::gpar(fill = color$gray, col = NA))
  grid::grid.lines(x, y, gp = grid::gpar(col = color$white, lwd = 4.2,
                                         lineend = "round",
                                         linejoin = "round"))
  grid::grid.circle(x[2:3], y[2:3], r = grid::unit(0.046, "inches"),
                    gp = grid::gpar(fill = color$ink, col = color$blue,
                                    lwd = 2.4))
  # Wordmark centered at 27% of the height, 62.5% of the width wide.
  scale <- 0.625 * width_in / diff(range(letters$x))
  grid::grid.path(
    x = grid::unit(width_in / 2 + (letters$x - mean(range(letters$x))) * scale,
                   "inches"),
    y = grid::unit(0.273 * height_in +
                     (letters$y - mean(range(letters$y))) * scale, "inches"),
    id = letters$id, rule = "winding",
    gp = grid::gpar(fill = color$white, col = NA)
  )
}

ragg::agg_png(file.path(output_dir, "logo.png"), width = width_in,
              height = height_in, units = "in", res = 300,
              background = "transparent")
draw_logo()
invisible(grDevices::dev.off())
svglite::svglite(file.path(output_dir, "logo.svg"), width = width_in,
                 height = height_in, bg = "transparent")
draw_logo()
invisible(grDevices::dev.off())
