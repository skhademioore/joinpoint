# Theme and color scales for joinpoint graphics

`theme_joinpoint()` is a clean publication theme with hairline
horizontal gridlines and recessive axes. `scale_color_joinpoint()` and
`scale_fill_joinpoint()` use an eight color qualitative palette whose
order keeps adjacent colors distinguishable under color vision
deficiency.

## Usage

``` r
theme_joinpoint(base_size = 11, base_family = "")

scale_color_joinpoint(..., aesthetics = "colour")

scale_fill_joinpoint(...)
```

## Arguments

- base_size:

  Base font size in points.

- base_family:

  Base font family.

- ...:

  Passed to
  [`ggplot2::discrete_scale()`](https://ggplot2.tidyverse.org/reference/discrete_scale.html).

- aesthetics:

  Aesthetics the color scale applies to.

## Value

A ggplot2 theme or scale.

## Examples

``` r
library(ggplot2)
ggplot(us_mortality, aes(year, rate, color = cause)) +
  geom_line() +
  scale_color_joinpoint() +
  theme_joinpoint()
```
