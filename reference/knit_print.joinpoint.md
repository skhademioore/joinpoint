# Printing joinpoint fits in R Markdown and Quarto

A fit printed in a chunk of an R Markdown or Quarto document shows the
same summary as in the console. With
`options(joinpoint.knit_table = TRUE)` it appears instead as the
publication table of
[`jp_table()`](https://skhademioore.github.io/joinpoint/reference/jp_table.md),
formatted for the output (HTML, LaTeX, Word or Markdown); the chunk
option `tab.cap` sets its caption.

## Usage

``` r
# S3 method for class 'joinpoint'
knit_print(x, options = NULL, ...)
```

## Arguments

- x:

  A `joinpoint` object.

- options:

  Chunk options (passed by knitr).

- ...:

  Further arguments to the table method.

## Value

The output of
[`knitr::knit_print()`](https://rdrr.io/pkg/knitr/man/knit_print.html).
