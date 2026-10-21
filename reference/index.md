# Package index

## Fit joinpoint models

Fit joinpoint models to one or more trend series, and control the grid
search, the model selection and the confidence intervals.

- [`joinpoint()`](https://skhademioore.github.io/joinpoint/reference/joinpoint.md)
  : Fit joinpoint regression models
- [`jp_control()`](https://skhademioore.github.io/joinpoint/reference/jp_control.md)
  : Control settings for joinpoint fitting

## Trend estimates

Segment APCs, AAPCs and joinpoint locations with their confidence
intervals, model selection statistics, regression coefficients and a
table for a manuscript.

- [`apc()`](https://skhademioore.github.io/joinpoint/reference/apc.md)
  [`aapc()`](https://skhademioore.github.io/joinpoint/reference/apc.md)
  [`joinpoints()`](https://skhademioore.github.io/joinpoint/reference/apc.md)
  : Trend estimates from a joinpoint fit
- [`jp_table()`](https://skhademioore.github.io/joinpoint/reference/jp_table.md)
  : Publication table of trends
- [`model_selection()`](https://skhademioore.github.io/joinpoint/reference/model_selection.md)
  : Model selection results
- [`model_estimates()`](https://skhademioore.github.io/joinpoint/reference/model_estimates.md)
  : Model estimates
- [`jump_estimates()`](https://skhademioore.github.io/joinpoint/reference/jump_estimates.md)
  : Jump and comparability ratio estimates

## Model methods

Standard methods for fitted joinpoint models.

- [`as.data.frame(`*`<joinpoint>`*`)`](https://skhademioore.github.io/joinpoint/reference/as.data.frame.joinpoint.md)
  : Observed and modeled values
- [`predict(`*`<joinpoint>`*`)`](https://skhademioore.github.io/joinpoint/reference/predict.joinpoint.md)
  : Predictions from a joinpoint fit
- [`vcov(`*`<joinpoint>`*`)`](https://skhademioore.github.io/joinpoint/reference/vcov.joinpoint.md)
  : Covariance of the segment parameters

## Comparing and clustering trends

Tests of parallelism and coincidence, direct age adjustment and
clustering of contiguous age groups.

- [`jp_compare()`](https://skhademioore.github.io/joinpoint/reference/jp_compare.md)
  : Pairwise comparison of joinpoint trends
- [`age_adjust()`](https://skhademioore.github.io/joinpoint/reference/age_adjust.md)
  : Directly age adjusted rates
- [`jp_cluster()`](https://skhademioore.github.io/joinpoint/reference/jp_cluster.md)
  : Multi-group clustering of age specific trends

## Graphics

Publication quality figures built with ggplot2.

- [`autoplot(`*`<joinpoint>`*`)`](https://skhademioore.github.io/joinpoint/reference/autoplot.joinpoint.md)
  [`plot(`*`<joinpoint>`*`)`](https://skhademioore.github.io/joinpoint/reference/autoplot.joinpoint.md)
  : Plot a joinpoint fit
- [`plot_apc()`](https://skhademioore.github.io/joinpoint/reference/plot_apc.md)
  [`plot_aapc()`](https://skhademioore.github.io/joinpoint/reference/plot_apc.md)
  : Plot segment APCs and AAPCs
- [`plot_models()`](https://skhademioore.github.io/joinpoint/reference/plot_models.md)
  : Plot all fitted models
- [`plot_selection()`](https://skhademioore.github.io/joinpoint/reference/plot_selection.md)
  : Plot model selection
- [`plot_residuals()`](https://skhademioore.github.io/joinpoint/reference/plot_residuals.md)
  : Residual diagnostics
- [`autoplot(`*`<jp_comparison>`*`)`](https://skhademioore.github.io/joinpoint/reference/autoplot.jp_comparison.md)
  : Plot a pairwise comparison
- [`theme_joinpoint()`](https://skhademioore.github.io/joinpoint/reference/theme_joinpoint.md)
  [`scale_color_joinpoint()`](https://skhademioore.github.io/joinpoint/reference/theme_joinpoint.md)
  [`scale_fill_joinpoint()`](https://skhademioore.github.io/joinpoint/reference/theme_joinpoint.md)
  : Theme and color scales for joinpoint graphics

## Export

- [`jp_export()`](https://skhademioore.github.io/joinpoint/reference/jp_export.md)
  : Export joinpoint results

## Data sets

- [`us_mortality`](https://skhademioore.github.io/joinpoint/reference/us_mortality.md)
  : Age adjusted death rates for six causes, United States 1900 to 1998
- [`ontario_bmi`](https://skhademioore.github.io/joinpoint/reference/ontario_bmi.md)
  : Prepregnancy body mass index in Ontario, 2012 to 2021
- [`std_populations`](https://skhademioore.github.io/joinpoint/reference/std_populations.md)
  : Standard populations for age adjustment
