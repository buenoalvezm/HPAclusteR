# A simple theme for visualizations

`theme_hc()` provides a minimalist ggplot2 theme with no background or
grid lines. Usefull for clean, publication-ready point, bar, line and
boxplots among others.

## Usage

``` r
theme_hc()
```

## Value

A ggplot2 theme object that can be added to plots with `+`.

## Examples

``` r
library(ggplot2)
ggplot2::ggplot(mpg, aes(displ, hwy)) +
  ggplot2::geom_point() +
  theme_hc()
```
