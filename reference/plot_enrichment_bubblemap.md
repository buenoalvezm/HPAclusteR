# Plot enrichment bubblemap for annotation results

Plot enrichment bubblemap for annotation results

## Usage

``` r
plot_enrichment_bubblemap(
  enrichment_results,
  odds_ratio_cap = 50,
  size_range = c(1, 4),
  fill_palette = c("white", "orangered"),
  facet_var = NULL
)
```

## Arguments

- enrichment_results:

  Data frame with annotation results (e.g. database_enrichment,
  kegg_enrichment, or go_enrichment)

- odds_ratio_cap:

  Maximum value for odds ratio to cap bubble size (default: 50)

- size_range:

  Numeric vector of length 2 for point size range (default: c(1, 4))

- fill_palette:

  Color palette for odds ratio (default: c("white", "orangered"))

- facet_var:

  Optional variable to facet by (default: NULL)

## Value

ggplot object (bubblemap)
