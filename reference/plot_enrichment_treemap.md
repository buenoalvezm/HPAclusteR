# Plot treemap for each cluster in reduced GO terms

Plot treemap for each cluster in reduced GO terms

## Usage

``` r
plot_enrichment_treemap(
  reduced_terms,
  ontology_pal = c(MF = "#FF6F00", CC = "#C71B00", BP = "#018EA0")
)
```

## Arguments

- reduced_terms:

  Data frame, e.g. simp_res\$reducedTerms

- ontology_pal:

  Named vector of ontology colors (default: c(MF="#FF6F00",
  CC="#C71B00", BP="#018EA0"))

## Value

Named list of ggplot treemaps, one per cluster
