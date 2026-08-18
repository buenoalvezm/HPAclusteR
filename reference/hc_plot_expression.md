# Plot per-cluster heatmaps for scaled and z-score expression

`hc_plot_expression()` creates heatmaps of gene expression for each
cluster identified in the AnnDatR object. It generates two types of
plots per cluster: one showing z-score normalized expression and another
showing scaled relative expression. Each heatmap is accompanied by a
confidence bar indicating the membership confidence of genes in the
respective cluster.

## Usage

``` r
hc_plot_expression(AnnDatR, scaled_limits = c(0, 1), show_sample_labels = TRUE)
```

## Arguments

- AnnDatR:

  AnnDatR object (genes x samples)

- scaled_limits:

  List with zscore and scaled limits (default: zscore = c(-2,2), scaled
  = c(0,1))

- show_sample_labels:

  Logical, whether to show sample labels on x-axis (default: TRUE)

## Value

List with two named lists: \$zscore and \$scaled, each a list of ggplot
objects per cluster

## Examples

``` r
# Run clustering pipeline
adata_res <- hc_pca(example_adata, components = 40)
#> 111 of 76518 values (0.15%) are missing; using na_action = 'impute'.
adata_res <- hc_distance(adata_res, components = 20)
adata_res <- hc_snn(adata_res, neighbors = 15)
#> Building SNN graph from the provided distance matrix.
adata_res <- hc_cluster_consensus(adata_res, resolution = 8, n_seeds = 20)
#> Running louvain clustering at resolution 8 across 20 seeds.
#> Median number of communities across seeds: 34 (range 32-36).
#> Iteration: 0 *** value: 950.475
#> Iteration: 1 *** value: 273.082
#> Iteration: 2 *** value: 101.189
#> Iteration: 3 *** value: 92.0887
#> Iteration: 4 *** value: 90.7377
#> Iteration: 5 *** value: 86.5719
#> Iteration: 6 *** value: 86.5719
#> Minimum: 86.5719

# Plot expression heatmaps
expression_plots <- hc_plot_expression(adata_res, show_sample_labels = FALSE)
expression_plots$zscore[["1"]]  # Z-score plot for cluster 1

expression_plots$scaled[["1"]]  # Scaled plot for cluster 1
```
