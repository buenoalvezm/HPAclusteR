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
adata_res <- hc_distance(adata_res, components = 20)
adata_res <- hc_snn(adata_res, neighbors = 15)
#> Building SNN based on a provided distance matrix
#> Computing SNN
adata_res <- hc_cluster_consensus(adata_res, resolution = 7)
#> Iteration: 0 *** value: 948.279
#> Iteration: 1 *** value: 70.5063
#> Iteration: 2 *** value: 22.246
#> Iteration: 3 *** value: 22.2447
#> Iteration: 4 *** value: 22.2447
#> Minimum: 22.2447
#> Joining with `by = join_by(cons_cluster)`

# Plot expression heatmaps
expression_plots <- hc_plot_expression(adata_res, show_sample_labels = FALSE)
expression_plots$zscore[["1"]]  # Z-score plot for cluster 1

expression_plots$scaled[["1"]]  # Scaled plot for cluster 1
```
