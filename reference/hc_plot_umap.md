# Plot UMAP with clusters and hulls

`hc_plot_umap()` creates a UMAP plot from an AnnData object, displaying
points, cluster hulls, or both.

## Usage

``` r
hc_plot_umap(AnnDatR, plot = c("points", "hulls", "both"))
```

## Arguments

- AnnDatR:

  An AnnData object containing UMAP coordinates and cluster hulls.

- plot:

  Type of plot to create: "points" for scatter plot, "hulls" for cluster
  hulls, or "both" for overlaying points on hulls. Default is "points".

## Value

A ggplot2 object representing the UMAP visualization.

## Examples

``` r
# Plot UMAP points from an AnnDatR object
adata_res <- hc_pca(example_adata, components = 40)
adata_res <- hc_distance(adata_res, components = 20)
adata_res <- hc_snn(adata_res, neighbors = 15)
#> Building SNN based on a provided distance matrix
#> Computing SNN
adata_res <- hc_cluster_consensus(adata_res, resolution = 6.3)
#> Iteration: 0 *** value: 970.898
#> Iteration: 1 *** value: 281.717
#> Iteration: 2 *** value: 100.383
#> Iteration: 3 *** value: 94.5865
#> Iteration: 4 *** value: 94.2944
#> Iteration: 5 *** value: 93.8117
#> Iteration: 6 *** value: 93.2751
#> Iteration: 7 *** value: 93.0317
#> Iteration: 8 *** value: 93.0109
#> Iteration: 9 *** value: 93.0109
#> Minimum: 93.0109
#> Joining with `by = join_by(cons_cluster)`
adata_res <- hc_umap(adata_res)
adata_res <- hc_cluster_hulls(adata_res)
hc_plot_umap(adata_res)
```
