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
adata_res <- hc_umap(adata_res)
#> 22:15:28 UMAP embedding parameters a = 0.9922 b = 1.112
#> 22:15:28 Initializing from normalized Laplacian + noise (using RSpectra)
#> 22:15:28 Commencing optimization for 500 epochs, with 49428 positive edges
#> 22:15:28 Using rng type: pcg
#> 22:15:30 Optimization finished
adata_res <- hc_cluster_hulls(adata_res)
hc_plot_umap(adata_res)
```
