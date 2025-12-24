# Calculate UMAP cluster hulls

`hc_cluster_hulls()` computes density-based cluster hulls on UMAP
embeddings stored within an AnnDatR object. It identifies subclusters
using DBSCAN, estimates density, filters based on cumulative density and
fraction limits, and generates concave hulls for each cluster.

## Usage

``` r
hc_cluster_hulls(
  AnnDatR,
  n = 300,
  cum_z_lim = 0.95,
  frac_lim = 0.05,
  plot_range = rep(c(min(c(AnnDatR[["obs"]][["UMAP1"]], AnnDatR[["obs"]][["UMAP2"]])),
    max(c(AnnDatR[["obs"]][["UMAP1"]], AnnDatR[["obs"]][["UMAP2"]]))), 2) * 1.05,
  poly_concavity = 1,
  poly_smoothing = 1,
  relative_bandwidth = 1/200
)
```

## Arguments

- AnnDatR:

  An AnnDatR object containing UMAP embeddings and cluster assignments.

- n:

  Number of grid points for density estimation (default is 300).

- cum_z_lim:

  Cumulative density limit for filtering (default is 0.95).

- frac_lim:

  Minimum fraction limit for subcluster inclusion (default is 0.05).

- plot_range:

  Range of the plot in the format c(xmin, xmax, ymin, ymax) (default is
  1.05 times the min and max of UMAP1 and UMAP2).

- poly_concavity:

  Controls the concavity of the polygons (default is 1). Higher values
  result in less detailed polygons.

- poly_smoothing:

  Controls the smoothing of the polygons (default is 1). Higher values
  result in less detailed polygons.

- relative_bandwidth:

  Relative bandwidth for density estimation as a fraction of the plot
  diagonal (default is 1/200).

## Value

AnnDatR object with UMAP cluster hulls within the AnnDatR object.

## Examples

``` r
# Calculate UMAP cluster hulls on an AnnDatR object after UMAP and clustering
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
adata_res <- hc_umap(adata_res)
adata_res <- hc_cluster_hulls(adata_res)
head(adata_res$uns$UMAP_hulls$hulls)
#> # A tibble: 6 × 7
#>   cluster sub_cluster sub_type landmass      X     Y polygon_id
#>   <chr>         <dbl> <chr>       <dbl>  <dbl> <dbl> <chr>     
#> 1 7                 1 primary         1 -0.565 0.432 7_1_1     
#> 2 7                 1 primary         1 -0.572 0.439 7_1_1     
#> 3 7                 1 primary         1 -0.579 0.446 7_1_1     
#> 4 7                 1 primary         1 -0.572 0.453 7_1_1     
#> 5 7                 1 primary         1 -0.565 0.46  7_1_1     
#> 6 7                 1 primary         1 -0.558 0.467 7_1_1     
```
