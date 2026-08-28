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
#> 22:15:16 UMAP embedding parameters a = 0.9922 b = 1.112
#> 22:15:16 Initializing from normalized Laplacian + noise (using RSpectra)
#> 22:15:16 Commencing optimization for 500 epochs, with 49428 positive edges
#> 22:15:16 Using rng type: pcg
#> 22:15:17 Optimization finished
adata_res <- hc_cluster_hulls(adata_res)
head(adata_res$uns$UMAP_hulls$hulls)
#> # A tibble: 6 × 7
#>   cluster sub_cluster sub_type landmass      X      Y polygon_id
#>   <chr>         <int> <chr>       <int>  <dbl>  <dbl> <chr>     
#> 1 32                1 primary         1 -0.929 -0.889 32_1_1    
#> 2 32                1 primary         1 -0.936 -0.882 32_1_1    
#> 3 32                1 primary         1 -0.949 -0.882 32_1_1    
#> 4 32                1 primary         1 -0.956 -0.889 32_1_1    
#> 5 32                1 primary         1 -0.956 -0.903 32_1_1    
#> 6 32                1 primary         1 -0.949 -0.909 32_1_1    
```
