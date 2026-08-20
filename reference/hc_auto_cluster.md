# One-call clustering pipeline

`hc_auto_cluster` performs a full clustering pipeline on an AnnDataR
object, including PCA, distance calculation, SNN graph construction,
consensus clustering, UMAP embedding, and cluster hull generation. It
prints a UMAP plot of the resulting clusters.

## Usage

``` r
hc_auto_cluster(
  AnnDatR,
  layer = NULL,
  pca_components = 50,
  pca_transform = "log1p",
  distance_components = NULL,
  distance_method = "spearman",
  snn_neighbors = 20,
  snn_prune = 1/15,
  cluster_resolution = 6,
  cluster_method = "louvain",
  cluster_seeds = 100,
  verbose = TRUE
)
```

## Arguments

- AnnDatR:

  An AnnDataR object.

- layer:

  Optional layer name within the AnnDatR object to use instead of the
  default data.

- pca_components:

  Number of principal components to compute (default is 50).

- pca_transform:

  Data transformation method to apply before PCA. Options are "log1p",
  "sqrt", or "none" (default is "log1p").

- distance_components:

  Number of principal components to be used. If NULL, uses Kaiser’s rule
  to determine the number of components.

- distance_method:

  Distance metric to use (default is 'spearman'). This must be one of
  "euclidean", "maximum", "manhattan", "canberra", "binary",
  "minkowski", "pearson", "spearman" or "kendall".

- snn_neighbors:

  Number of nearest neighbors to consider (default is 20).

- snn_prune:

  Pruning threshold for the SNN graph (default is 1/15).

- cluster_resolution:

  Resolution parameter for clustering (default is 6).

- cluster_method:

  Clustering method to use: "louvain" (default) or "leiden".

- cluster_seeds:

  Number of random seeds to aggregate over in the consensus step
  (default is 100).

- verbose:

  Logical indicating whether to print progress messages (default is
  TRUE).

## Value

An AnnDataR object with clustering and intermediate results added. The
UMAP plot of the clusters is also printed.

## Examples

``` r
adata_res <- hc_auto_cluster(
  example_adata,
  cluster_resolution = 8,
  cluster_seeds = 20
)
#> 111 of 76518 values (0.15%) are missing; using na_action = 'impute'.
#> Kaiser's rule is above 80% variation. Success.
#> Building SNN graph from the provided distance matrix.
#> Running louvain clustering at resolution 8 across 20 seeds.
#> Median number of communities across seeds: 31 (range 29-33).
#> Iteration: 0 *** value: 948.105
#> Iteration: 1 *** value: 226.66
#> Iteration: 2 *** value: 70.0129
#> Iteration: 3 *** value: 69.8398
#> Iteration: 4 *** value: 69.8181
#> Iteration: 5 *** value: 69.8181
#> Minimum: 69.8181
#> 11:39:46 UMAP embedding parameters a = 0.9922 b = 1.112
#> 11:39:47 Initializing from normalized Laplacian + noise (using RSpectra)
#> 11:39:47 Commencing optimization for 500 epochs, with 55050 positive edges
#> 11:39:47 Using rng type: pcg
#> 11:39:48 Optimization finished

head(adata_res$uns$consensus_clustering)
#> # A tibble: 6 × 2
#>   gene            cluster
#>   <chr>           <chr>  
#> 1 ENSG00000002745 23     
#> 2 ENSG00000004660 30     
#> 3 ENSG00000006047 4      
#> 4 ENSG00000006059 3      
#> 5 ENSG00000006453 11     
#> 6 ENSG00000006740 8      
```
