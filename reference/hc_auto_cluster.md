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
  snn_similarity = "euclidean",
  cluster_resolution = 6,
  cluster_method = "louvain",
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

- snn_similarity:

  Distance metric to use for nearest neighbor search (default is
  "euclidean"). Other options include "cosine", "manhattan", and
  "hamming".

- cluster_resolution:

  Resolution parameter for clustering (default is 6).

- cluster_method:

  Clustering method to use: "louvain" (default) or "leiden".

- verbose:

  Logical indicating whether to print progress messages (default is
  TRUE).

## Value

An AnnDataR object with clustering and intermediate results added. The
UMAP plot of the clusters is also printed.

## Examples

``` r
adata_res <- hc_auto_cluster(example_adata)
#> [1] "Kaiser's rule is above 80% variation. Success"
#> Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
#> Also defined by ‘spam’
#> Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
#> Also defined by ‘spam’
#> Building SNN based on a provided distance matrix
#> Computing SNN
#> Iteration: 0 *** value: 968.488
#> Iteration: 1 *** value: 236.854
#> Iteration: 2 *** value: 60.801
#> Iteration: 3 *** value: 56.147
#> Iteration: 4 *** value: 55.579
#> Iteration: 5 *** value: 55.5776
#> Iteration: 6 *** value: 55.5776
#> Minimum: 55.5776
#> Joining with `by = join_by(cons_cluster)`
#> Downloading uv...
#> Done!
#> The 'umap-learn' Python package is not installed. Installing it now...
#> Warning: An ephemeral virtual environment managed by 'reticulate' is currently in use.
#> To add more packages to your current session, call `py_require()` instead
#> of `py_install()`. Running:
#>   `py_require(c("umap-learn"))`
#> Warning: There was 1 warning in `dplyr::mutate()`.
#> ℹ In argument: `density_data = lapply(...)`.
#> ℹ In group 1: `cluster = "21"`, `sub_cluster = 1`, `sub_type = "primary"`.
#> Caused by warning:
#> ! The `x` argument of `as_tibble.matrix()` must have unique column names if
#>   `.name_repair` is omitted as of tibble 2.0.0.
#> ℹ Using compatibility `.name_repair`.
#> ℹ The deprecated feature was likely used in the HPAclusteR package.
#>   Please report the issue at
#>   <https://github.com/buenoalvezm/HPAclusteR/issues>.
#> Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> ℹ The deprecated feature was likely used in the HPAclusteR package.
#>   Please report the issue at
#>   <https://github.com/buenoalvezm/HPAclusteR/issues>.

head(adata_res$uns$consensus_clustering)
#> # A tibble: 6 × 2
#>   gene            cluster
#>   <chr>           <chr>  
#> 1 ENSG00000000003 23     
#> 2 ENSG00000000005 23     
#> 3 ENSG00000000419 3      
#> 4 ENSG00000000457 32     
#> 5 ENSG00000000460 22     
#> 6 ENSG00000000938 18     
```
