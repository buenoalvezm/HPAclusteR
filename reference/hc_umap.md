# Create UMAP embeddings from SNN graph

`hc_umap()` generates UMAP embeddings based on the Shared Nearest
Neighbors (SNN) graph stored in the AnnDatR object.

## Usage

``` r
hc_umap(
  AnnDatR,
  n_epochs = NULL,
  min_dist = 0.3,
  spread = 1,
  seed = 42,
  verbose = TRUE
)
```

## Arguments

- AnnDatR:

  An AnnDatR object containing the data with SNN graph results.

- n_epochs:

  Number of epochs for UMAP optimization. If `NULL`, it will be set to
  200 for datasets with more than 10,000 observations and 500 otherwise
  (default is `NULL`).

- min_dist:

  Minimum distance between embedded points (default is 0.3).

- spread:

  Effective scale of the embedded points (default is 1).

- seed:

  Random seed for reproducibility (default is 42).

- verbose:

  Logical indicating whether to print progress messages (default is
  `TRUE`).

## Value

An AnnDatR object with the raw embedding in `obsm$X_umap_raw`, the
rescaled embedding in `obsm$X_umap`, and `UMAP1`/`UMAP2` columns added
to `obs`.

## Details

The embedding is produced by
[`uwot::optimize_graph_layout()`](https://jlmelville.github.io/uwot/reference/optimize_graph_layout.html),
a native R and C++ implementation of UMAP's simplicial set embedding.
Earlier versions called the Python `umap-learn` package through
`reticulate` by way of `Seurat::RunUMAP()`, which required a working
Python installation with pinned `numpy` and `umap-learn` versions and
was the main source of installation and continuous integration failures.

The algorithm and its parameters are the same, but the two
implementations use different random number generators and optimisers,
so embeddings are equivalent in structure rather than numerically
identical to those produced by earlier versions.

## Examples

``` r
# Generate UMAP embeddings on an AnnDatR object after SNN graph computation
adata_res <- hc_pca(example_adata, components = 40)
#> 111 of 76518 values (0.15%) are missing; using na_action = 'impute'.
adata_res <- hc_distance(adata_res, components = 20)
adata_res <- hc_snn(adata_res, neighbors = 15)
#> Building SNN graph from the provided distance matrix.
adata_res <- hc_umap(adata_res)
#> 22:15:33 UMAP embedding parameters a = 0.9922 b = 1.112
#> 22:15:33 Initializing from normalized Laplacian + noise (using RSpectra)
#> 22:15:33 Commencing optimization for 500 epochs, with 49428 positive edges
#> 22:15:33 Using rng type: pcg
#> 22:15:35 Optimization finished
head(adata_res$obsm$X_umap)
#>                      UMAP_1     UMAP_2
#> ENSG00000002745  0.79382229  0.9476842
#> ENSG00000004660  0.10275772  0.3765000
#> ENSG00000006047 -0.83528043 -0.8373416
#> ENSG00000006059  0.67569219  0.4034760
#> ENSG00000006453  0.57238511  0.3686830
#> ENSG00000006740  0.05058434  0.3516931
```
