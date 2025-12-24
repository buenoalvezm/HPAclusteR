# Compute distance matrix from PCA scores

`hc_distance()` computes a distance matrix based on PCA scores stored in
the AnnDatR object.

## Usage

``` r
hc_distance(AnnDatR, components = NULL, method = "spearman")
```

## Arguments

- AnnDatR:

  An AnnDatR object containing the data with PCA results.

- components:

  Number of principal components to be used. If NULL, it will be set to
  the first component with over 80% explained variance.

- method:

  Distance metric to use (default is 'spearman'). This must be one of
  "euclidean", "maximum", "manhattan", "canberra", "binary",
  "minkowski", "pearson", "spearman" or "kendall".

## Value

Distance matrix stored within the AnnDatR object.

## Examples

``` r
# Compute distance matrix on an AnnDatR object after PCA
adata_res <- hc_pca(example_adata, components = 40)
adata_res <- hc_distance(adata_res, components = 20)
head(adata_res$uns$distance)
#> [1] 0.8857143 1.3804511 1.0812030 0.9609023 1.2120301 1.2451128
```
