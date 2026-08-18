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

  Number of principal components to be used. If `NULL`, it will be set
  to the first component with over 80% explained variance.

- method:

  Distance metric to use (default is `"spearman"`). This must be one of
  `"euclidean"`, `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"`,
  `"minkowski"`, `"pearson"`, `"spearman"` or `"kendall"`.

## Value

An AnnDatR object with the distance matrix stored in `uns$distance`.

## Examples

``` r
# Compute distance matrix on an AnnDatR object after PCA
adata_res <- hc_pca(example_adata, components = 40)
#> 111 of 76518 values (0.15%) are missing; using na_action = 'impute'.
adata_res <- hc_distance(adata_res, components = 20)
head(adata_res$uns$distance)
#> [1] 1.0285714 1.2345865 0.8436090 0.8045113 1.2827068 1.1759398
```
