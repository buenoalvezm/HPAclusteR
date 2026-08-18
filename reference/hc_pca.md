# Perform PCA on an AnnDatR object

`hc_pca()` performs Principal Component Analysis (PCA) on the data
contained within an AnnDatR object.

## Usage

``` r
hc_pca(
  AnnDatR,
  components = 50,
  transform = c("log1p", "sqrt", "none"),
  scale_by = "sample",
  layer = NULL,
  na_action = c("impute", "nipals", "omit", "error"),
  verbose = TRUE
)
```

## Arguments

- AnnDatR:

  An AnnDatR object containing the data to perform PCA on.

- components:

  Number of principal components to compute (default is 50).

- transform:

  Data transformation method to apply before PCA. Options are `"log1p"`,
  `"sqrt"` or `"none"` (default is `"log1p"`).

- scale_by:

  Scaling method to apply to the data before PCA. Options include
  `"sample"`, `"gene"`, `"pareto"`, `"vector"`, `"min-max"`,
  `"max_sequencial"` or `"max"` (default is `"sample"`).

- layer:

  Optional layer name within the AnnDatR object to use for PCA. If
  `NULL`, the main data matrix is used (default is `NULL`).

- na_action:

  How to deal with missing values in the expression matrix. One of
  `"impute"` (default, replace missing values by their column mean and
  run a fast SVD), `"nipals"` (iterative PCA that handles missing values
  directly; exact but much slower, and requires the `pcaMethods`
  package), `"omit"` (drop genes with any missing value) or `"error"`
  (stop).

- verbose:

  Logical indicating whether to report on missing values (default is
  `TRUE`).

## Value

An AnnDatR object with the PCA fit stored in `uns$pca` (a list with
`scores`, `loadings`, `sdev`, `r2cum` and `method`) and the scores in
`obsm$X_pca`.

## Details

Missing values are the main performance trap in this step. Earlier
versions of this package delegated to
[`pcaMethods::pca()`](https://rdrr.io/pkg/pcaMethods/man/pca.html),
which silently fell back to the iterative NIPALS algorithm whenever it
encountered an `NA`. NIPALS is orders of magnitude slower than an SVD on
a complete matrix, and nothing in the output indicated why the call had
become slow.

`hc_pca()` now reports how many values are missing and, by default,
imputes them with their column mean so that the fast SVD path is taken.
Set `na_action = "nipals"` to reproduce the previous behaviour exactly.

## Examples

``` r
# Run PCA on an AnnDatR object
adata_res <- hc_pca(example_adata, components = 40)
#> 111 of 76518 values (0.15%) are missing; using na_action = 'impute'.

# PCA scores
head(adata_res$obsm$X_pca[, 1:5])
#>                       PC1       PC2        PC3        PC4        PC5
#> ENSG00000002745  2.525134  1.244945  4.9524275  2.9628522  0.5502585
#> ENSG00000004660 -3.978987  4.237604  1.0125002  2.3806353  0.1701373
#> ENSG00000006047 -3.335958 -5.230821  0.2785071 -0.8539639 -1.0674333
#> ENSG00000006059  0.384914 -1.662885  4.2574866  1.7776202  0.3784697
#> ENSG00000006453  4.327132 -2.237923  2.0932705 -0.3790134 -5.5485171
#> ENSG00000006740 -3.499337  4.827055 -1.0947394 -0.6904189 -3.5299211

# Cumulative explained variance
head(adata_res$uns$pca$r2cum)
#> [1] 0.2353205 0.4247893 0.5235696 0.6031322 0.6678779 0.7212278
```
