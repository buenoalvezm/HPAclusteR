# Determine number of components using Kaiser's rule

`hc_kaisers_rule()` determines the number of principal components to
retain based on Kaiser's rule, which suggests keeping components with
eigenvalues greater than 1. If the cumulative explained variance at this
point is less than 80%, it can optionally suggest an alternative number
of components that achieves at least 80% explained variance.

## Usage

``` r
hc_kaisers_rule(AnnDatR, with_alternative = TRUE)
```

## Arguments

- AnnDatR:

  An AnnDatR object containing the data with PCA results.

- with_alternative:

  Logical indicating whether to suggest an alternative number of
  components if the explained variance at Kaiser's rule is below 80%
  (default is TRUE).

## Value

Number of principal components to retain based on Kaiser's rule.

## Examples

``` r
# Determine number of components using Kaiser's rule
adata_res <- hc_pca(example_adata, components = 40)
n_components <- hc_kaisers_rule(adata_res)
#> [1] "Kaiser's rule is above 80% variation. Success"
print(n_components)
#> PC16 
#>   16 
```
