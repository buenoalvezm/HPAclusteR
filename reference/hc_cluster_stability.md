# Calculate ARI scores for multiple clustering seeds

`hc_cluster_stability()` computes the Adjusted Rand Index (ARI) scores
between clustering results obtained from multiple random seeds.

## Usage

``` r
hc_cluster_stability(AnnDatR)
```

## Arguments

- AnnDatR:

  AnnDatR An AnnDatR object containing clustering results.

## Value

A tibble containing pairs of seeds and their corresponding ARI scores.

## Examples

``` r
# Perform clustering with multiple seeds and calculate ARI scores
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
hc_cluster_stability(adata_res)
#> $cluster_stability_hist

#> 
#> $ari_scores
#> # A tibble: 4,950 × 3
#> # Groups:   seed1, seed2 [4,950]
#>    seed1 seed2   ARI
#>    <int> <int> <dbl>
#>  1     1     2 0.979
#>  2     1     3 0.962
#>  3     1     4 0.990
#>  4     1     5 0.949
#>  5     1     6 0.954
#>  6     1     7 0.963
#>  7     1     8 0.953
#>  8     1     9 0.974
#>  9     1    10 0.939
#> 10     1    11 0.979
#> # ℹ 4,940 more rows
#> 
```
