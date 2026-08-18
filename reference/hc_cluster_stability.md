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
hc_cluster_stability(adata_res)
#> $cluster_stability_hist

#> 
#> $ari_scores
#> # A tibble: 190 × 3
#>    seed1 seed2   ARI
#>    <int> <int> <dbl>
#>  1     1     2 0.856
#>  2     1     3 0.869
#>  3     1     4 0.861
#>  4     1     5 0.901
#>  5     1     6 0.875
#>  6     1     7 0.917
#>  7     1     8 0.816
#>  8     1     9 0.893
#>  9     1    10 0.899
#> 10     1    11 0.872
#> # ℹ 180 more rows
#> 
```
