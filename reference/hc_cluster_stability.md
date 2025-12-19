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
adata_res <- hc_cluster_consensus(adata_res, resolution = 6.3)
#> Iteration: 0 *** value: 970.898
#> Iteration: 1 *** value: 281.717
#> Iteration: 2 *** value: 100.383
#> Iteration: 3 *** value: 94.5865
#> Iteration: 4 *** value: 94.2944
#> Iteration: 5 *** value: 93.8117
#> Iteration: 6 *** value: 93.2751
#> Iteration: 7 *** value: 93.0317
#> Iteration: 8 *** value: 93.0109
#> Iteration: 9 *** value: 93.0109
#> Minimum: 93.0109
#> Joining with `by = join_by(cons_cluster)`
hc_cluster_stability(adata_res)
#> $cluster_stability_hist

#> 
#> $ari_scores
#> # A tibble: 4,950 × 3
#> # Groups:   seed1, seed2 [4,950]
#>    seed1 seed2   ARI
#>    <int> <int> <dbl>
#>  1     1     2 0.904
#>  2     1     3 0.918
#>  3     1     4 0.903
#>  4     1     5 0.932
#>  5     1     6 0.841
#>  6     1     7 0.882
#>  7     1     8 0.903
#>  8     1     9 0.884
#>  9     1    10 0.879
#> 10     1    11 0.865
#> # ℹ 4,940 more rows
#> 
```
