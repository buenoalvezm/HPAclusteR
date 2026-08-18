# Compare clusters using hypergeometric test

`hc_cluster_compare()` compares clusters from two AnnDatR objects using
a hypergeometric test to assess the significance of overlap between
clusters.

## Usage

``` r
hc_cluster_compare(
  AnnDatR_A,
  AnnDatR_B,
  graph_type = c("bipartite", "optimized")
)
```

## Arguments

- AnnDatR_A:

  First AnnDatR object

- AnnDatR_B:

  Second AnnDatR object

- graph_type:

  Type of graph to return: "bipartite" for simple overlap graph,
  "optimized" for optimized matching graph (default: "bipartite")

## Value

A list containing:

- matches: A tibble of cluster comparisons with overlap statistics and
  p-values.

- heatmap: A ggplot2 heatmap visualizing cluster overlaps.

- match_score (if `graph_type` is "optimized"): A numeric score
  representing the total normalized overlap percentage.

## Details

If `graph_type` is set to "optimized", the function will run a cluster
matching algorithm to ensure no many-to-many matches between clusters,
while still allowing one-to-many or many-to-one matches. The algorithm
will:

1.  Filter significant matches based on adjusted p-value threshold
    (0.05)

2.  Sort matches by overlap size in descending order

3.  Iteratively add matches to the result set, ensuring no many-to-many
    matches occur.

4.  Calculate a match score representing the total normalized overlap
    percentage.

## Examples

``` r
# Run clustering A
adata_res <- hc_pca(example_adata, components = 40)
#> 111 of 76518 values (0.15%) are missing; using na_action = 'impute'.
adata_res <- hc_distance(adata_res, components = 20)
adata_res <- hc_snn(adata_res, neighbors = 15)
#> Building SNN graph from the provided distance matrix.
adata_res_a <- hc_cluster_consensus(adata_res, resolution = 10, n_seeds = 20)
#> Running louvain clustering at resolution 10 across 20 seeds.
#> Median number of communities across seeds: 44 (range 42-47).
#> Iteration: 0 *** value: 955.881
#> Iteration: 1 *** value: 335.09
#> Iteration: 2 *** value: 130.357
#> Iteration: 3 *** value: 103.718
#> Iteration: 4 *** value: 103.171
#> Iteration: 5 *** value: 103.046
#> Iteration: 6 *** value: 102.761
#> Iteration: 7 *** value: 102.621
#> Iteration: 8 *** value: 102.621
#> Minimum: 102.621

# Run clustering B
adata_res_b <- hc_cluster_consensus(adata_res, resolution = 8, n_seeds = 20)
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

# Compare clusters between the two results
hc_cluster_compare(adata_res_a, adata_res_b)
#> 0.00% of genes do not match between the two datasets.
#> $matches
#> # A tibble: 47 × 8
#>    cluster_A cluster_B n_genes_A n_genes_B n_overlap percentage_overlap    p_val
#>    <chr>     <chr>         <int>     <int>     <int>              <dbl>    <dbl>
#>  1 1         16               34         9         4               44.4 1.35e- 4
#>  2 1         20               34        43        28               82.4 1.43e-37
#>  3 10        29               20        19        19              100   4.17e-39
#>  4 11        2                19        59        16               84.2 3.01e-18
#>  5 12        33               16        53        16              100   4.77e-22
#>  6 13        14               19        16        15               93.8 1.20e-28
#>  7 13        32               19        36         4               21.0 4.00e- 3
#>  8 14        31               28        33        26               92.9 1.48e-42
#>  9 15        7                36        40        36              100   1.30e-61
#> 10 16        1                24        29        24              100   1.55e-43
#> # ℹ 37 more rows
#> # ℹ 1 more variable: adj_p_val <dbl>
#> 
#> $heatmap

#> 
#> $network

#> 
```
