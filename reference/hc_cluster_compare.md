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
adata_res <- hc_distance(adata_res, components = 20)
adata_res <- hc_snn(adata_res, neighbors = 15)
#> Building SNN based on a provided distance matrix
#> Computing SNN
adata_res_a <- hc_cluster_consensus(adata_res, resolution = 6)
#> Iteration: 0 *** value: 970.251
#> Iteration: 1 *** value: 311.124
#> Iteration: 2 *** value: 144.245
#> Iteration: 3 *** value: 110.026
#> Iteration: 4 *** value: 106.622
#> Iteration: 5 *** value: 106.162
#> Iteration: 6 *** value: 105.387
#> Iteration: 7 *** value: 103.152
#> Iteration: 8 *** value: 97.0298
#> Iteration: 9 *** value: 96.4117
#> Iteration: 10 *** value: 96.3528
#> Iteration: 11 *** value: 96.3275
#> Iteration: 12 *** value: 96.3135
#> Iteration: 13 *** value: 96.3135
#> Minimum: 96.3135
#> Joining with `by = join_by(gene)`
#> Joining with `by = join_by(new_cluster)`
#> Joining with `by = join_by(cons_cluster)`

# Run clustering B
adata_res_b <- hc_cluster_consensus(adata_res, resolution = 7)
#> Iteration: 0 *** value: 973.302
#> Iteration: 1 *** value: 290.071
#> Iteration: 2 *** value: 127.083
#> Iteration: 3 *** value: 106.654
#> Iteration: 4 *** value: 105.703
#> Iteration: 5 *** value: 104.763
#> Iteration: 6 *** value: 104.583
#> Iteration: 7 *** value: 104.554
#> Iteration: 8 *** value: 104.549
#> Iteration: 9 *** value: 104.545
#> Iteration: 10 *** value: 104.545
#> Minimum: 104.545
#> Joining with `by = join_by(gene)`
#> Joining with `by = join_by(new_cluster)`
#> Joining with `by = join_by(cons_cluster)`

# Compare clusters between the two results
hc_cluster_compare(adata_res_a, adata_res_b)
#> 0.00% of genes do not match between the two datasets.
#> $matches
#> # A tibble: 43 × 8
#>    cluster_A cluster_B n_genes_A n_genes_B n_overlap percentage_overlap    p_val
#>    <chr>     <chr>         <int>     <int>     <int>              <dbl>    <dbl>
#>  1 1         1                29        31        29              100   6.19e-54
#>  2 10        15               37         7         2               28.6 2.49e- 2
#>  3 10        9                37        35        35              100   1.26e-62
#>  4 11        10               53        41        41              100   2.05e-62
#>  5 11        11               53         6         4               66.7 9.79e- 5
#>  6 11        42               53         7         7              100   7.93e-10
#>  7 12        41               19        20        19              100   2.89e-39
#>  8 13        12               28        27        27              100   4.34e-52
#>  9 14        13               22        22        22              100   1.42e-45
#> 10 15        14               34        32        32              100   2.44e-58
#> # ℹ 33 more rows
#> # ℹ 1 more variable: adj_p_val <dbl>
#> 
#> $heatmap

#> 
#> $network

#> 
```
