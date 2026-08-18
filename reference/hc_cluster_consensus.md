# Consensus clustering of genes

Clusters genes using multiple runs of a clustering algorithm with
different seeds and resolutions, then derives a consensus clustering.

## Usage

``` r
hc_cluster_consensus(
  AnnDatR,
  resolution = 6,
  method = c("louvain", "leiden"),
  n_seeds = 100,
  seed = 42,
  min_k = 30,
  max_k = 110,
  verbose = TRUE
)
```

## Arguments

- AnnDatR:

  An AnnDatR object containing the data with SNN graph.

- resolution:

  Resolution parameter for clustering (default is 6).

- method:

  Clustering method to use: "louvain" (default) or "leiden".

- n_seeds:

  Number of different random seeds to use for clustering (default is
  100).

- seed:

  Random seed for reproducibility (default is 42).

- min_k:

  Lower bound on the median number of communities across seeds (default
  is 30). Resolutions that produce fewer are rejected.

- max_k:

  Upper bound on the median number of communities across seeds (default
  is 110).

- verbose:

  Logical indicating whether to print progress messages (default is
  TRUE).

## Value

Consensus clustering results stored within the AnnDatR object.

## Details

The resolution is validated against `min_k` and `max_k`: if the median
number of communities across seeds falls outside that window, the
function stops and asks for a different `resolution`. These bounds were
previously hard-coded.

Re-running this function on an object that already carries a clustering
replaces the previous `cluster` and `cluster_colors` columns in `obs`.

## Examples

``` r
# Perform consensus clustering on an AnnDatR object after SNN graph computation
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
head(adata_res$uns$consensus_clustering)
#> # A tibble: 6 × 2
#>   gene            cluster
#>   <chr>           <chr>  
#> 1 ENSG00000002745 17     
#> 2 ENSG00000004660 22     
#> 3 ENSG00000006047 4      
#> 4 ENSG00000006059 16     
#> 5 ENSG00000006453 8      
#> 6 ENSG00000006740 19     
head(adata_res$obs)
#> # A tibble: 6 × 9
#>   ensembl_id      gene_symbol chromosome_name gene_biotype   start_position
#>   <chr>           <chr>       <chr>           <chr>                   <dbl>
#> 1 ENSG00000002745 WNT16       7               protein_coding      121325367
#> 2 ENSG00000004660 CAMKK1      17              protein_coding        3860315
#> 3 ENSG00000006047 YBX2        17              protein_coding        7288263
#> 4 ENSG00000006059 KRT33A      17              protein_coding       41346092
#> 5 ENSG00000006453 BAIAP2L1    7               protein_coding       98291650
#> 6 ENSG00000006740 ARHGAP44    17              protein_coding       12789498
#> # ℹ 4 more variables: end_position <dbl>, length <dbl>, cluster <chr>,
#> #   cluster_colors <chr>
```
