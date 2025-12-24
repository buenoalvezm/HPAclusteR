# Consensus clustering of genes

Clusters genes using multiple runs of a clustering algorithm with
different seeds and resolutions, then derives a consensus clustering.

## Usage

``` r
hc_cluster_consensus(
  AnnDatR,
  resolution = 6,
  method = "louvain",
  n_seeds = 100,
  seed = 42,
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

- verbose:

  Logical indicating whether to print progress messages (default is
  TRUE).

## Value

Consensus clustering results stored within the AnnDatR object.

## Examples

``` r
# Perform consensus clustering on an AnnDatR object after SNN graph computation
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
head(adata_res$uns$consensus_clustering)
#> # A tibble: 6 × 2
#>   gene            cluster
#>   <chr>           <chr>  
#> 1 ENSG00000002745 13     
#> 2 ENSG00000004660 22     
#> 3 ENSG00000006047 12     
#> 4 ENSG00000006059 19     
#> 5 ENSG00000006453 19     
#> 6 ENSG00000006740 22     
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
