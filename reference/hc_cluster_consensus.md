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
head(adata_res$uns$consensus_clustering)
#> # A tibble: 6 × 2
#>   gene            cluster
#>   <chr>           <chr>  
#> 1 ENSG00000000003 21     
#> 2 ENSG00000000005 5      
#> 3 ENSG00000000419 8      
#> 4 ENSG00000000457 6      
#> 5 ENSG00000000460 22     
#> 6 ENSG00000000938 35     
head(adata_res$obs)
#> # A tibble: 6 × 13
#>   ensembl_id      gene_symbol chromosome_name gene_biotype   start_position
#>   <chr>           <chr>       <chr>           <chr>                   <dbl>
#> 1 ENSG00000000003 TSPAN6      X               protein_coding      100627108
#> 2 ENSG00000000005 TNMD        X               protein_coding      100584936
#> 3 ENSG00000000419 DPM1        20              protein_coding       50934867
#> 4 ENSG00000000457 SCYL3       1               protein_coding      169849631
#> 5 ENSG00000000460 C1orf112    1               protein_coding      169662007
#> 6 ENSG00000000938 FGR         1               protein_coding       27612064
#> # ℹ 8 more variables: end_position <dbl>, length <dbl>, n_detected <dbl>,
#> #   sum_nTPM_counts <dbl>, detected <lgl>, variance <dbl>, cluster <chr>,
#> #   cluster_colors <chr>
```
