# Run GO enrichment analysis (BP, MF, CC) for all clusters using clusterProfiler::enrichGO

Run GO enrichment analysis (BP, MF, CC) for all clusters using
clusterProfiler::enrichGO

## Usage

``` r
run_go_enrichment(AnnDatR, universe = NULL, verbose = TRUE)
```

## Arguments

- AnnDatR:

  AnnDatR object containing clustering results

- universe:

  Character vector of background genes (default: NULL, all genes in
  clustering data)

- verbose:

  Logical indicating whether to print progress messages (default is
  TRUE).

## Value

Tibble with enrichment results for all clusters and all GO ontologies,
with Ensembl gene IDs
