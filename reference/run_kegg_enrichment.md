# Run KEGG enrichment analysis (with Ensembl-to-Entrez conversion)

Run KEGG enrichment analysis (with Ensembl-to-Entrez conversion)

## Usage

``` r
run_kegg_enrichment(AnnDatR, universe = NULL, verbose = TRUE)
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

Tibble with enrichment results for all clusters
