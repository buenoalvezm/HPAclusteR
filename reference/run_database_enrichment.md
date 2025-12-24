# Run enrichment analysis for all annotation databases

Run enrichment analysis for all annotation databases

## Usage

``` r
run_database_enrichment(
  annotation_terms,
  AnnDatR,
  universe = NULL,
  verbose = TRUE
)
```

## Arguments

- annotation_terms:

  Tibble from build_annotation_terms_tibble() (db_id, ensg_id, term,
  term_id)

- AnnDatR:

  AnnDatR object containing clustering results

- universe:

  Character vector of background genes (default: NULL, all genes in
  clustering data)

- verbose:

  Logical indicating whether to print progress messages (default is
  TRUE).

## Value

Tibble with enrichment results for all clusters and databases
